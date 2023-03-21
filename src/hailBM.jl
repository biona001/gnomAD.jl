"""
Currently a `HailBlockMatrix` is just a python wrapper, where data reading is 
achieved within python by the Hail software, then the result is passed into
Julia via PyCall.jl. The `bm` field behaves as a `BlockMatrix` from hail. 
In the future, a native Julia reader that directly reads
the compreseed format is desired. 

The Hail Block Matrix format is described here
https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#blockmatrix

A more specific file description is provided here
https://discuss.hail.is/t/blockmatrix-specification/3118
"""
struct HailBlockMatrix <: AbstractMatrix{Float64}
    bm_files::String
    chr::Vector{String}
    pos::Vector{Int}
    ref::Vector{String}
    alt::Vector{String}
    af::Vector{Float64}
    build::String
    bm::PyObject
end

# do not allow getindex for a single dimension
Base.getindex(x::HailBlockMatrix, i::Int) = 
    error("Can only access HailBlockMatrix using 2 indices (i, j)")
Base.getindex(x::HailBlockMatrix, ::AbstractRange) = 
    error("Can only access HailBlockMatrix using 2 dimension index, e.g. bm[1:10, 1:10]")

function Base.getindex(x::HailBlockMatrix, i::Int, j::Int)
    i ≤ 0 || j ≤ 0 && error("attempt to access matrix at index [$i, $j]")
    i > j ? 0.0 : x.bm[i, j]
end
function Base.getindex(x::HailBlockMatrix, idx1::UnitRange, idx2::UnitRange)
    _get_block(x, first(idx1), last(idx1), first(idx2), last(idx2))
end
function _get_block(x::HailBlockMatrix, 
    x_start::Int, x_end::Int, y_start::Int, y_end::Int)
    1 ≤ x_start ≤ size(x, 1) && 1 ≤ y_start ≤ size(x, 1) && 
        1 ≤ x_end ≤ size(x, 1) && 1 ≤ y_end ≤ size(x, 1) || 
        error("x_start out of bounds")
    # convert function uses python indexing, i.e. 0 based
    return py"convert"(x.bm, x_start - 1, x_end, y_start - 1, y_end)
end

"""
    get_block(bm::HailBlockMatrix, chr, start_bp, end_bp)

Reads in a block of LD matrix from chromosome `chr` between basepairs 
`start_bp` and `end_bp`. The inputs `chr`, `start_bp`, `end_bp` can be Int
or String. The result will always be a matrix even if `start_bp == end_bp`.

If `start_bp` or `end_bp` is not in the LD matrix, we will return
the smallest region that does NOT include them. For example, if 
`(start_bp, end_bp) = (555, 777)`
and SNP positions in the LD panel are 
`bm.pos    = [..., 400,  500,  600,  700,  800,  900, ...]`
`snp_names = [..., snp4, snp5, snp6, snp7, snp8, snp9, ..]`
then we will return the LD matrix for `[snp6, snp7]`

# Note
Make sure `start_bp` and `end_bp` is from the same human genome build as the 
LD matrices. Pan-UKBB and genomAD both use hg38.
"""
function get_block(bm::HailBlockMatrix, chr::Union{String, Int}, 
    start_bp::Union{String, Int}, end_bp::Union{String, Int})
    # convert start_bp/end_bp to Int and chr to String
    start_bp = typeof(start_bp) == Int ? start_bp : parse(Int, start_bp)
    end_bp = typeof(end_bp) == Int ? end_bp : parse(Int, end_bp)
    chr = typeof(chr) == String ? chr : string(chr)
    # search for starting/ending position
    chr_range = findall(x -> x == chr, bm.chr)
    chr_pos = @view(bm.pos[chr_range])
    idx_start = searchsortedfirst(chr_pos, start_bp)
    idx_end = searchsortedlast(chr_pos, end_bp)
    # check for errors
    idx_end < chr_pos[1] && error("end_bp occurs before first SNP in bm")
    idx_start > chr_pos[end] && error("start_bp occurs after last SNP in bm")
    return bm[idx_start:idx_end, idx_start:idx_end]
end

Base.size(x::HailBlockMatrix, k::Int) = 
    k == 1 ? x.bm.n_rows : k == 2 ? x.bm.n_cols : k > 2 ? 1 : error("Dimension k out of range")
Base.size(x::HailBlockMatrix) = (size(x, 1), size(x, 2))

"""
    hail_block_matrix(bm_files::String, ht_files::String)

Creates a `HailBlockMatrix` which allows reading chunks of LD matrix data into 
memory. It is a subtype of `AbstractMatrix`, so operations like indexing and 
`size()` works, but only accessing its elements in contiguous chunks is "fast". 
I may support more functionalities depending on interest. 

# Examples

```julia
using EasyLD
datadir = "/Users/biona001/.julia/dev/EasyLD/data"
bm_file = joinpath(datadir, "UKBB.EUR.ldadj.bm")
ht_file = joinpath(datadir, "UKBB.EUR.ldadj.variant.ht")
bm = hail_block_matrix(bm_file, ht_file); # need a ';' to avoid displaying a few entries of bm, which takes ~0.1 seconds per entry

# get matrix dimension
size(bm) # returns (23960350, 23960350)

# read a single entry
bm[1, 1] # returns 0.999959841549239

# read first 10k by 10k block into memory (takes roughly 7 seconds)
bm[1:10000, 1:10000]

# arbitrary slicing works but is very slow
bm[1:3, 1:2:100] # ~22 seconds

# read a specific chromosome region
chr = 1
start_pos = 11063
end_pos = 91588
sigma = get_block(bm, chr, start_pos, end_pos)
```
"""
function hail_block_matrix(bm_files::String, ht_files::String)
    isdir(bm_files) || error("Directory $bm_files does not exist")
    chr, pos, ref, alt, af, build = get_row_info(ht_files)
    return HailBlockMatrix(bm_files, chr, pos, ref, alt, af, build,
        hail_linalg.BlockMatrix.read(bm_files))
end

"""
    read_variant_index_tables(ht_file::String)

Read variant index hail tables into a DataFrame. The first time this function
gets called, we will read the original `.ht` files into memory and write the
result to a tab-separated value `.tsv` file into the same directory as `ht_file`
"""
function read_variant_index_tables(ht_file::String)
    isdir(ht_file) || error("$ht_file is not a directory")
    tsv_file = joinpath(ht_file, "variant.ht.tsv")
    if !isfile(tsv_file)
        println("Exporting raw hail table data to file $tsv_file"); flush(stdout)
        hail_table = hail.read_table(ht_file)
        hail_table.export(tsv_file)
    end
    return CSV.read(tsv_file, DataFrame)
end

"""
    get_row_info(ht_file::String)

Reads chromosome/position/reference-allele/alternate-allele/alternate-allele-frequency
of each variable in the variant index files, and the human genome build. 
Some panels (e.g. Pan-UKB) also contain rsID, but others (e.g. gnomAD) does not,
so rsID is not saved by default. 
"""
function get_row_info(ht_file::String)
    df = read_variant_index_tables(ht_file)
    locus = split.(df[!, "locus"], ':')
    chr = [locus[i][1] for i in eachindex(locus)] |> Vector{String}
    pos = [parse(Int, locus[i][2]) for i in eachindex(locus)]
    ref, alt = _extract_ref_alt_alleles(df[!, "alleles"])
    build = extract_genome_build(ht_file)
    af = "AF" in names(df) ? df[!, "AF"] : extract_alternate_allele_freq(ht_file)
    # check that within each chr, basepair positions are sorted
    for c in unique(chr)
        idx = findall(x -> x == c, chr)
        issorted(@view(pos[idx])) || 
            error("Basepair positions in chr $chr not sorted!")
    end
    return chr, pos, ref, alt, af, build
end

"""
    extract_alternate_allele_freq(ht_file::String)

Extract alternate allele frequency from the `ht_file` assuming this information 
is present in the variant index files as JSON format with header `header_name`

For example, in gnomAD, we have
```julia
julia> ht_file = "gnomad.genomes.r2.1.1.nfe_test.common.adj.ld.variant_indices.ht"
julia> df = read_variant_index_tables(ht_file)
julia> names(df)
4-element Vector{String}:
 "locus"
 "alleles"
 "pop_freq"
 "idx"
julia> df[1, "pop_freq"] # first entry with header `pop_freq` that contains AF
"{\"AC\":861,\"AF\":0.40922053231939165,\"AN\":2104,\"homozygote_count\":166}"
```
Then the AF field for all variants will be extracted
"""
function extract_alternate_allele_freq(ht_file::String; header_name = "pop_freq")
    df = read_variant_index_tables(ht_file)
    af = zeros(Float64, size(df, 1))
    for (i, s) in enumerate(df[!, header_name])
        d = JSON.parse(s)
        af[i] = d["AF"]
    end
    return af
end

"""
    _extract_ref_alt_alleles(x::Vector{String})

Helper function to extract the ref/alt alleles from hail tables. Users should call
`get_row_info` instead of this function. 
"""
function _extract_ref_alt_alleles(alleles::Vector{String})
    ref, alt = String[], String[]
    for allele in alleles
        # idx = findall(x -> x == '"', allele)
        idx1 = findfirst(x -> x == '"', allele)
        idx2 = findnext(x -> x == '"', allele, idx1+1)
        idx3 = findnext(x -> x == '"', allele, idx2+1)
        idx4 = findlast(x -> x == '"', allele)
        push!(ref, allele[idx1+1:idx2-1])
        push!(alt, allele[idx3+1:idx4-1])
    end
    # check ref/alt is different
    for (i, (r, a)) in enumerate(zip(ref, alt))
        r == a && error("Ref allele is equal to alt allele at position $i")
    end
    return ref, alt
end
# function _extract_ref_alt_alleles(alleles::Vector{String}) # too slow
#     regex = r"\"(.*?)\""
#     refs, alts = String[], String[]
#     for allele in alleles
#         A1, A2 = split(s, ',')
#         push!(refs, match(regex, A1).captures[1])
#         push!(alts, match(regex, A2).captures[1])
#     end
#     return refs, alts
# end

"""
    extract_human_genome_build(ht_file::String)

Extracts the human genome build information from the variant index folder. We 
assume this information is stored in a file called `metadata.json` or 
`metadata.json.gz` for format `Locus(XXX)` where `XXX` is the reference build.
If this format is not detected, the genome build information will be "Unknown"
"""
function extract_genome_build(ht_file::String)
    locus_dat = ""
    if isfile(joinpath(ht_file, "metadata.json"))
        locus_dat = read_metadata(joinpath(ht_file, "metadata.json"))["table_type"]
    elseif isfile(joinpath(ht_file, "metadata.json.gz"))
        locus_dat = read_metadata(joinpath(ht_file, "metadata.json.gz"))["table_type"]      
    end
    # search for pattern Locus(XXX)
    regex = r"Locus\((.*?)\)"
    return occursin(regex, locus_dat) ? match(regex, locus_dat).captures[1] : "Unknown"
end

"""
    read_metadata(metafile::String)

Reads the metadata file and returns a dictionary. `metafile` should end with 
`.json` or `.json.gz`
"""
function read_metadata(metafile::String)
    if endswith(metafile, ".json")
        s = read(metafile) |> String
    elseif endswith(metafile, ".json.gz")
        open(GzipDecompressorStream, metafile) do io
            s = read(io) |> String
        end
    else
        error("metafile should end with .json or .json.gz")
    end
    return JSON.parse(s)
end
