"""
Currently a `HailBlockMatrix` is just a python wrapper, where data reading is 
achieved within python by the Hail software, then the result is passed into
Julia via PyCall.jl. The `bm` field behaves as a `BlockMatrix` from hail. The
`info` fields stores the `chr/pos/ref/alt/AF` information, where `pos` uses
coordinates in genome build `build`.

The Hail Block Matrix format is described here
https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#blockmatrix

A more specific file description is provided here
https://discuss.hail.is/t/blockmatrix-specification/3118
"""
struct HailBlockMatrix <: AbstractMatrix{Float64}
    bm_files::String
    info::DataFrame
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
    get_block(bm::HailBlockMatrix, chr, start_bp, end_bp; [min_maf], 
        [snps_to_keep], [min_eigval], [enforce_psd], [rk])

Reads in a block of LD matrix from chromosome `chr` between basepairs 
`start_bp` and `end_bp`. The inputs `chr`, `start_bp`, `end_bp` can be Int
or String. The result will always be a matrix even if `start_bp == end_bp`. If 
`start_bp` or `end_bp` is not in the LD matrix, we will return
the smallest region that does NOT include them. For example, if 
`(start_bp, end_bp) = (555, 777)`
and SNP positions in the LD panel are 
`positions = [..., 400,  500,  600,  700,  800,  900, ...]`
with
`snp_names = [..., snp4, snp5, snp6, snp7, snp8, snp9, ..]`, 
then we will return the LD matrix for `[snp6, snp7]`

# Inputs
+ `bm`: A `HailBlockMatrix`
+ `chr`: Chromosome number, can be an Int or String (e.g. `1` or `"1"`)
+ `start_bp`: Starting basepair, can be an Int or String
+ `end_bp`: Ending basepair, can be an Int or String

# Optional inputs
+ `min_maf`: Minimum minor allele frequency. Only variants with alternate allele
    frequency between `[min_maf, 1-min_maf]` is kept. Default `min_maf=0.01`
+ `snps_to_keep`: Vector of SNP positions to import. If both `snps_to_keep` and
    `min_maf` are specified, only SNPs whose position is listed in `snps_to_keep`
    whose minor allele frequency exceeds `min_maf` will be kept. 
+ `min_eigval`: Smallest eigenvalue allowed in Σ. All eigenvalues smaller than 
    `min_eigval` will be set to `min_eigval` (default `1e-5`). This option is only
    used if `enforce_psd=true`.
+ `enforce_psd`: LD data stored in Pan-UKB or gnomAD LD panels only includes the
    upper triangular portion. If `enforce_psd` is true, we will copy the upper 
    triangular portion to the lower triangular portion, force eigenvalues to be
    above `min_eigval`, and then scale the covariance matrix into a correlation
    matrix (default `true`)
+ `rk`: An `Int`. If specified and `enforce_psd=true`, we will keep only the top
    `rk` eigenvalues of imported `Σ` before truncating the rest to `min_eigval`.
    This experiments with the "truncated-SVD" approach of 
    `https://www.sciencedirect.com/science/article/pii/S0002929717303919#bib19`
    (default `Inf`)

# Note
Make sure `start_bp` and `end_bp` is from the same human genome build as the 
LD matrices. One can verify build information via `bm.build`. Both Pan-UKBB and 
gnomAD uses hg19 (i.e. GRCh37).
"""
function get_block(bm::HailBlockMatrix, chr::Union{String, Int}, 
    start_bp::Union{String, Int}, end_bp::Union{String, Int};
    min_maf::Real = 0.01, snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing,
    min_eigval = 1e-5, enforce_psd::Bool=true, rk=Inf,
    )
    # convert start_bp/end_bp to Int and chr to String
    start_bp = typeof(start_bp) == Int ? start_bp : parse(Int, start_bp)
    end_bp = typeof(end_bp) == Int ? end_bp : parse(Int, end_bp)
    chr = typeof(chr) == String ? chr : string(chr)
    # search for starting/ending position
    chr_range = findall(x -> x == chr, bm.info[!, "chr"])
    chr_pos = @view(bm.info[chr_range, "pos"])
    end_bp < chr_pos[1] && error("end_bp occurs before first SNP in chr$chr of bm")
    start_bp > chr_pos[end] && error("start_bp occurs after last SNP in chr$chr of bm")
    # import LD block and snp information
    chr_offset = chr_range[1] - 1
    idx_start = searchsortedfirst(chr_pos, start_bp) + chr_offset
    idx_end = searchsortedlast(chr_pos, end_bp) + chr_offset
    sigma = bm[idx_start:idx_end, idx_start:idx_end]
    df = bm.info[idx_start:idx_end, :]
    # keep SNPs with MAF above threshold and only SNPs that are requested to keep
    idx = findall(x -> 1-min_maf ≥ x ≥ min_maf, bm.info[idx_start:idx_end, "AF"])
    if !isnothing(snps_to_keep) 
        intersect!(idx, indexin(snps_to_keep, df[!, "pos"]))
    end
    Sigma = sigma[idx, idx]
    df = df[idx, :]
    # ensure Sigma is PSD
    if enforce_psd
        LinearAlgebra.copytri!(Sigma, 'U') # copy upper triangular part to lower triangular
        evals, evecs = eigen(Sigma)
        if length(evals) > rk
            evals[1:length(evals)-rk] .= min_eigval
        end
        evals[findall(x -> x < min_eigval, evals)] .= min_eigval
        Sigma = evecs * Diagonal(evals) * evecs'
        cov2cor!(Sigma, sqrt.(diag(Sigma))) # scale to correlation matrix
    end
    return Sigma, df
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
    df = read_variant_index_tables(ht_files)
    build = _extract_genome_build(ht_files)
    bm = HailBlockMatrix(bm_files, df, build,
        hail_linalg.BlockMatrix.read(bm_files))

    if size(bm, 1) != size(bm.info, 1)
        @warn "Dimension of variant index files does not agree " * 
            "with block matrix files, retrying the pre-processing step..."
        df = read_variant_index_tables(ht_files, reprocess=true)
        bm = HailBlockMatrix(bm_files, df, build,
            hail_linalg.BlockMatrix.read(bm_files))
        size(bm, 1) != size(bm.info, 1) && error(
            "Dimension of variant index files still does not agree " * 
            "with block matrix files, please try re-downloading the files."
            )
    end

    return bm
end

"""
    read_variant_index_tables(ht_file::String; [alt_allele_header], [reprocess])

Read variant index hail tables into a `DataFrame`. The first time this function
gets called, we will read the original `.ht` files into memory and write the
result to a tab-separated value `.tsv` file into the same directory as `ht_file`, 
and we will also create a comma separated file ending in `.csv` which pre-process
the `.tsv` file for easier reading. If `reprocess=true`, we will re-run the
preprocessing step.

# Note
+ `rsIDs`: Some panels (e.g. Pan-UKB) also contain rsID, but others (e.g. gnomAD) does not,
    so rsID may not be available. 
+ `AF`: Alternate Allele frequency is sometimes available in the `.tsv` file in
    its own column (e.g. Pan-UKB), but other times it's hidden inside a column 
    of the `.tsv` file (e.g. in gnomAD it is under the header `pop_freq` which
    includes `AF` but also other information). In the later case, we will assume
    the column in `.tsv` has header name `alt_allele_header` and is JSON
    formatted, and further that alt-allele freq is available as `AF` key. See
    [`_extract_alternate_allele_freq`](@ref)
"""
function read_variant_index_tables(ht_file::String; 
    alt_allele_header = "pop_freq", reprocess::Bool=false)
    isdir(ht_file) || error("$ht_file is not a directory")
    tsv_file = joinpath(ht_file, "variant.ht.tsv")
    csv_file = joinpath(ht_file, "variant.ht.csv")
    if reprocess || !isfile(csv_file)
        # first export basic hail table info
        println("Exporting hail table data to file $csv_file"); flush(stdout)
        hail_table = hail.read_table(ht_file)
        hail_table.export(tsv_file)
        # read chr/pos/ref/alt
        df = CSV.read(tsv_file, DataFrame)
        locus = split.(df[!, "locus"], ':')
        chr = [locus[i][1] for i in eachindex(locus)] |> Vector{String}
        pos = [parse(Int, locus[i][2]) for i in eachindex(locus)]
        ref, alt = _extract_ref_alt_alleles(df[!, "alleles"])
        # check that within each chr, basepair positions are sorted
        for c in unique(chr)
            idx = findall(x -> x == c, chr)
            issorted(@view(pos[idx])) || 
                error("Basepair positions in chr $chr not sorted!")
        end
        # alt allele freq
        if !("AF" in names(df))
            df[!, :AF] = 
                _extract_alternate_allele_freq(df, header_name=alt_allele_header)
        end
        # save csv file
        df[!, :chr] = chr
        df[!, :pos] = pos
        df[!, :ref] = ref
        df[!, :alt] = alt
        select!(df, Not([:idx, :locus, :alleles]))
        if alt_allele_header in names(df)
            select!(df, Not(Symbol(alt_allele_header)))
        end
        CSV.write(csv_file, df)
    end
    return CSV.read(csv_file, DataFrame)
end

"""
    _extract_alternate_allele_freq(df::DataFrame; header_name = "pop_freq")

Internal helper function that extracts alternate allele frequency from the 
`ht_file` assuming this information is present in the variant index files 
as JSON format with header `header_name`. Only `gnomAD` panel will actually
reach here since Pan-UKB has `AF` as a separate field in its variant index files

For example, in gnomAD, we have
```julia
julia> tsv_file = "gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht/variant.ht.tsv"
julia> df = CSV.read(tsv_file, DataFrame)
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
function _extract_alternate_allele_freq(df::DataFrame; header_name = "pop_freq")
    af = zeros(Float64, size(df, 1))
    for (i, s) in enumerate(df[!, header_name])
        af[i] = JSON.parse(s)["AF"]
    end
    return af
end

"""
    _extract_ref_alt_alleles(alleles::Vector{String})

Internal helper function to extract the ref/alt alleles from hail tables. Each 
element of `alleles` should be a `string` that has 2 alleles wrapped in quotations,
as in the case for Pan-UKB and gnomAD. The first is assumed to be reference and 
the second is alt. As an example, `["CCTAA","C"]` with give `ref = "CCTAA"` and 
`alt = "C"`. 
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

"""
    _extract_genome_build(ht_file::String)

Internal helper function to extracts the human genome build information from 
the variant index folder. We assume this information is stored in a file called 
`metadata.json` or `metadata.json.gz` in the following format `Locus(XXX)` where
`XXX` is the reference build. If this format is not detected, the genome build
information will be "Unknown"
"""
function _extract_genome_build(ht_file::String)
    locus_dat = ""
    if isfile(joinpath(ht_file, "metadata.json"))
        locus_dat = _read_metadata(joinpath(ht_file, "metadata.json"))["table_type"]
    elseif isfile(joinpath(ht_file, "metadata.json.gz"))
        locus_dat = _read_metadata(joinpath(ht_file, "metadata.json.gz"))["table_type"]      
    end
    # search for pattern Locus(XXX)
    regex = r"Locus\((.*?)\)"
    return occursin(regex, locus_dat) ? match(regex, locus_dat).captures[1] : "Unknown"
end

"""
    _read_metadata(metafile::String)

Internal helpder function that reads the metadata file and returns a dictionary.
`metafile` should end with `.json` or `.json.gz`
"""
function _read_metadata(metafile::String)
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
