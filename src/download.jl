"""
    download_gnomad_LD_matrices(population::String, outdir::String, [start_from], [num_files])

Downloads the LD matrices from gnomAD panel for specified population and saves result to
`outdir`. See `https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium`.

# Inputs
+ `population`: Options include
    - `afr` for African/African American population
    - `amr` for Latino/Admixed American population
    - `asj` for Ashkenazi Jewish population
    - `eas` for East Asian population
    - `fin` for European (Finnish) population
    - `nfe` for European (non-Finnish) population
    - `est` for Estonian population
    - `nwe` for North-western European population
    - `seu` for Southern European population
+ `outdir`: Directory for which files will be saved.
+ `start_from`: Int specifying which file to start downloading at
+ `num_files`: Int specifying how many files to download

# Multiple download sessions at once
Since there are lots of files to download, one can specify the number
of files to download by modifying the `num_files` files option. This function
will download the first `num_files` file counting from the `start_from` element 
of the output of `get_all_filenames()`. For example, if 
`myfiles = get_all_filenames("afr")`, then the function will download files 
from `myfiles[start_from]` through `myfiles[start_from+num_files-1]`
"""
function download_gnomad_LD_matrices(
    population::String,
    outdir::String;
    start_from::Int = 1,
    num_files = Inf
    )
    # initialize local directories
    datadir = joinpath(outdir, "gnomad.genomes.r2.1.1.$population.common.adj.ld.bm")
    partsdir = joinpath(datadir, "parts")
    isdir(partsdir) || mkpath(partsdir)
    # remote directories
    bucket = "gnomad-public-us-east-1/release"
    bmpath = "2.1.1/ld/gnomad.genomes.r2.1.1.$population.common.adj.ld.bm"
    # download metadata and success indicator file
    if start_from == 1
        for file in ["metadata.json", "_SUCCESS"]
            s3_get_file(bucket, joinpath(bmpath, file), joinpath(datadir, file))
        end
    end
    # .bm (LD matrix Hail Block Matrix) files are in bucket/bmpath/parts/*
    bmfiles = readdir(S3Path(joinpath("s3://", bucket, bmpath, "parts") * "/"))
    last_file = start_from+num_files-1 > length(bmfiles) ? 
        length(bmfiles) : start_from+num_files-1
    @showprogress for file in bmfiles[start_from:last_file]
        s3_get_file(bucket, 
            joinpath(bmpath, "parts", file), 
            joinpath(partsdir, file))
    end
    return nothing
end

"""
    download_ukb_LD_matrices(population::String, outdir::String, [start_from], [num_files])

Downloads the LD matrices from pan-UK-biobank for specified population and saves
result to `outdir`. See 
`https://pan-dev.ukbb.broadinstitute.org/docs/hail-format/index.html`.

# Inputs
+ `population`: Several options are available. These codes refer only to ancestry 
    groupings used in GWAS, not necessarily other demographic or self-reported data.
    - `EUR` for European ancestry
    - `CSA` for Central/South Asian ancestry
    - `AFR` for African ancestry
    - `EAS` for East Asian ancestry
    - `MID` for Middle Eastern ancestry
    - `AMR` for Admixed American ancestry
+ `outdir`: Directory for which files will be saved.
+ `start_from`: Int specifying which file to start downloading at
+ `num_files`: Int specifying how many files to download

# Multiple download sessions at once
Since there are lots of files to download, one can specify the number
of files to download by modifying the `num_files` files option. This function
will download the first `num_files` file counting from the `start_from` element 
of the output of `get_all_filenames()`. For example, if 
`myfiles = get_all_filenames("afr")`, then the function will download files 
from `myfiles[start_from]` through `myfiles[start_from+num_files-1]`
"""
function download_ukb_LD_matrices(
    population::String,
    outdir::String;
    start_from::Int = 1,
    num_files = Inf
    )
    # initialize local directories
    datadir = joinpath(outdir, "UKBB.$population.ldadj.bm")
    partsdir = joinpath(datadir, "parts")
    isdir(partsdir) || mkpath(partsdir)
    # remote directories
    bucket = "pan-ukb-us-east-1/ld_release"
    bmpath = "UKBB.$population.ldadj.bm"
    # download metadata and success indicator file
    if start_from == 1
        for file in ["metadata.json", "_SUCCESS"]
            s3_get_file(bucket, joinpath(bmpath, file), joinpath(datadir, file))
        end
    end
    # .bm (LD matrix Hail Block Matrix) files are in bucket/bmpath/parts/*
    bmfiles = readdir(S3Path(joinpath("s3://", bucket, bmpath, "parts") * "/"))
    last_file = start_from+num_files-1 > length(bmfiles) ? 
        length(bmfiles) : start_from+num_files-1
    @showprogress for file in bmfiles[start_from:last_file]
        s3_get_file(bucket, 
            joinpath(bmpath, "parts", file), 
            joinpath(partsdir, file))
    end
    return nothing
end

"""
    get_gnomad_filenames(population::String)

Returns a list of file names that would be downloaded for the specified population.
This function will not download anything. If `join=true`, the full path will be returned
"""
function get_gnomad_filenames(population::String; join::Bool=false)
    bucket = "gnomad-public-us-east-1/release"
    path = "2.1.1/ld/gnomad.genomes.r2.1.1.$population.common.adj.ld.bm/parts"
    dir = joinpath("s3://", bucket, path)
    return readdir(S3Path(dir * "/"), join=join)
end

"""
    get_ukb_filenames(population::String)

Returns a list of file names that would be downloaded for the specified population.
This function will not download anything. If `join=true`, the full path will be returned
"""
function get_ukb_filenames(population::String; join::Bool=false)
    bucket = "pan-ukb-us-east-1/ld_release"
    path = "UKBB.$population.ldadj.bm/parts"
    dir = joinpath("s3://", bucket, path)
    return readdir(S3Path(dir * "/"), join=join)
end

"""
    read_metadata(metafile::String)

Reads the metadata file and returns a dictionary
"""
function read_metadata(metafile::String)
    s = read(metafile) |> String
    return JSON.parse(s)
end

"""
    download_ukb_variant_index_tables(population::String, outdir::String)

Downloads variant index hail tables from pan-UK-biobank for specified population 
and saves result to `outdir`. `population` accepts the same population labels as 
`download_ukb_LD_matrices`. These files are typically at most a few GB.
"""
function download_ukb_variant_index_tables(population::String, outdir::String)
    isdir(outdir) || mkpath(outdir)
    bucket = "pan-ukb-us-east-1/ld_release"
    foldername = "UKBB.$population.ldadj.variant.ht"
    run(`aws s3 cp s3://$(joinpath(bucket, foldername)) $outdir --recursive`)
    println("success")
end

"""
    download_gnomad_variant_index_tables(population::String, outdir::String)

Downloads variant index hail tables from gnomAD panel for specified population 
and saves result to `outdir`. `population` accepts the same population labels as 
`download_gnomad_LD_matrices`. These files are typically at most a few GB.
"""
function download_gnomad_variant_index_tables(population::String, outdir::String)
    isdir(outdir) || mkpath(outdir)
    bucket = "gnomad-public-us-east-1/release"
    path = "2.1.1/ld/gnomad.genomes.r2.1.1.$population.common.adj.ld.variant_indices.ht"
    run(`aws s3 cp s3://$(joinpath(bucket, path)) $outdir --recursive`)
    println("success")
end
