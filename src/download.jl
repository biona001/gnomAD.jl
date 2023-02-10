"""
    download_LD_matrices(population::String, outdir::String, num_files = Inf)

Downloads the LD matrices for specified population and saves result to
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
function download_LD_matrices(
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
    for file in ["metadata.json", "_SUCCESS"]
        s3_get_file(bucket, joinpath(bmpath, file), joinpath(datadir, file))
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
    get_all_filenames(population::String)

Returns a list of file names that would be downloaded for the specified population.
This function will not download anything. If `join=true`, the full path will be returned
"""
function get_all_filenames(population::String; join::Bool=false)
    bucket = "gnomad-public-us-east-1/release"
    path = "2.1.1/ld/gnomad.genomes.r2.1.1.$population.common.adj.ld.bm/parts"
    dir = joinpath("s3://", bucket, path)
    return readdir(S3Path(dir * "/"), join=join)
end
