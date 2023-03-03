# Downloading and Parsing gnomAD LD matrices

This is a Julia package that helps users download the [gnomAD LD matrices](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium) and read chunks of it into memory. 

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, execute the following
```julia
using Pkg
Pkg.add("https://github.com/biona001/gnomAD.jl")
```

## Quick start

Downloading the LD matrices for a specific population is achieved with `download_LD_matrices`
function. `start_from` and `num_files` keywords allow one to control how many matrices to
download (not specifying will download all matrices). A progress meter is automatically displayed. 
```julia
julia> using gnomAD
julia> population = "nfe"
julia> outdir = "/Users/biona001/.julia/dev/gnomAD/data"
julia> download_LD_matrices(population, outdir, start_from=1, num_files=3)

Progress: 100%|█████████████████████████████████████████| Time: 0:00:08
```
The result would be saved into `outdir` directory with name `gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm` which itself is a directory. 

To check how many file there are:
```julia
julia> get_all_filenames(population, join=false)

10023-element Vector{String}:
 "part-00000-14-0-0-7e282c91-9060-db76-1d4a-a02877b62910"
 "part-00001-14-1-0-2bb43364-b9d9-7048-c4e0-74c951357d81"
 "part-00002-14-2-0-998bf016-ab2f-aaf9-ce80-7888e4ee7524"
 "part-00003-14-3-0-5ad64a0b-ecb1-cf91-69ab-e3ee2fcdef51"
 "part-00004-14-4-0-7a03381b-dc38-5f21-f8dd-ad869dd1f340"
 "part-00005-14-5-0-d78f2159-d71c-adc2-47eb-5af6a6c34015"
 "part-00006-14-6-0-dc8b5638-9246-15b8-579e-d14ff6a69646"
 "part-00007-14-7-0-f41e5872-2ed6-36d2-e625-bb6d6261c6bc"
 "part-00008-14-8-0-27a9b5b4-b329-d7a8-9d86-f3184776cb09"
 "part-00009-14-9-0-68655a89-424f-0cb7-b35c-583b52c859e3"
 "part-00010-14-10-0-a60938ee-75ec-91c9-91f3-fbcffa73a20e"
 "part-00011-14-11-0-3be840ef-ae93-3775-47fe-fc33eedb9911"
 "part-00012-14-12-0-8cbce5b3-023a-b252-52c2-cc9c717d0ec5"
 ⋮
 "part-10011-14-10011-0-1ed25fc8-64a3-7084-60be-717f4c447fac"
 "part-10012-14-10012-0-342cfc1b-d185-4c2d-ad6e-d9595c5b9071"
 "part-10013-14-10013-0-004a90a8-ec7c-0664-d592-47686c3df576"
 "part-10014-14-10014-0-fdad432a-59ea-ef55-aa6c-dfebf128272c"
 "part-10015-14-10015-0-5216889b-97ba-8e71-8d77-27c2e41a7a38"
 "part-10016-14-10016-0-955e1d70-ef0f-5dc7-cf3f-55c2d7626c53"
 "part-10017-14-10017-0-14ea97f3-00af-93ba-ffa7-f6a0ff9c6541"
 "part-10018-14-10018-0-2cc71e56-18a8-de9e-85c8-1513ec63cceb"
 "part-10019-14-10019-0-16de6180-8034-d811-4875-a91b2da0dad8"
 "part-10020-14-10020-0-1a65d50c-ffe5-2f10-c510-60e0cfa5b186"
 "part-10021-14-10021-0-fc83fd8d-e05e-6955-fc58-bcf4296985ba"
 "part-10022-14-10022-0-216231b4-e554-eac9-d165-ffa8a90fef05"

```

Read block of gnomAD data into memory: 

```julia
julia> using gnomAD
julia> data = "/Users/biona001/.julia/dev/gnomAD/data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm"
julia> bm = hail_block_matrix(data); # need a ';' to avoid displaying a few entries of bm, which takes ~0.1 seconds per entry

julia> size(bm) # matrix dimension
(14207204, 14207204)

julia> bm[1, 1] # read a single entry (note: julia uses 1-based indexing)
0.9999999999999998

julia> bm[1:10000, 1:10000] # read first 10k by 10k block (takes roughly 7 seconds)
10000×10000 Matrix{Float64}:
 1.0  0.0326361  -0.0181891   …   0.0         0.0         0.0
 0.0  1.0        -0.00041394      0.0         0.0         0.0
 0.0  0.0         1.0             0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 0.0  0.0         0.0         …   0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 0.0  0.0         0.0             0.0         0.0         0.0
 ⋮                            ⋱                          
 0.0  0.0         0.0            -0.146398    0.193441   -0.0275316
 0.0  0.0         0.0            -0.0631539   0.0803796   0.0113878
 0.0  0.0         0.0            -0.0397994   0.0521778  -0.00356272
 0.0  0.0         0.0            -0.0309998   0.0109042  -0.00833901
 0.0  0.0         0.0         …  -0.112318   -0.285862   -0.0291041
 0.0  0.0         0.0            -0.152088    0.362752   -0.0176148
 0.0  0.0         0.0             1.0         0.22565     0.0399238
 0.0  0.0         0.0             0.0         1.0        -0.0670657
 0.0  0.0         0.0             0.0         0.0         1.0


# arbitrary slicing works but is very slow
 julia> bm[1:3, 1:2:100] # ~22 seconds
 3×50 Matrix{Float64}:
 1.0  -0.0181891   0.00802797  -0.0421512  …  0.0143534     0.0172358
 0.0  -0.00041394  0.00232164  -0.0112066     5.39059e-5   -0.0188801
 0.0   1.0         0.0         -0.0282011     0.000733417  -0.0277042
```

## Requirements
+ The [PyCall package](https://github.com/JuliaPy/PyCall.jl) to call Python
+ A working local installation of Python packages `hail` and `numpy` already installed.
    After installing `PyCall.jl`, I had to change 
    `ENV["PYTHON"] = "PATH OF PYTHON EXECUTABLE"` and rebuild `PyCall.jl` via
    `using Pkg; Pkg.build("PyCall")`

We currently rely on the Hail python package to read the matrix and use the PyCall 
package to pass the result into Julia, but in the future, we will consider developing 
a native Julia parser to avoid the reliance on Python and Hail altogether. 

## Trouble shooting

+ `Intel MKL FATAL ERROR: Cannot load libmkl_intel_thread.dylib`
Try following suggestions in [this post](https://github.com/JuliaPy/PyPlot.jl/issues/315). In particular, the following worked for me:
```julia
# within Julia
using Conda
Conda.rm("mkl")
Conda.add("nomkl")
Conda.add("hail")
```
+ `OutOfMemoryError: Java heap space`: Try increasing the java heap size by `export _JAVA_OPTIONS="-Xmx24g"` to your `.bash_profile` file. Here the `-Xmx24g` implies maximum heap size of 24 GB, which you are free to adjust. 
