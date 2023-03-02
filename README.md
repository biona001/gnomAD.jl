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
julia> population = "afr"
julia> outdir = "/Users/biona001/.julia/dev/gnomAD/data"
julia> download_LD_matrices(population, outdir, start_from=1, num_files=3)

Progress: 100%|█████████████████████████████████████████| Time: 0:00:08
```

To check how many file there are, we provide a convenient function
```julia
julia> get_all_filenames(population, join=false)

22432-element Vector{String}:
 "part-00000-14-0-0-4eff4018-c6cf-a264-f2e5-e17fb42bfd6e"
 "part-00001-14-1-0-c4c9820f-a915-487a-4056-c855f7442286"
 "part-00002-14-2-0-8af33c53-352d-4d45-cc3b-59a7956fad1a"
 "part-00003-14-3-0-ae814c01-2af9-238c-f100-485d9cb6379c"
 "part-00004-14-4-0-17d481bd-5570-5237-2390-4b1afc51a404"
 "part-00005-14-5-0-3f7ad064-938c-32a9-1ca5-3690f3f95f79"
 "part-00006-14-6-0-189ff387-742c-54c1-1b52-686e3e686b34"
 "part-00007-14-7-0-3d71f087-403a-ef35-4b87-ecbde776a0c2"
 "part-00008-14-8-0-757a13d6-caa1-d0a9-6999-807af9e749f5"
 "part-00009-14-9-0-d599eec5-cb98-fb75-baf7-8869481f131f"
 "part-00010-14-10-0-0779bd48-d897-affc-35de-5f77861cce90"
 "part-00011-14-11-0-0fe8305a-ff43-e5fe-2190-c037ff5a1377"
 "part-00012-14-12-0-b0cafccc-9ad8-f8c9-5534-874c0d8af261"
 ⋮
 "part-22404-14-22404-0-0db2ff54-cb75-e7ab-7ee3-308b22b623c3"
 "part-22405-14-22405-0-6cfed31e-d521-cbff-d756-04a6844b4ff2"
 "part-22406-14-22406-0-7d248ee1-b887-eb00-c7a2-37e66573096c"
 "part-22407-14-22407-0-130eff98-f6e3-848e-9654-32307de33f53"
 "part-22408-14-22408-0-06354349-b83f-2265-39eb-9187e0b7d50b"
 "part-22409-14-22409-0-caff5e31-d811-c65d-4d5f-a249ba444a03"
 "part-22410-14-22410-0-1c5003b0-fe39-0c78-005c-8b9978cf5427"
 "part-22411-14-22411-0-99d9b653-fd96-0a28-0c18-506ee744d098"
 "part-22412-14-22412-0-9fb04ec2-034d-9a41-9b46-3f6b80b32a2d"
 "part-22413-14-22413-0-fa53677f-dfa7-413c-4190-1273c8057717"
 "part-22414-14-22414-0-66a0f19c-d8ab-876a-b70e-4e6741411f69"
 "part-22415-14-22415-0-39a245e9-c9f3-8275-ba4e-ffd24fe3c891"
```

Finally, here is how to read a block of gnomAD data into memory. We currently rely
on the Hail python package to read the matrix and use the PyCall package to pass
it into Julia, but in the future, we will consider developing a native Julia 
parser to avoid the reliance on Python and Hail altogether. 

```julia
julia> using gnomAD
julia> data = "gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm"
julia> bm = hail_block_matrix(data);

julia> bm[1, 1] # read first entry (note: julia uses 1-based indexing)
0.9999999999999998

julia> bm[1:1000, 1:1000] # takes roughly 0.5 seconds
1000×1000 Matrix{Float64}:
 1.0  0.0326361  -0.0181891   0.0166341   …   0.00767459    0.0104462
 0.0  1.0        -0.00041394  0.00540477     -0.10473      -0.00184638
 0.0  0.0         1.0         0.0177243      -0.000167982  -0.00929269
 0.0  0.0         0.0         1.0            -2.1684e-19   -0.00329412
 0.0  0.0         0.0         0.0             1.2618e-5    -0.000859109
 0.0  0.0         0.0         0.0         …   2.98156e-19   0.0507852
 0.0  0.0         0.0         0.0             2.88077e-5   -0.0117684
 0.0  0.0         0.0         0.0            -0.000116338  -0.0104104
 0.0  0.0         0.0         0.0             2.67669e-5   -0.00728978
 0.0  0.0         0.0         0.0             2.331e-5     -0.00253932
 0.0  0.0         0.0         0.0         …   0.00150567   -0.00353229
 0.0  0.0         0.0         0.0             0.00509611    0.0109818
 0.0  0.0         0.0         0.0            -0.0126056     0.0126193
 ⋮                                        ⋱                
 0.0  0.0         0.0         0.0            -0.0182409    -0.0106091
 0.0  0.0         0.0         0.0            -0.00116818    0.00235988
 0.0  0.0         0.0         0.0         …  -0.00204165   -0.0130619
 0.0  0.0         0.0         0.0             0.00840741   -0.00771116
 0.0  0.0         0.0         0.0             0.983214      0.00734632
 0.0  0.0         0.0         0.0             0.0156904    -0.00610452
 0.0  0.0         0.0         0.0             0.00396482   -0.0110234
 0.0  0.0         0.0         0.0         …  -0.00290343   -0.0109804
 0.0  0.0         0.0         0.0             0.00723443   -0.013669
 0.0  0.0         0.0         0.0             0.99998       0.00738016
 0.0  0.0         0.0         0.0             1.0           0.00737766
 0.0  0.0         0.0         0.0             0.0           1.0

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
