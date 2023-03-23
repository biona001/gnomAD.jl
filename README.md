# Downloading and Parsing Hail LD matrices

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biona001.github.io/EasyLD.jl/dev/) | [![build Actions Status](https://github.com/biona001/EasyLD.jl/workflows/CI/badge.svg)](https://github.com/biona001/EasyLD.jl/actions) | [![codecov](https://codecov.io/gh/biona001/EasyLD.jl/branch/master/graph/badge.svg?token=YyPqiFpIM1)](https://codecov.io/gh/biona001/EasyLD.jl) |

This is a Julia package that helps users download and read blocks of LD (linkage disequilibrium) matrices stored as [HailBlockMatrix format](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#blockmatrix) into memory. This package is inspired by an [existing pipeline](https://github.com/aaronsossin/gnomAD_LD_Easy_Querying) curated by Aaron Sossin. 

We tested this package for processing the following 

+ [gnomAD LD matrices](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium)
+ [Pan-UKBB LD matrices](https://pan-dev.ukbb.broadinstitute.org/docs/hail-format/index.html)

## Installation

This resource uses [AWSS3.jl](https://github.com/JuliaCloud/AWSS3.jl) to download files from Amazon servers. Thus, one needs to register an account with AWS first. 

To install this package, download [Julia](https://julialang.org/downloads/). Within Julia, execute the following
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/EasyLD.jl")
```

Requirements:
+ A working local installation of Python package `hail`
+ Java version 8 or 11, according to [hail documentation](https://hail.is/docs/0.2/getting_started_developing.html). 

We currently rely on the Hail python package to read the matrix and use the PyCall 
package to pass the result into Julia, but in the future, we will consider developing 
a native Julia parser to avoid the reliance on Python and Hail altogether. 

## Trouble shooting

+ `PyCall.jl` cannot find `hail` and `numpy` even though I installed them? Julia's `PyCall.jl` by default installs its own local copy of python that is different than the system default. Thus, after installing `PyCall.jl`, within Julia I had to do `ENV["PYTHON"] = "PATH OF PYTHON EXECUTABLE"` and rebuild `PyCall.jl` via `using Pkg; Pkg.build("PyCall")`. This will point the python version in `PyCall.jl` to be the default python on your system.
+ `Intel MKL FATAL ERROR: Cannot load libmkl_intel_thread.dylib`
Try following suggestions in [this post](https://github.com/JuliaPy/PyPlot.jl/issues/315). In particular, the following worked for me:
```julia
# within Julia
using Conda
Conda.rm("mkl")
Conda.add("nomkl")
Conda.add("hail")
```
+ `OutOfMemoryError: Java heap space`: Try increasing the java heap size by adding `export _JAVA_OPTIONS="-Xmx24g"` to your `.bash_profile` file. Here the `-Xmx24g` implies maximum heap size of 24 GB, which you are free to adjust. 
