# EasyLD.jl

This is a Julia package that helps users download and read blocks of LD (linkage disequilibrium) matrices stored as [Hail Block Matrix format](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#blockmatrix) into memory. This package is inspired by an [existing pipeline](https://github.com/aaronsossin/gnomAD_LD_Easy_Querying) curated by Aaron Sossin.

We tested this package for processing the following

1. [gnomAD LD matrices](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium)
2. [Pan-UKBB LD matrices](https://pan-dev.ukbb.broadinstitute.org/docs/hail-format/index.html)

## Installation

This resource uses [AWSS3.jl](https://github.com/JuliaCloud/AWSS3.jl) to download files from Amazon servers. Thus, one needs to register an account with AWS first.

To install this package, download [Julia](https://julialang.org/downloads/). Within Julia, execute the following

```julia
using Pkg
pkg"add https://github.com/biona001/EasyLD.jl"
```

## Manual Outline

```@contents
Pages = [
    "man/api.md"
    "man/pan_UKB.md"
    "man/gnomAD.md"
]
Depth = 2
```
