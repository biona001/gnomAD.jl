# Downloading and Parsing Hail LD matrices

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biona001.github.io/EasyLD.jl/dev/) | [![build Actions Status](https://github.com/biona001/EasyLD.jl/workflows/CI/badge.svg)](https://github.com/biona001/EasyLD.jl/actions) | [![codecov](https://codecov.io/gh/biona001/EasyLD.jl/branch/master/graph/badge.svg?token=YyPqiFpIM1)](https://codecov.io/gh/biona001/EasyLD.jl) |

This is a Julia package that helps users download and read blocks of LD (linkage disequilibrium) matrices into memory. This package is inspired by an [existing pipeline](https://github.com/aaronsossin/gnomAD_LD_Easy_Querying) curated by Aaron Sossin. 

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

## gnomAD: downloading LD matrices and Variant index files

Downloading the LD matrices for a specific population:
```julia
julia> using EasyLD
julia> population = "nfe"
julia> outdir = "/Users/biona001/.julia/dev/EasyLD/data"
julia> download_gnomad_LD_matrices(population, outdir, start_from=1, num_files=3)

Progress: 100%|█████████████████████████████████████████| Time: 0:00:08
```
1. The result would be saved into `outdir` directory with name `gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm` which itself is a directory. 
2. For multithreaded downloads, one can use `start_from` and `num_files` keywords to control how many matrices to download (not specifying will download all matrices). A progress meter is automatically displayed. 

To check how many file there are:
```julia
julia> get_gnomad_filenames(population, join=false)

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

Downloading variant index hail tables: (this tells us chr/pos/ref/alt/SNP-name/alt-allele-frequency)

```julia
julia> download_gnomad_variant_index_tables(population, outdir)
```

The result will be stored as `gnomad.genomes.r2.1.1.$population.common.adj.ld.variant_indices.ht` in `outdir`

## Pan-UKBB: downloading LD matrices and Variant index files

Downloading the LD matrices for a specific population:
```julia
julia> using EasyLD
julia> population = "EUR"
julia> outdir = "/Users/biona001/.julia/dev/EasyLD/data"
julia> download_ukb_LD_matrices(population, outdir, start_from=1, num_files=10)

Progress: 100%|█████████████████████████████████████████| Time: 0:03:48
```
1. The result would be saved into `outdir` directory with name `UKBB.EUR.ldadj.bm` which itself is a directory. 
2. For multithreaded downloads, one can use `start_from` and `num_files` keywords to control how many matrices to download (not specifying will download all matrices). A progress meter is automatically displayed. 


To check how many file there are:
```julia
julia> get_ukb_filenames("EUR", join=false)

124281-element Vector{String}:
 "part-000000-44-0-0-22f828c8-17c3-7c3c-1fa1-1fc113144aca"
 "part-000001-44-1-0-d35bc353-c9b7-28e7-315b-f1bcd8d9e50b"
 "part-000002-44-2-2-7dc4e8fd-75f7-26d0-4b9c-65628f25cf34"
 "part-000003-44-3-0-74be7ee3-ba85-f9ce-4add-9638cf387f8a"
 "part-000004-44-4-0-b11cbd6a-50e2-4e5e-ae70-932ff3aebd5f"
 "part-000005-44-5-0-104f3537-b7e4-0959-555b-7e65ca2b00ae"
 "part-000006-44-6-0-f99bebaa-7361-9c3e-65a0-2700858c2b28"
 "part-000007-44-7-0-90293a19-7ae6-f6a6-1141-48f62ee6f8b1"
 "part-000008-44-8-0-847d7052-f784-0c5f-c4f8-d16df8745a23"
 "part-000009-44-9-0-2b61d351-aeda-4c03-52aa-cb380599219d"
 "part-000010-44-10-0-90b9d808-fa98-2727-ca3e-d63676464185"
 "part-000011-44-11-0-7a15dbc9-b5c7-3500-8258-886bbf569e35"
 "part-000012-44-12-0-73b49602-7c17-3f3b-860c-6aaa8b25f36b"
 ⋮
 "part-124226-44-124226-0-9948291d-f1c5-775a-4ed3-7d75c9e7d007"
 "part-124227-44-124227-0-43ffd695-b0c7-577e-0d5a-8ea5f3ed84d0"
 "part-124228-44-124228-0-d22a4e30-04f2-82aa-dfd2-4d3ae95e57b4"
 "part-124229-44-124229-0-12730aef-d105-61a3-5342-f20ed3bf3645"
 "part-124230-44-124230-0-984e649b-fc68-5397-3b47-52b76e87c2c6"
 "part-124231-44-124231-0-9fd5803f-f5c7-1c4c-a5fa-77b76dc1f52d"
 "part-124232-44-124232-0-b891bf57-4536-661b-d679-10aac2532957"
 "part-124233-44-124233-0-1c6f8f3f-9337-fc25-53ef-63e76c5a37cf"
 "part-124234-44-124234-0-ad9c21cc-9c7c-0ccb-9335-921cc35ae152"
 "part-124235-44-124235-0-0f31fad5-1ba4-79b3-dff0-cb17d6b694ed"
 "part-124236-44-124236-0-98383aba-a5c1-8093-c27f-001aae4730d3"
 "part-124237-44-124237-1-36053042-3e3b-0495-e362-d7e756ac8f5a"
```

Downloading variant index hail tables: (this tells us chr/pos/ref/alt/SNP-name/alt-allele-frequency)

```julia
julia> download_ukb_variant_index_tables(population, outdir)
```

The result will be stored as `UKBB.$population.ldadj.variant.ht` in `outdir`

## Reading BlockMatrix into memory

Read block of data into memory: 

```julia
julia> using EasyLD
julia> bm_file = "/Users/biona001/.julia/dev/EasyLD/data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm"
julia> ht_file = "/Users/biona001/.julia/dev/EasyLD/data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht"
julia> bm = hail_block_matrix(bm_file, ht_file); # the ';' avoids displaying a few entries of bm, which takes ~0.1 seconds per entry

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

# read a specific region using chr and basepair positions
julia> chr = "1"
julia> start_pos = 10146
julia> end_pos = 10181
julia> sigma = get_block(bm, chr, start_pos, end_pos)
5×5 Matrix{Float64}:
 1.0  0.0326361  -0.0181891   0.0166341    0.00802797
 0.0  1.0        -0.00041394  0.00540477   0.00232164
 0.0  0.0         1.0         0.0177243    0.0
 0.0  0.0         0.0         1.0         -0.0073965
 0.0  0.0         0.0         0.0          1.0
```

## Reading variant index files into a DataFrame

If one wishes, the `read_variant_index_tables` function reads the variant index files into a `DataFrame`. Because variant index files are small (at most a few GB), we will export it as a human-readable `.tsv` file to the `.ht` directory for faster reading. This is only done once. 

```julia
julia> ht_file = "gnomad.genomes.r2.1.1.nfe_test.common.adj.ld.variant_indices.ht"
julia> df = read_variant_index_tables(ht_file)

14207204×4 DataFrame
      Row │ locus        alleles                            pop_freq                           idx      
          │ String15     String                             String                             Int64    
──────────┼─────────────────────────────────────────────────────────────────────────────────────────────
        1 │ 1:10146      ["AC","A"]                         {"AC":861,"AF":0.409220532319391…         0
        2 │ 1:10151      ["TA","T"]                         {"AC":3,"AF":0.00810810810810810…         1
        3 │ 1:10177      ["A","C"]                          {"AC":13,"AF":0.1444444444444444…         2
        4 │ 1:10178      ["CCTAA","C"]                      {"AC":2,"AF":0.02325581395348837…         3
        5 │ 1:10181      ["A","T"]                          {"AC":2,"AF":0.00617283950617283…         4
        6 │ 1:10250      ["A","C"]                          {"AC":3,"AF":0.03260869565217391…         5
        7 │ 1:10257      ["A","C"]                          {"AC":10,"AF":0.0595238095238095…         6
        8 │ 1:10327      ["T","C"]                          {"AC":23,"AF":0.1982758620689655…         7
        9 │ 1:10329      ["AC","A"]                         {"AC":8,"AF":0.05,"AN":160,"homo…         8
       10 │ 1:10333      ["CT","C"]                         {"AC":5,"AF":0.00720461095100864…         9
       11 │ 1:10347      ["AACCCT","A"]                     {"AC":11,"AF":0.0066105769230769…        10
       12 │ 1:10397      ["CCCCTAA","C"]                    {"AC":34,"AF":0.0072930072930072…        11
       13 │ 1:10403      ["ACCCTAACCCTAACCCTAACCCTAACCCTA…  {"AC":1751,"AF":0.24689791314156…        12
    ⋮     │      ⋮                       ⋮                                  ⋮                     ⋮
 14207193 │ X:155259328  ["G","C"]                          {"AC":83,"AF":0.0166132906325060…  14207192
 14207194 │ X:155259339  ["G","A"]                          {"AC":144,"AF":0.018191005558362…  14207193
 14207195 │ X:155259571  ["T","G"]                          {"AC":490,"AF":0.338397790055248…  14207194
 14207196 │ X:155259572  ["G","A"]                          {"AC":24,"AF":0.0074906367041198…  14207195
 14207197 │ X:155259589  ["T","G"]                          {"AC":7,"AF":0.00527108433734939…  14207196
 14207198 │ X:155259601  ["T","A"]                          {"AC":15,"AF":0.0103163686382393…  14207197
 14207199 │ X:155259762  ["G","A"]                          {"AC":44,"AF":0.1571428571428571…  14207198
 14207200 │ X:155259767  ["T","TA"]                         {"AC":35,"AF":0.1,"AN":350,"homo…  14207199
 14207201 │ X:155259827  ["T","G"]                          {"AC":208,"AF":0.338762214983713…  14207200
 14207202 │ X:155259894  ["AG","A"]                         {"AC":23,"AF":0.0592783505154639…  14207201
 14207203 │ X:155259894  ["AGGGGTTAG","A"]                  {"AC":85,"AF":0.2190721649484536…  14207202
 14207204 │ X:155259920  ["AG","A"]                         {"AC":2,"AF":0.02040816326530612…  14207203
```

## Requirements
+ A working local installation of Python packages `hail` and `numpy`.
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
