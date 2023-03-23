
# gnomAD data

+ Link: https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium
+ `nfe` population totals ~555GB


```julia
using Revise
using EasyLD
using CSV
using DataFrames
using Statistics
using LinearAlgebra
```

## Downloading 
First check how many files there are:


```julia
get_gnomad_filenames("nfe", join=false)
```




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




```julia
population = "nfe"
outdir = "/Users/biona001/.julia/dev/EasyLD/data"
download_gnomad_LD_matrices(population, outdir, start_from=1, num_files=10)
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:00:29[39m


+ The result would be saved into outdir directory with name `gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm` which itself is a directory.

+ For multithreaded downloads, one can use `start_from` and `num_files` keywords to control how many matrices to download (not specifying will download all matrices). A progress meter is automatically displayed.

One also needs the variant index files which tells us chr/pos/ref/alt/SNP-name/alt-allele-frequency


```julia
population = "nfe"
outdir = "/Users/biona001/.julia/dev/EasyLD/data"
@time download_gnomad_variant_index_tables(population, outdir)
```

The result will be stored as `gnomad.genomes.r2.1.1.$population.common.adj.ld.variant_indices.ht` in `outdir`

## Reading LD panel with Matrix interface

The first time `hail_block_matrix` gets called, we do some pre-processing to the variant index files so subsequent calls will be faster. 


```julia
bm_file = "/Users/biona001/.julia/dev/EasyLD/data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm"
ht_file = "/Users/biona001/.julia/dev/EasyLD/data/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht"
@time bm = hail_block_matrix(bm_file, ht_file); # the ';' avoids displaying a few entries of bm, which takes ~0.1 seconds per entry
```

     33.650270 seconds (477.52 M allocations: 10.556 GiB, 12.67% gc time, 18.13% compilation time)


Check size of matrix


```julia
size(bm)
```




    (14207204, 14207204)



Read first 10000 by 10000 block into memory


```julia
Sigma = bm[1:10000, 1:10000]
```

    [Stage 0:=================================================>         (5 + 1) / 6]




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
     0.0  0.0         0.0         …   0.0         0.0         0.0
     0.0  0.0         0.0             0.0         0.0         0.0
     0.0  0.0         0.0             0.0         0.0         0.0
     ⋮                            ⋱                          
     0.0  0.0         0.0             0.368793    0.36933    -0.0120351
     0.0  0.0         0.0             0.789835    0.10943     0.0197234
     0.0  0.0         0.0         …  -0.167317    0.206943    0.0121849
     0.0  0.0         0.0            -0.146398    0.193441   -0.0275316
     0.0  0.0         0.0            -0.0631539   0.0803796   0.0113878
     0.0  0.0         0.0            -0.0397994   0.0521778  -0.00356272
     0.0  0.0         0.0            -0.0309998   0.0109042  -0.00833901
     0.0  0.0         0.0         …  -0.112318   -0.285862   -0.0291041
     0.0  0.0         0.0            -0.152088    0.362752   -0.0176148
     0.0  0.0         0.0             1.0         0.22565     0.0399238
     0.0  0.0         0.0             0.0         1.0        -0.0670657
     0.0  0.0         0.0             0.0         0.0         1.0



Check if the given block is PSD by computing its eigenvalues


```julia
eigvals(Symmetric(Sigma)) # Symmetric uses upper triangular portion of data
```




    10000-element Vector{Float64}:
      -2.5798836838860755
      -1.415006812481892
      -0.7674440501660297
      -0.6247787622218666
      -0.584928280772804
      -0.5435038109298787
      -0.5187701863140765
      -0.43166908059139875
      -0.41532154277852423
      -0.38813217585569226
      -0.3700017799653646
      -0.35667047980202393
      -0.3485121536359235
       ⋮
     109.56227091251577
     113.33023962897241
     127.38408283290424
     130.99148352578467
     132.82965433833886
     145.9473736360883
     155.27371456569082
     177.7213297373829
     194.1048422669041
     231.15057521884452
     293.44132970720796
     317.6410079580756



Arbitrary slicing works but is very slow


```julia
@time bm[1:3, 1:2:100] # read columns 1, 3, 5, ...
```

     30.278933 seconds (759 allocations: 13.219 KiB)





    3×50 Matrix{Float64}:
     1.0  -0.0181891   0.00802797  -0.0421512  …  0.0143534     0.0172358
     0.0  -0.00041394  0.00232164  -0.0112066     5.39059e-5   -0.0188801
     0.0   1.0         0.0         -0.0282011     0.000733417  -0.0277042



## Read in a block with `get_block`

One can also extract a block by specifying the chromosome and starting/ending basepair


```julia
chr = 1
start_pos = 57292
end_pos = 59193
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.0)
sigma
```




    8×8 Matrix{Float64}:
     1.0  0.0135236  -0.0249515  -0.0617588  …  -0.00321314  -0.00718484
     0.0  1.0         0.0167518   0.0209895      0.00877351   0.012266
     0.0  0.0         1.0        -0.0125702     -0.0074448   -0.0373812
     0.0  0.0         0.0         1.0            0.0120626   -0.0486212
     0.0  0.0         0.0         0.0           -0.0186309   -0.0454317
     0.0  0.0         0.0         0.0        …   0.900573     0.00806292
     0.0  0.0         0.0         0.0            1.0          0.00246611
     0.0  0.0         0.0         0.0            0.0          1.0




```julia
# SNP information of this block
@show df;
```

    df = 8×5 DataFrame
     Row │ AF         chr      pos    ref     alt
         │ Float64    String3  Int64  String  String
    ─────┼───────────────────────────────────────────
       1 │ 0.0259568  1        57292  C       T
       2 │ 0.992971   1        57952  A       C
       3 │ 0.0586652  1        58176  G       A
       4 │ 0.36321    1        58866  C       G
       5 │ 0.0921872  1        59040  T       C
       6 │ 0.010808   1        59108  G       A
       7 │ 0.0112479  1        59121  G       T
       8 │ 0.0370268  1        59193  T       G


When importing blocks, one can filter for minimum minor allele frequency


```julia
# keep SNPs with MAF > 0.05
chr = 1
start_pos = 57292
end_pos = 59193
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.05)
sigma
```




    3×3 Matrix{Float64}:
     1.0  -0.0125702  -0.0323548
     0.0   1.0        -0.0178985
     0.0   0.0         1.0




```julia
@show df;
```

    df = 3×5 DataFrame
     Row │ AF         chr      pos    ref     alt
         │ Float64    String3  Int64  String  String
    ─────┼───────────────────────────────────────────
       1 │ 0.0586652  1        58176  G       A
       2 │ 0.36321    1        58866  C       G
       3 │ 0.0921872  1        59040  T       C


One can also provide a list of SNP positions, and we will only keep SNPs that have those position which also pass the `min_maf` filter


```julia
# keep SNPs with MAF > 0.05 and only include a list of SNPs with known positions
chr = 1
start_pos = 57292
end_pos = 59193
snps_to_keep = [58176, 58866]
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.05, snps_to_keep=snps_to_keep)
sigma
```




    2×2 Matrix{Float64}:
     1.0  -0.0125702
     0.0   1.0




```julia
@show df;
```

    df = 2×5 DataFrame
     Row │ AF         chr      pos    ref     alt
         │ Float64    String3  Int64  String  String
    ─────┼───────────────────────────────────────────
       1 │ 0.0586652  1        58176  G       A
       2 │ 0.36321    1        58866  C       G


## SNP information

The `HailBlockMatrix` struct contains the SNP information for each row/column of the LD matrix. 


```julia
@show bm.info[1:10, :]; # print first 10 columns
```

    bm.info[1:10, :] = 10×5 DataFrame
     Row │ AF          chr      pos    ref     alt
         │ Float64     String3  Int64  String  String
    ─────┼────────────────────────────────────────────
       1 │ 0.409221    1        10146  AC      A
       2 │ 0.00810811  1        10151  TA      T
       3 │ 0.144444    1        10177  A       C
       4 │ 0.0232558   1        10178  CCTAA   C
       5 │ 0.00617284  1        10181  A       T
       6 │ 0.0326087   1        10250  A       C
       7 │ 0.0595238   1        10257  A       C
       8 │ 0.198276    1        10327  T       C
       9 │ 0.05        1        10329  AC      A
      10 │ 0.00720461  1        10333  CT      C


Alternatively, one can import this dataframe directly with the following syntax


```julia
df = read_variant_index_tables(ht_file)
```
