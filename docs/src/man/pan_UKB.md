
# Pan-UKB data

+ Link: https://pan-dev.ukbb.broadinstitute.org/docs/hail-format/index.html#extracting-a-subset-of-ld-matrix
+ European samples are the largest, totaling 14.1T of data


```julia
using EasyLD
using CSV
using DataFrames
using Statistics
using LinearAlgebra
```

## Downloading

First check how many files there are:


```julia
get_ukb_filenames("EUR", join=false)
```




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



Try downloading the first 10 files (~1GB)


```julia
population = "EUR"
outdir = "/Users/biona001/.julia/dev/EasyLD/data"
download_ukb_LD_matrices(population, outdir, start_from=1, num_files=10)
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:00:49[39m


+ The result would be saved into outdir directory with name `UKBB.EUR.ldadj.bm` which itself is a directory.

+ For multithreaded downloads, one can use `start_from` and `num_files` keywords to control how many matrices to download (not specifying will download all matrices). A progress meter is automatically displayed.

Finally, download variant index hail table


```julia
population = "EUR"
outdir = "/Users/biona001/.julia/dev/EasyLD/data"
download_ukb_variant_index_tables(population, outdir)
```

The result will be stored as `UKBB.$population.ldadj.variant.ht` in `outdir`

## Reading LD panel with Matrix interface

The first time `hail_block_matrix` gets called, we do some pre-processing to the variant index files so subsequent calls will be faster. 


```julia
bm_file = "/Users/biona001/.julia/dev/EasyLD/data/UKBB.EUR.ldadj.bm"
ht_file = "/Users/biona001/.julia/dev/EasyLD/data/UKBB.EUR.ldadj.variant.ht"
@time bm = hail_block_matrix(bm_file, ht_file); # the ';' avoids displaying a few entries of bm, which takes ~0.1 seconds per entry
```

    130.413312 seconds (978.68 M allocations: 23.643 GiB, 50.65% gc time, 4.97% compilation time)


Check size of matrix


```julia
size(bm)
```




    (23960350, 23960350)



Read first 10000 by 10000 block into memory


```julia
Sigma = bm[1:10000, 1:10000]
```

    [Stage 1:=================================================>         (5 + 1) / 6]




    10000×10000 Matrix{Float64}:
     0.99996  0.00098559  -0.000830543  …  -0.000520665   0.00422904
     0.0      0.999942    -0.00104453      -0.00494589    0.00050959
     0.0      0.0          0.999963        -0.00129536    0.00216722
     0.0      0.0          0.0             -0.00663755    0.000142521
     0.0      0.0          0.0              0.0104731     0.000936687
     0.0      0.0          0.0          …  -0.00202249   -0.000497824
     0.0      0.0          0.0             -0.00550819    0.00283421
     0.0      0.0          0.0             -0.00145102   -0.000761935
     0.0      0.0          0.0              0.00287458   -0.0010154
     0.0      0.0          0.0             -0.000450895  -0.00130722
     0.0      0.0          0.0          …  -0.00101125    0.000451512
     0.0      0.0          0.0              0.00590455   -0.00110266
     0.0      0.0          0.0             -0.00274252   -0.000817142
     ⋮                                  ⋱                
     0.0      0.0          0.0             -0.109458     -0.00866031
     0.0      0.0          0.0              0.0158602     0.0203563
     0.0      0.0          0.0          …   0.00745596    0.0184919
     0.0      0.0          0.0              0.983906     -0.0122066
     0.0      0.0          0.0              0.0155851     0.0197608
     0.0      0.0          0.0              0.965104     -0.0159355
     0.0      0.0          0.0             -0.129534      0.088277
     0.0      0.0          0.0          …   0.622303      0.059004
     0.0      0.0          0.0              0.965347     -0.0160028
     0.0      0.0          0.0              0.0925609     0.999286
     0.0      0.0          0.0              0.999972      0.0926374
     0.0      0.0          0.0              0.0           0.999974



Check if the given block is PSD by computing its eigenvalues


```julia
eigvals(Symmetric(Sigma)) # Symmetric uses upper triangular portion of data
```




    10000-element Vector{Float64}:
      -4.496303144059796e-15
      -1.762987426134658e-15
      -1.040030990442209e-15
      -4.2089306816541705e-16
      -3.960324780237165e-16
      -3.9060196040323796e-16
      -2.985734642577166e-16
      -2.8601497287222336e-16
      -2.5553650845256305e-16
      -1.998130772222011e-16
      -1.8449274904356797e-16
      -1.7028015609869147e-16
      -1.3192630645503314e-16
       ⋮
     116.57001725005479
     124.37732317876345
     125.11477060743002
     135.80425399688727
     142.93300594996916
     149.24000792252738
     161.47927858350278
     181.82233513932786
     228.3300216753156
     235.961635959381
     280.6540947869515
     322.07381500507836



## Read in a block with `get_block`

One can also extract a block by specifying the chromosome and starting/ending basepair


```julia
chr = 1
start_pos = 11063
end_pos = 91588
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.0)
sigma
```




    9×9 Matrix{Float64}:
     0.99996  0.00098559  -0.000830543  …  -0.00045888    0.000221516
     0.0      0.999942    -0.00104453      -0.000332668  -0.000439517
     0.0      0.0          0.999963        -0.000551892  -0.000786793
     0.0      0.0          0.0             -0.000503539  -0.000712834
     0.0      0.0          0.0             -0.000294435   0.00575015
     0.0      0.0          0.0          …  -0.000112519  -0.00015419
     0.0      0.0          0.0             -0.000264108   0.00413586
     0.0      0.0          0.0              0.999978     -0.000225057
     0.0      0.0          0.0              0.0           0.99995




```julia
# SNP information of this block
@show df;
```

    df = 9×6 DataFrame
     Row │ rsid         AF          chr      pos    ref     alt
         │ String       Float64     String3  Int64  String  String
    ─────┼─────────────────────────────────────────────────────────
       1 │ rs561109771  4.7982e-5   1        11063  T       G
       2 │ rs562993331  0.00027798  1        13259  G       A
       3 │ rs578081284  0.00083096  1        17641  G       A
       4 │ rs576081345  0.00065859  1        57222  T       C
       5 │ rs570371753  0.00024023  1        58396  T       C
       6 │ rs561430336  2.7728e-5   1        63668  G       A
       7 │ rs2531267    0.00018542  1        69569  T       C
       8 │ rs557418932  8.1599e-5   1        79192  T       G
       9 │ rs554639997  0.00015868  1        91588  G       A


When importing blocks, one can filter for minimum minor allele frequency


```julia
# keep SNPs with MAF > 0.0001
chr = 1
start_pos = 11063
end_pos = 91588
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.0001)
sigma
```




    6×6 Matrix{Float64}:
     0.999942  -0.00104453  -0.000939341  …  -0.000502128  -0.000439517
     0.0        0.999963    -0.00164293      -0.000874333  -0.000786793
     0.0        0.0          0.999972        -0.00079952   -0.000712834
     0.0        0.0          0.0             -0.00044672    0.00575015
     0.0        0.0          0.0              0.999986      0.00413586
     0.0        0.0          0.0          …   0.0           0.99995




```julia
@show df;
```

    df = 6×6 DataFrame
     Row │ rsid         AF          chr      pos    ref     alt
         │ String       Float64     String3  Int64  String  String
    ─────┼─────────────────────────────────────────────────────────
       1 │ rs562993331  0.00027798  1        13259  G       A
       2 │ rs578081284  0.00083096  1        17641  G       A
       3 │ rs576081345  0.00065859  1        57222  T       C
       4 │ rs570371753  0.00024023  1        58396  T       C
       5 │ rs2531267    0.00018542  1        69569  T       C
       6 │ rs554639997  0.00015868  1        91588  G       A


One can also provide a list of SNP positions, and we will only keep SNPs that have those position which also pass the `min_maf` filter


```julia
# keep SNPs with MAF > 0.0001 and only include a list of SNPs with known positions
chr = 1
start_pos = 11063
end_pos = 91588
snps_to_keep = [13259, 58396, 91588]
sigma, df = get_block(bm, chr, start_pos, end_pos; min_maf=0.0001, snps_to_keep=snps_to_keep)
sigma
```

    [Stage 2:>                                                          (0 + 1) / 1]




    3×3 Matrix{Float64}:
     0.999942  -0.00056507  -0.000439517
     0.0        0.999982     0.00575015
     0.0        0.0          0.99995




```julia
@show df;
```

    df = 3×6 DataFrame
     Row │ rsid         AF          chr      pos    ref     alt
         │ String       Float64     String3  Int64  String  String
    ─────┼─────────────────────────────────────────────────────────
       1 │ rs562993331  0.00027798  1        13259  G       A
       2 │ rs570371753  0.00024023  1        58396  T       C
       3 │ rs554639997  0.00015868  1        91588  G       A


## SNP information

SNP information can be accessed from the `bm.info` field, where `AF` corresponds to alternate allele frequency. 


```julia
@show bm.info[1:10, :]; # first 10 rows
```

    bm.info[1:10, :] = 10×6 DataFrame
     Row │ rsid         AF          chr      pos     ref     alt
         │ String       Float64     String3  Int64   String  String
    ─────┼──────────────────────────────────────────────────────────
       1 │ rs561109771  4.7982e-5   1         11063  T       G
       2 │ rs562993331  0.00027798  1         13259  G       A
       3 │ rs578081284  0.00083096  1         17641  G       A
       4 │ rs576081345  0.00065859  1         57222  T       C
       5 │ rs570371753  0.00024023  1         58396  T       C
       6 │ rs561430336  2.7728e-5   1         63668  G       A
       7 │ rs2531267    0.00018542  1         69569  T       C
       8 │ rs557418932  8.1599e-5   1         79192  T       G
       9 │ rs554639997  0.00015868  1         91588  G       A
      10 │ rs575442534  0.00066441  1        533573  G       A


Alternatively, one can read SNP information by giving the variant index folder


```julia
df = read_variant_index_tables(ht_file)
```
