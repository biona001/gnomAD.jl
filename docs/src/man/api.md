
# API

Here is a list of available function calls. A detailed description can be found below. 
```@index
Pages = ["api.md"]
```

## Exported functions
```@docs
  download_gnomad_LD_matrices
  download_ukb_LD_matrices
  download_gnomad_variant_index_tables
  download_ukb_variant_index_tables
  get_gnomad_filenames
  get_ukb_filenames 
  hail_block_matrix
  read_variant_index_tables
  get_block
```

## Internal helper functions and structs

```@docs
  EasyLD.HailBlockMatrix
  EasyLD._extract_alternate_allele_freq
  EasyLD._extract_ref_alt_alleles
  EasyLD._extract_genome_build
  EasyLD._read_metadata
```
