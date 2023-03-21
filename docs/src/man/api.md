
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

## Internal helpder functions
```@docs
  _extract_alternate_allele_freq
  _extract_ref_alt_alleles
  _extract_genome_build
  _read_metadata
```