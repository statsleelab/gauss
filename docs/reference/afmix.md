# Calculate population weights using RAF

Calculate population weights using RAF

## Usage

``` r
afmix(
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  interval = NULL
)
```

## Arguments

- input_file:

  file name of input data containing rsid, chr, bp, a1, a2, and af1

- reference_index_file:

  file name of reference panel index data

- reference_data_file:

  file name of reference panel data

- reference_pop_desc_file:

  file name of reference panel population description data

- interval:

  number of non-overlapping SNP sets used in calculating population
  weights

## Value

R data frame containing population IDs and weights
