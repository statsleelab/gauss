# Calculate population weights using association Z-scores

This function mirrors prep_zmix5 but computes correlations at the
superpopulation level.

## Usage

``` r
prep_zmix5_sup(
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  percentile = NULL,
  interval = NULL
)
```

## Arguments

- input_file:

  file name of input data containing rsid, chr, bp, a1, a2, and z

- reference_index_file:

  file name of reference panel index data

- reference_data_file:

  file name of reference panel data

- reference_pop_desc_file:

  file name of reference panel population description data

- percentile:

  percentile cutoff for normalized variance

- interval:

  stepping distance within the SNP vector for selecting the first SNP of
  each pair.

## Value

R data frame containing superpopulation IDs and weights
