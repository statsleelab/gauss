# Estimate ancestry proportions from Z-scores

Runs the zmix preparation and constrained quadratic programming to
estimate ancestry proportions at the population or superpopulation
level.

## Usage

``` r
zmix(
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  percentile = 0.9,
  interval = 10,
  level = c("population", "superpopulation")
)
```

## Arguments

- input_file:

  File name of input data containing rsid, chr, bp, a1, a2, and z.

- reference_index_file:

  File name of reference panel index data.

- reference_data_file:

  File name of reference panel data.

- reference_pop_desc_file:

  File name of reference panel population description data.

- percentile:

  Percentile cutoff for normalized variance.

- interval:

  Interval for selecting SNPs.

- level:

  Resolution of ancestry proportions: "population" or "superpopulation".

## Value

A data.frame of ancestry proportions at the requested level.
