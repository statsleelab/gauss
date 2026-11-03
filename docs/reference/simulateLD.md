# Simulate LD for user provided SNPs from mixed ethnicity cohorts

Simulate LD for user provided SNPs from mixed ethnicity cohorts

## Usage

``` r
simulateLD(
  chr,
  start_bp,
  end_bp,
  pop_wgt_df,
  sim_size,
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  af1_cutoff = NULL
)
```

## Arguments

- chr:

  chromosome number

- start_bp:

  start base pair position of prediction window

- end_bp:

  end base pair position of prediction window

- pop_wgt_df:

  R data frame with the population abbreviation in the first column and
  its mixing proportion in the second, read positionally. Proportions
  must be non-negative; they are rescaled to sum to 1, which the mixture
  covariance assumes. See `data(PGC2_SCZ_ANC_Prop)` for the expected
  shape. Note that
  [`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
  returns three columns, so its first column must be dropped first.

- sim_size:

  Number of simulated subjects

- input_file:

  file name of GWAS summary statistics data containing rsid, chr, bp,
  a1, a2, af1, and z

- reference_index_file:

  file name of reference panel index data

- reference_data_file:

  file name of reference panel data

- reference_pop_desc_file:

  file name of reference panel population description data

- af1_cutoff:

  cutoff of reference allele, a1, frequency

## Value

R list containing a data frame containing rsid, chr, bp, a1, a2 and
af1mix and a correlation matrix
