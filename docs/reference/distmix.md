# Direct imputation of summary statistics for unmeasured SNPs from mixed ethnicity cohorts

Direct imputation of summary statistics for unmeasured SNPs from mixed
ethnicity cohorts

## Usage

``` r
distmix(
  chr,
  start_bp,
  end_bp,
  wing_size,
  pop_wgt_df,
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

- wing_size:

  the size of the area flanking the left and right of the prediction
  window

- pop_wgt_df:

  R data frame with the population abbreviation in the first column and
  its mixing proportion in the second, read positionally. Proportions
  must be non-negative; they are rescaled to sum to 1, which the mixture
  covariance assumes. See `data(PGC2_SCZ_ANC_Prop)` for the expected
  shape. Note that
  [`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
  returns three columns, so its first column must be dropped first.

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

R data frame containing rsid, chr, bp, a1, a2, af1mix, z, pval, info,
type
