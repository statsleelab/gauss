# Compute LD for measured SNPs from mixed ethnicity cohorts

Compute LD for measured SNPs from mixed ethnicity cohorts

## Usage

``` r
computeLD(
  chr,
  start_bp,
  end_bp,
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

- pop_wgt_df:

  R data frame containing population IDs and weights

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
