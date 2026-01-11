# Prepare datasets for Recessive Imputation

Prepare datasets for Recessive Imputation

## Usage

``` r
prep_recessive_impute(
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

A List containing: - snplist: A data frame of SNPs in the prediction
window with columns rsid, chr, bp, a1, a2, af1mix, z, and type, - zvec:
A numeric vector of Z-scores for measured SNPs in the extended window, -
cormat: A correlation matrix (additive-coded) among measured SNPs in the
extended window, - cormat_add: A correlation matrix between
additive-coded predicted SNPs and additive-coded measured SNPs, -
cormat_dom: A correlation matrix between dominant-coded predicted SNPs
and additive-coded measured SNPs, - cormat_rec: A correlation matrix
between recessive-coded predicted SNPs and additive-coded measured SNPs.
