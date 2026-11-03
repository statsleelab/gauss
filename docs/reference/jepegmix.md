# JEPEGMIX: gene-level joint analysis of functional SNPs in cosmopolitan cohorts

JEPEGMIX: gene-level joint analysis of functional SNPs in cosmopolitan
cohorts

## Usage

``` r
jepegmix(
  pop_wgt_df,
  input_file,
  annotation_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  af1_cutoff = NULL
)
```

## Arguments

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
  a1, a2 and z

- annotation_file:

  file name of the SNP annotation data: a whitespace-delimited table
  with a header line and the columns `rsid`, `chr`, `bp`, `a1`, `a2`,
  `geneid`, `categ`, `wgt`. Variants are matched to the study data by
  `(chr, bp, a1, a2)` in either allele order, never by `rsid`, so `rsid`
  may be given as `.`. `categ` must be one of `PROTEIN`, `TFBS`,
  `WTH_HAIR`, `WTH_TARGET`, `CIS_EQTL` or `TRANS_EQTL`, reported in the
  result as `PFS`, `TFB`, `STR`, `TAR`, `CIS` and `TRN` respectively;
  rows carrying any other value are skipped and counted in a warning.
  `wgt` is the variant's weight within its category.

- reference_index_file:

  file name of reference panel index data

- reference_data_file:

  file name of reference panel data

- reference_pop_desc_file:

  file name of reference panel population description data

- af1_cutoff:

  cutoff of reference allele, a1, frequency

## Value

R dataframe containing geneid, chisq, df, jepeg_pval, num_snp,
top_categ, top_categ_pval, top_snp, and top_snp_pval
