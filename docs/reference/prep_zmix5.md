# Calculate population weights using association Z-scores

This function first finds ancestry-informed SNPs by computing a
normalized variance, ratio = var(af1)/(mean(af1)\*(1-mean(af1))) and
choose top 1 ratio and use these SNPs in the zmix_prep procedure. It
uses an interval parameter to determine how SNPs are selected for
analysis. If interval is not provided, it defaults to 1, meaning every
SNP is considered. SNP Pairing: Pairs each SNP with every other SNP in
the subvector (snp_subvec) If Interval=1000, it pairs SNP_0 with
SNP_1000, SNP_0 with SNP_2000 ... and SNP_1000 with SNP_2000, SNP_1000
with SNP_3000 ...

## Usage

``` r
prep_zmix5(
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

R data frame containing population IDs and weights
