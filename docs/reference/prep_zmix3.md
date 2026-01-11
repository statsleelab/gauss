# Calculate population weights using association Z-scores

Pairing Logic Based on Steps: This function pairs each SNP with a
specific number of subsequent SNPs as defined by the 'steps' parameter.
This function is a variant of ZMIX. While ZMIX2 pairs one SNP in the
snp_vec with all the remaining SNPs, ZMIX3 limits the number of SNPs to
be paired by using the 'steps' argument. For example, if 'steps' is set
to 5, each SNP in the snp_subvec is paired with its next 5 subsequent
SNPs.

## Usage

``` r
prep_zmix3(
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  interval = NULL,
  steps = NULL
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

- interval:

  stepping distance within the SNP vector for selecting the first SNP of
  each pair.

- steps:

  number of subsequent SNPs to pair with each SNP.

## Value

R data frame containing population IDs and weights
