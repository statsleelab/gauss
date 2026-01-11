# Calculate population weights using association Z-scores

Interval and Offset Logic: Instead of considering all pairwise
combinations, zmix2 selects SNP pairs based on a fixed offset. For
example, if offset is 3 and and interval is 1000, it pairs SNP_0 with
SNP_3, SNP_1000 with SNP_1003 ... Reduced Pairing Scope: This approach
reduces the total number of SNP pairs considered compared to zmix,
focusing on pairs separated by a specific distance (offset) in the SNP
vector.

## Usage

``` r
prep_zmix2(
  input_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  interval = NULL,
  offset = NULL
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

- offset:

  distance between the two SNPs in a pair

## Value

R data frame containing population IDs and weights
