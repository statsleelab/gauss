# Calculate population weights using association Z-scores

ZMIX4 employs a more selective approach by using a fixed offset to
determine the SNP pairs. For instance, if the offset is set to 3 and the
interval is specified as 1000, the function will pair the first SNP
(SNP_0) with the fourth one (SNP_3), then SNP_1000 with SNP_1003, and so
forth. In subsequent iterations, this pattern continues in a similar
manner - SNP_1 gets paired with SNP_4, SNP_1001 with SNP_1004, and so
on, effectively creating pairs at specific intervals with a consistent
offset.

## Usage

``` r
prep_zmix4(
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
