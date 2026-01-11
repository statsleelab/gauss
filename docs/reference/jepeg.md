# Joint effect on phenotype of eQTLs/functional SNPs associated with a gene

Joint effect on phenotype of eQTLs/functional SNPs associated with a
gene

## Usage

``` r
jepeg(
  study_pop,
  input_file,
  annotation_file,
  reference_index_file,
  reference_data_file,
  reference_pop_desc_file,
  af1_cutoff = NULL
)
```

## Arguments

- study_pop:

  study population group

- input_file:

  file name of GWAS summary statistics data containing rsid, chr, bp,
  a1, a2 and z

- reference_index_file:

  file name of reference panel index data

- reference_data_file:

  file name of reference panel data

- reference_pop_desc_file:

  file name of reference panel population description data

- af1_cutoff:

  cutoff of reference allele, a1, frequency

- annotation:

  file name of the SNP annotation data set

## Value

R dataframe containing geneid, chisq, df, jepeg_pval, num_snp,
top_categ, top_categ_pval, top_snp, and top_snp_pval
