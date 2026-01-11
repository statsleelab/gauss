# Estimating ancestry proportions

## Introduction

As genome-wide association studies (GWAS) continue to expand in scope
and diversity, the precise quantification of a study’s ancestral makeup
has become increasingly important. The
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function of the GAUSS package, provides users with a robust method for
accurately estimating these ancestry proportions within multi-ethnic
GWAS, utilizing solely summary statistics data. This vignette offers a
step-by-step guide on how to effectively employ the
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function.

## Load necessary packages

``` r

library(gauss)
library(dplyr)
library(data.table)
library(kableExtra)
```

## Overview of afmix()

The [`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function requires four specific arguments:

- `input_file`: File name of the allele frequency data
- `reference_index_file`: File name of reference panel index data
- `reference_data_file`: File name of reference panel data
- `reference_pop_desc_file`: File name of reference panel population
  description data

Upon execution,
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function return a data.frame comprising three columns:

- `sup.pop`: super population name
- `pop`: population name
- `wgt`: population proportion

## Preparing the Input Data

### Allele Frequency Data

The allele frequency data file should be formatted as a space-delimited
text file, containing six specific columns: `rsid` (SNP ID), `chr`
(chromosome number), `bp` (base pair position), `a1` (reference allele),
`a2` (alternative allele), and `af1` (reference allele frequency). In
this vignette, we will be using allele frequency data of Chromosome 22
from the Psychiatric Genomic Consortium’s Phase 2 Schizophrenia (PGC
SCZ2) GWAS.

Below, we specify the path to the allele frequency data file.

``` r

# Path to the input file
input_file <- "../data/PGC2_Chr22_ilmn1M_AF1.txt"

# Input file should include six columns (rsid, chr, bp, a1, a2, and af1)
input.data <- fread(input_file, header = TRUE)
head(input.data)
#>         rsid   chr       bp     a1     a2       af1
#>       <char> <int>    <int> <char> <char>     <num>
#> 1: rs1000427    22 36890105      A      G 0.1159800
#> 2: rs1000470    22 24026845      A      C 0.1240050
#> 3: rs1000539    22 20202729      A      G 0.2660650
#> 4:   rs10009    22 22051709      A      G 0.6319500
#> 5: rs1001022    22 26403488      T      C 0.9570250
#> 6: rs1001213    22 34131736      A      G 0.0667805
```

### Reference Panel Data Files

Next, we assign the paths to reference panel data files. In this
example, we use [the 33KG reference
panel](https://statsleelab.github.io/gauss/articles/ref_33KG.html).

``` r

# Paths to the reference files (replace these with your actual paths)
reference_index_file <-"../ref/33KG/33kg_index.gz"
reference_data_file <- "../ref/33KG/33kg_geno.gz"
reference_pop_desc_file<-"../ref/33KG/33kg_pop_desc.txt"
```

## Running afmix()

With the necessary arguments and data files in place, we are ready to
run the
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function to compute ancestry proportion of `PGC SCZ2 GWAS`.

``` r

wgt.df <- afmix(input_file=input_file,
                reference_index_file = reference_index_file,
                reference_data_file = reference_data_file,
                reference_pop_desc_file = reference_pop_desc_file)
```

## Results: Estimated Ancestry Proportions

Here, we display the estimated ancestry proportions in a HTML table.

``` r

wgt.df %>% kable("html")
```

| sup.pop | pop |   wgt |
|:--------|:----|------:|
| AFR     | ACB | 0.006 |
| AFR     | ASW | 0.036 |
| SAS     | BEB | 0.005 |
| ASN     | CCE | 0.008 |
| ASN     | CCS | 0.004 |
| ASN     | CDX | 0.018 |
| EUR     | CEU | 0.165 |
| AMR     | CLM | 0.025 |
| ASN     | CNE | 0.003 |
| ASN     | CSE | 0.012 |
| EUR     | FIN | 0.138 |
| EUR     | GBR | 0.165 |
| SAS     | GIH | 0.006 |
| EUR     | IBS | 0.099 |
| ASN     | JPT | 0.011 |
| ASN     | KHV | 0.017 |
| AMR     | MXL | 0.030 |
| EUR     | ORK | 0.166 |
| SAS     | PJL | 0.016 |
| AMR     | PUR | 0.045 |
| EUR     | TSI | 0.086 |

## Summarizing Ancestry Proportions by Super Population

Here, we calculate the total proportion for each super population and
present the results in a table.

``` r

wgt.df %>%
  group_by(sup.pop) %>%
  summarise(wgt=sum(wgt), .groups="drop") %>%
  kable("html")
```

| sup.pop |   wgt |
|:--------|------:|
| AFR     | 0.042 |
| AMR     | 0.100 |
| ASN     | 0.073 |
| EUR     | 0.819 |
| SAS     | 0.027 |

## References

- Lee et al. DISTMIX: direct imputation of summary statistics for
  unmeasured SNPs from mixed ethnicity cohorts. Bioinformatics.
  <https://doi.org/10.1093/bioinformatics/btv348>
