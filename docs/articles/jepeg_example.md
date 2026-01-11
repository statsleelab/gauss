# Transcriptome wide association study with GAUSS

## Introduction

Transcriptome-Wide Association Studies (TWAS) serve as a powerful
framework for identifying the functional relationships between genetic
variations and complex traits. By integrating results from Genome-Wide
Association Studies (GWAS) with rich functional annotations, TWAS
provides a more comprehensive view of the genetic architecture of traits
and diseases. GAUSS offers a robust suite of tools to streamline and
enhance TWAS analyses. Two of its specialized functions,
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md) and
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md),
offer tailored solutions for homogeneous and multi-ethnic cohorts,
respectively. This vignette aims to be your comprehensive guide, walking
you through the hands-on application of
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md) and
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
in conducting TWAS.

## Load necessary packages

``` r

# Load necessary packages
library(gauss)
library(tidyverse)
library(data.table)
library(kableExtra)
```

## Preparing the Input Data

Both
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md) and
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
functions necessitate three core input datasets:

- Association Z-score Data
- SNP Annotation Data
- Reference Panel Data

### Association Z-score Data

This data should be in the form of a space-delimited text file including
six columns with the following names:

- `rsid`: SNP ID
- `chr`: Chromosome number
- `bp`: Base pair position
- `a1`: Reference allele
- `a2`: Alternative allele
- `z`: Association Z-score

For this example, we will use the Psychiatric Genomic Consortium’s Phase
2 Schizophrenia (PGC SCZ2) GWAS dataset.

``` r

# Path to the input file
input_file <- "../data/PGC2_Chr22_ilmn1M_Z.txt"

# Input file should include six columns (rsid, chr, bp, a1, a2, and z)
input.data <- fread(input_file, header = TRUE)
head(input.data)
#>         rsid   chr       bp     a1     a2          z
#>       <char> <int>    <int> <char> <char>      <num>
#> 1: rs1000427    22 36890105      A      G -1.4969741
#> 2: rs1000470    22 24026845      A      C  1.8407377
#> 3: rs1000539    22 20202729      A      G  2.1683772
#> 4:   rs10009    22 22051709      A      G -0.0556759
#> 5: rs1001022    22 26403488      T      C  2.9694844
#> 6: rs1001213    22 34131736      A      G -0.8957244
```

### SNP Annotation Data

The SNP Annotation Data should be placed in an accessible directory, as
shown below:

``` r

# Paths to the annotation file (replace these with your actual paths)
annotation_file <- "../data/JEPEG_SNP_Annotation.v1.0.txt"
```

### Reference Panel Data

We will be using [the 33KG reference
panel](https://statsleelab.github.io/gauss/articles/ref_33KG.html). Make
sure to replace the file paths below with those of your actual reference
panel files.

``` r

reference_index_file <-"../ref/33KG/33kg_index.gz"
reference_data_file <- "../ref/33KG/33kg_geno.gz"
reference_pop_desc_file<-"../ref/33KG/33kg_pop_desc.txt"
```

## The `jepeg()` function

The [`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md)
function enables researchers to perform TWAS by assessing the joint
effect of expression quantitative trait loci (eQTLs) and other
functional variants within a specific gene on phenotype. It requires an
SNP annotation file that provides a linkage between SNPs and their
corresponding gene-level functional characteristics. By integrating this
data, jepeg() formulates gene-level test statistics, aiming to pinpoint
genes whose expression levels might be modulated by genetic variations.”

### Arguments

The [`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md)
function requires the following arguments:

- `study_pop`: Study population group
- `input_file`: File name of the association Z-score data
- `annotation`: File name of the SNP annotation data set
- `reference_index_file`: File name of reference panel index data
- `reference_data_file`: File name of reference panel data
- `reference_pop_desc_file`: File name of reference panel population
  description data
- `af1_cutoff`: Cutoff of reference allele (a1) frequency

### Outputs

The [`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md)
function returns a data frame with following columns:

- `geneid`: Gene name.
- `chisq`: Chi-square test statistic value
- `df`: Degrees of freedom
- `jepeg_pval`: JEPEG p-value
- `num_snp`: Number of functional SNPs associated with gene
- `top_categ`: Top functional category
- `top_categ_pval`: Top category p-value
- `top_snp`: Top functional SNP ID
- `top_snp_pval`: Top functional SNP p-value

### Example Usage

In this example, we will use the
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md)
function to perform a TWAS analysis utilizing summary statistics from
the PGC SCZ2 GWAS dataset. Assuming our study cohort consists of
participants of British descent, we will specify `study_pop = "GBR"`.
This will ensure the function employs genotype data for individuals from
Britain, specifically England and Scotland, as sourced from the 33KG
reference panel.

``` r

af1_cutoff = 0.001

res <- jepeg(study_pop = "GBR",
             input_file=input_file,
             annotation_file=annotation_file,
             reference_index_file = reference_index_file,
             reference_data_file = reference_data_file,
             reference_pop_desc_file = reference_pop_desc_file)

res <- res[order(res$jepeg_pval),]
```

### `jepeg()` Results

``` r

head(res) %>% kable("html")
```

|  | geneid | chisq | df | jepeg_pval | num_snp | top_categ | top_categ_pval | top_snp | top_snp_pval |
|:---|:---|---:|---:|---:|---:|:---|---:|:---|---:|
| 85 | DPYD | 38.41841 | 1 | 0.00e+00 | 1 | TRN | 0.00e+00 | rs3788568 | 0.0e+00 |
| 72 | CXCL14 | 33.98061 | 1 | 0.00e+00 | 1 | TRN | 0.00e+00 | rs133047 | 0.0e+00 |
| 93 | EP300 | 29.29304 | 1 | 1.00e-07 | 1 | PFS | 1.00e-07 | rs20551 | 0.0e+00 |
| 328 | WBP2NL | 24.71184 | 1 | 7.00e-07 | 2 | PFS | 7.00e-07 | rs2301521 | 1.0e-07 |
| 186 | NDUFA6 | 24.39774 | 1 | 8.00e-07 | 1 | PFS | 8.00e-07 | rs1801311 | 2.0e-07 |
| 332 | ZBED4 | 19.38566 | 1 | 1.07e-05 | 1 | PFS | 1.07e-05 | rs910799 | 3.9e-06 |

## The `jepegmix()` function

For studies that involve diverse ethnic backgrounds,
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
extends the utility of
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md) to
accommodate the complexities introduced by ethnic heterogeneity.

### Arguments

The
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
function accepts most of the same arguments as
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md).
However, it replaces the `study_pop` argument with `pop_wgt_df`:

- `pop_wgt_df`: An R data frame containing population names and
  corresponding ancestry proportions.

### Outputs

The output of
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
is identical to that of
[`jepeg()`](https://github.com/statsleelab/gauss/reference/jepeg.md).

### Example Usage

Before using
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md),
you need to prepare ancestry proportion data. This data should be
structured as a data frame containing two columns:

- `pop`: Population abbreviation
- `wgt`: Corresponding proportion for that population

To estimate these ancestry proportions, the
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
function can be employed. For a step-by-step guide on this process,
refer to the [`afmix()`
vignette](https://statsleelab.github.io/gauss/articles/afmix_example.html).

Here, we load pre-generated ancestry proportion data:

``` r

# Load the ancestry proportion data
data("PGC2_SCZ_ANC_Prop") # data frame name: PGC2_SCZ_ANC_Prop
head(PGC2_SCZ_ANC_Prop)
#>   pop   wgt
#> 1 ACB 0.006
#> 2 ASW 0.036
#> 3 BEB 0.005
#> 4 CCE 0.008
#> 5 CCS 0.004
#> 6 CDX 0.018
```

Now, let’s proceed to use
[`jepegmix()`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
to conduct TWAS. Here, we assume the study cohort includes participants
from diverse ethnic backgrounds.

``` r

af1_cutoff = 0.001

res <- jepegmix(pop_wgt_df = PGC2_SCZ_ANC_Prop,
                input_file=input_file,
                annotation_file=annotation_file,
                reference_index_file = reference_index_file,
                reference_data_file = reference_data_file,
                reference_pop_desc_file = reference_pop_desc_file,
                af1_cutoff = af1_cutoff)

res <- res[order(res$jepeg_pval),]
```

### `jepegmix()` Results

``` r

head(res) %>% kable("html")
```

|  | geneid | chisq | df | jepeg_pval | num_snp | top_categ | top_categ_pval | top_snp | top_snp_pval |
|:---|:---|---:|---:|---:|---:|:---|---:|:---|---:|
| 92 | DPYD | 38.41841 | 1 | 0.00e+00 | 1 | TRN | 0.00e+00 | rs3788568 | 0.0e+00 |
| 75 | CXCL14 | 33.81352 | 1 | 0.00e+00 | 2 | TRN | 0.00e+00 | rs133047 | 0.0e+00 |
| 101 | EP300 | 29.29304 | 1 | 1.00e-07 | 1 | PFS | 1.00e-07 | rs20551 | 0.0e+00 |
| 361 | WBP2NL | 24.71140 | 1 | 7.00e-07 | 2 | PFS | 7.00e-07 | rs2301521 | 1.0e-07 |
| 206 | NDUFA6 | 24.39774 | 1 | 8.00e-07 | 1 | PFS | 8.00e-07 | rs1801311 | 2.0e-07 |
| 365 | ZBED4 | 19.38566 | 1 | 1.07e-05 | 1 | PFS | 1.07e-05 | rs910799 | 3.9e-06 |

## References

- Lee et al. JEPEG: a summary statistics based tool for gene-level joint
  testing of functional variants. Bioinformatics, Volume 31, Issue 8,
  April 2015, Pages 1176–1182,
  <https://doi.org/10.1093/bioinformatics/btu816>
- Lee et al. JEPEGMIX: gene-level joint analysis of functional SNPs in
  cosmopolitan cohorts. Bioinformatics. 2016 Jan 15; 32(2): 295–297.
  doi: 10.1093/bioinformatics/btv567
