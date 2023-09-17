
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GAUSS

<!-- badges: start -->
<!-- badges: end -->

GAUSS (Genome Analysis Using Summary Statistics) is a versatile R
package designed for the analysis of multi-ethnic GWAS summary
statistics. It leverages a large multi-ethnic reference panel consisting
of 32,953 genomes ([33KG
Panel](https://statsleelab.github.io/gauss/articles/ref_33KG.html)),
including 20,281 Europeans, 10,800 East Asians, 522 South Asians, 817
Africans, and 533 Native Americans. The GAUSS package offers a robust
set of tools that enable researchers to: i) accurately estimate the
ancestry proportions within a multi-ethnic GWAS (`afmix()` and
`zimx()`), ii) calculate ancestry-informed linkage disequilibrium (LD)
metrics (`computeLD()`), iii) impute Z-scores for genetic variants not
originally measured (`dist()` and `distmix()`), iv) conduct in-depth
transcriptome-wide association studies (`jepeg()` and `jepegmix()`), and
v) apply corrections for the ‘Winner’s Curse’ biases (`fiqt()`).

## Installing GSL

The GAUSS package depends on the GNU scientifiy library (GSL), so the
users should install the gsl labrary first and then install the gauss
pacakge. You can find the detailed explanation for this in this
[link](https://statsleelab.github.io/gauss/articles/gsl_installation.html):

## GAUSS Installation

You can install the development version of gauss from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("statsleelab/gauss")
```

## Suggested Workflow

![GAUSS Flowchart](GAUSS_Flowchart_bioinfo.png)

GAUSS functions utilize five core input datasets: i) Association Z-score
Data (highlighted in blue), consisting of six columns of variables: SNP
ID (rsid), chromosome number (chr), base pair position (bp), reference
allele (a1), alternative allele (a2), and Z-score (z); ii) Allele
Frequency Data (colored in red), encapsulating rsid, chr, bp, a1, a2,
and reference allele frequency (af1); iii) Ethnic Mixing Proportion Data
(colored in green), comprising population name abbreviations and
corresponding mixing proportions; iv) SNP Annotation Data (colored in
orange), offering functional annotation information of SNPs; v)
Reference Panel Data (colored in yellow), e.g., [33KG
Panel](https://statsleelab.github.io/gauss/articles/ref_33KG.html)
offering genotype information for a significant pool of 32,953
individuals from 29 diverse ethnic groups.

Initially, users should determine whether their input summary statistics
data originate from mixed ethnicity cohorts. If they do not (suggesting
the data comes from a homogeneous cohort), users can use a variety of
analysis functions specific to homogeneous cohorts, including
`computeLD()` for LD estimation, `dist()` for imputing Z-scores of
missing SNPs, and `jepeg()` for conducting trascriptome-wide associaton
study (TWAS). These functions primarily require Association Z-score Data
and Reference Panel Data.

Conversely, if the data derive from mixed ethnicity cohorts and the
ethnic mixing proportions are known, GAUSS provides a range of tools for
multi-ethnic cohort analysis. These encompass `computeLD()` for LD
estimation, `distmix()` for Z-score imputation, and `jepegmix()` for
TWAS. These tools require Association Z-score Data, Reference Panel
Data, and Mixing Proportion Data. If ethnic mixing proportions are not
known, they should be estimated using either the `afmix()` function (if
Allele Frequency Data is available) or the `zmix()` function (if it is
not). Upon determining the ethnic mixing proportions, users can apply
the genome analysis functions designed for multi-ethnic cohorts.

Small colored squares on arrows pointing to GAUSS functions denote the
required input data sets for each function. For instance, an arrow
pointing to the afmix() function is marked by two colored squares (red
and yellow), indicating that this function utilizes both Allele
Frequency Data and Reference Panel Data as inputs.
