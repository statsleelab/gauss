# PGC2 Schizophrenia GWAS Ancestry Proportion

Estimated ancestry proportions of the PGC SCZ2 cohort over the 33KG
reference panel, in the shape every mixture analysis expects for
`pop_wgt_df`: column 1 the population abbreviation, column 2 the mixing
proportion.

## Usage

``` r
data(PGC2_SCZ_ANC_Prop)
```

## Format

A `data.frame` with 21 rows and 2 columns:

- pop:

  33KG population abbreviation.

- wgt:

  Mixing proportion; non-negative, summing to 1.

## Details

Estimated with
[`afmix`](https://github.com/statsleelab/gauss/reference/afmix.md) from
the study allele frequencies in `data/PGC2_Chr22_ilmn1M_AF1.txt`; 21 of
the panel's 29 populations carry a positive weight. The proportions are
rescaled to sum to 1, which the mixture covariance assumes;
[`afmix()`](https://github.com/statsleelab/gauss/reference/afmix.md)
itself does not rescale its output, so the values as first published
summed to 1.061. Rescaling changes only the scale, not any population's
share relative to the others. See `data-raw/make_PGC2_SCZ_ANC_Prop.R`.

## References

Schizophrenia Working Group of Psychiatric Genomics Consortium (2014).
Biological insights from 108 schizophrenia associated genetic loci.
Nature [Link to article](https://www.nature.com/articles/nature13595)

## Examples

``` r
if (FALSE) { # \dontrun{
data(PGC2_SCZ_ANC_Prop)
str(PGC2_SCZ_ANC_Prop)
} # }
```
