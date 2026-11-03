# PGC2 Schizophrenia Allele Frequencies, Chromosome 22

Study reference-allele frequencies for the 13,362 chromosome 22 variants
of the Psychiatric Genomics Consortium's second schizophrenia GWAS wave
(PGC SCZ2) that were measured on the Illumina 1M array.

## Usage

``` r
data(PGC2_Chr22_ilmn1M_AF1)
```

## Format

A `data.frame` with 13,362 rows and 6 columns:

- rsid:

  dbSNP rs identifier.

- chr:

  Chromosome; `22` throughout.

- bp:

  Base pair position.

- a1:

  Reference allele.

- a2:

  Alternative allele.

- af1:

  Frequency of `a1` in the study cohort.

## Details

This is the example input for
[`afmix`](https://github.com/statsleelab/gauss/reference/afmix.md),
which estimates a cohort's ancestry proportions by matching its allele
frequencies against the reference panel. It covers the same variants, in
the same order, as
[`PGC2_Chr22_ilmn1M_Z`](https://github.com/statsleelab/gauss/reference/PGC2_Chr22_ilmn1M_Z.md),
and the proportions estimated from it are shipped as
[`PGC2_SCZ_ANC_Prop`](https://github.com/statsleelab/gauss/reference/PGC2_SCZ_ANC_Prop.md).

## References

Schizophrenia Working Group of the Psychiatric Genomics Consortium
(2014). Biological insights from 108 schizophrenia-associated genetic
loci. *Nature* 511(7510):421-427.
[doi:10.1038/nature13595](https://doi.org/10.1038/nature13595)

## Examples

``` r
if (FALSE) { # \dontrun{
data(PGC2_Chr22_ilmn1M_AF1)
str(PGC2_Chr22_ilmn1M_AF1)
} # }
```
