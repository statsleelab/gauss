# PGC2 Schizophrenia Z-scores, Chromosome 22

Association Z-scores for the 13,362 chromosome 22 variants of the
Psychiatric Genomics Consortium's second schizophrenia GWAS wave (PGC
SCZ2) that were measured on the Illumina 1M array, spanning
16,054,667-51,213,613.

## Usage

``` r
data(PGC2_Chr22_ilmn1M_Z)
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

- z:

  Association Z-score.

## Details

This is the example input for
[`jepeg`](https://github.com/statsleelab/gauss/reference/jepeg.md). It
covers the same variants, in the same order, as
[`PGC2_Chr22_ilmn1M_AF1`](https://github.com/statsleelab/gauss/reference/PGC2_Chr22_ilmn1M_AF1.md).

## References

Schizophrenia Working Group of the Psychiatric Genomics Consortium
(2014). Biological insights from 108 schizophrenia-associated genetic
loci. *Nature* 511(7510):421-427.
[doi:10.1038/nature13595](https://doi.org/10.1038/nature13595)

## Examples

``` r
if (FALSE) { # \dontrun{
data(PGC2_Chr22_ilmn1M_Z)
str(PGC2_Chr22_ilmn1M_Z)
} # }
```
