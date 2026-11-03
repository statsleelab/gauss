# PGC2 Schizophrenia Z-scores, 3 Mb Window on Chromosome 10

Association Z-scores for 721 variants in a 2.99 Mb window on chromosome
10 (103,008,802-105,997,104), from the Psychiatric Genomics Consortium's
second schizophrenia GWAS wave (PGC SCZ2).

## Usage

``` r
data(PGC2_3Mb)
```

## Format

A `data.frame` with 721 rows and 6 columns:

- rsid:

  dbSNP rs identifier.

- chr:

  Chromosome; `10` throughout.

- bp:

  Base pair position.

- a1:

  Reference allele.

- a2:

  Alternative allele.

- z:

  Association Z-score.

## Details

The small window makes this the example input for the functions whose
cost grows with the square of the number of variants in the region:
[`computeLD`](https://github.com/statsleelab/gauss/reference/computeLD.md)
and [`dist`](https://github.com/statsleelab/gauss/reference/dist.md).

## References

Schizophrenia Working Group of the Psychiatric Genomics Consortium
(2014). Biological insights from 108 schizophrenia-associated genetic
loci. *Nature* 511(7510):421-427.
[doi:10.1038/nature13595](https://doi.org/10.1038/nature13595)

## Examples

``` r
if (FALSE) { # \dontrun{
data(PGC2_3Mb)
str(PGC2_3Mb)
} # }
```
