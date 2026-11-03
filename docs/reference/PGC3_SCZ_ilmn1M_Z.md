# PGC3 Schizophrenia Z-scores, Genome-wide

Association Z-scores for 791,425 variants across chromosomes 1-22, from
the Psychiatric Genomics Consortium's third schizophrenia GWAS wave (PGC
SCZ3), restricted to variants measured on the Illumina 1M array.

## Usage

``` r
data(PGC3_SCZ_ilmn1M_Z)
```

## Format

A `data.frame` with 791,425 rows and 6 columns:

- rsid:

  dbSNP rs identifier.

- chr:

  Chromosome, `1`-`22`.

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
[`fiqt`](https://github.com/statsleelab/gauss/reference/fiqt.md), whose
FDR inverse quantile transformation corrects Winner's Curse bias and
therefore needs a genome-wide set of Z-scores rather than one region.

## References

Trubetskoy V, Pardinas AF, Qi T, et al. (2022). Mapping genomic loci
implicates genes and synaptic biology in schizophrenia. *Nature*
604(7906):502-508.
[doi:10.1038/s41586-022-04434-5](https://doi.org/10.1038/s41586-022-04434-5)

## Examples

``` r
if (FALSE) { # \dontrun{
data(PGC3_SCZ_ilmn1M_Z)
str(PGC3_SCZ_ilmn1M_Z)
} # }
```
