# JEPEG Functional SNP Annotation, version 1.0

The SNP-to-gene functional annotation that
[`jepeg`](https://github.com/statsleelab/gauss/reference/jepeg.md) and
[`jepegmix`](https://github.com/statsleelab/gauss/reference/jepegmix.md)
read through their `annotation_file` argument. Each row ties one variant
to one gene under one functional category and carries the weight that
category's test applies to it.

## Usage

``` r
data(JEPEG_SNP_Annotation.v1.0)
```

## Format

A `data.frame` with 509,570 rows and 8 columns:

- rsid:

  dbSNP rs identifier, or `.` where none is recorded.

- chr:

  Chromosome, `1`-`22`, `X` or `Y`.

- bp:

  Base pair position.

- a1:

  Reference allele.

- a2:

  Alternative allele.

- geneid:

  Gene symbol.

- categ:

  Functional category: `PROTEIN` (260,518 rows), `TRANS_EQTL` (94,925),
  `CIS_EQTL` (75,118), `TFBS` (58,695), `WTH_TARGET` (19,660, within a
  microRNA target site) or `WTH_HAIR` (654, within a microRNA hairpin).

- wgt:

  Signed weight for this variant in this category.

## Details

A variant appears once per (gene, category) pair it belongs to, so the
509,570 rows cover 393,303 distinct variants over 21,458 genes. 1,651
rows carry `.` in place of an rsid; those variants are matched on
position and alleles instead.

The weight's scale and sign are category-specific — the JEPEG method
uses a different proxy measure for each, a deleteriousness score for
protein function and a free-energy change for microRNA binding, among
others — so weights are comparable within a category and not across
categories. `ReadAnnotation` recognises exactly the six category names
below and skips any row with another, reporting the count.

## References

Lee D, Williamson VS, Bigdeli TB, Riley BP, Webb BT, Fanous AH, Kendler
KS, Vladimirov VI, Bacanu SA (2015). JEPEG: a summary statistics based
tool for gene-level joint testing of functional variants.
*Bioinformatics* 31(8):1176-1182.
[doi:10.1093/bioinformatics/btu816](https://doi.org/10.1093/bioinformatics/btu816)

## Examples

``` r
if (FALSE) { # \dontrun{
data(JEPEG_SNP_Annotation.v1.0)
str(JEPEG_SNP_Annotation.v1.0)
} # }
```
