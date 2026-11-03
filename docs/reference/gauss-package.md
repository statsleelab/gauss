# gauss: Genome Analysis Using Summary Statistics

Analyses summary statistics from genome-wide association studies of
cohorts with mixed ancestry, against an external reference panel of
32,953 samples drawn from 29 populations. The panel is read directly off
disk by position, so it is never loaded into memory and is not
distributed with the package.

## Details

The entry points fall into four groups:

- Imputation of unmeasured variants:

  [`dist`](https://github.com/statsleelab/gauss/reference/dist.md) for a
  single population or superpopulation,
  [`distmix`](https://github.com/statsleelab/gauss/reference/distmix.md)
  for an arbitrary ancestry mixture.

- Gene-level association:

  [`jepeg`](https://github.com/statsleelab/gauss/reference/jepeg.md) and
  [`jepegmix`](https://github.com/statsleelab/gauss/reference/jepegmix.md),
  which aggregate functional annotations into a single statistic per
  gene.

- Region-level association:

  [`qcat`](https://github.com/statsleelab/gauss/reference/qcat.md) and
  [`qcatmix`](https://github.com/statsleelab/gauss/reference/qcatmix.md).

- Ancestry and linkage disequilibrium:

  [`afmix`](https://github.com/statsleelab/gauss/reference/afmix.md) and
  [`zmix`](https://github.com/statsleelab/gauss/reference/zmix.md)
  estimate mixing proportions from the study's own allele frequencies or
  Z-scores;
  [`computeLD`](https://github.com/statsleelab/gauss/reference/computeLD.md)
  returns the ancestry-informed LD matrix the other functions use
  internally.

Functions whose name ends in `mix` take a `pop_wgt_df` of mixing
proportions and are read positionally: the first column must name the
population and the second must carry its weight. See
`data(PGC2_SCZ_ANC_Prop)` for the expected shape.

## References

Lee D, Williamson VS, Bigdeli TB, Riley BP, Fanous AH, Vladimirov VI,
Bacanu SA (2015). DIST: direct imputation of summary statistics for
unmeasured SNPs. *Bioinformatics* 31(19):3099-3104.
[doi:10.1093/bioinformatics/btv348](https://doi.org/10.1093/bioinformatics/btv348)

## See also

Useful links:

- <https://github.com/statsleelab/gauss>

- Report bugs at <https://github.com/statsleelab/gauss/issues>

## Author

**Maintainer**: Donghyung Lee <leed13@miamioh.edu>

Authors:

- Donghyung Lee <leed13@miamioh.edu>
