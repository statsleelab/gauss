# gauss 1.1.0

## Bug fixes

* **`CalWgtCov()` now computes the mixture covariance the DISTMIX paper
  specifies, which changes the results of every mixture analysis.** Affected:
  `distmix()`, `computeLD()`, `qcatmix()`, `jepegmix()` and
  `prep_recessive_impute()`. Single-population and superpopulation analyses
  (`dist()`, `qcat()`, `jepeg()`), ancestry estimation (`afmix()`, `cpw2()`,
  `zmix()`) and `fiqt()` are unaffected.

  The ancestry-informed covariance is Equation (2) of Lee et al.,
  *Bioinformatics* 31(19):3099-3104 (doi:10.1093/bioinformatics/btv348),

      C = sum_k w_k C_k + sum_k w_k (mu_k - mu)(mu_k - mu)^T

  where `C_k` is the variance-covariance matrix of reference population *k* and
  `w_k` its mixing proportion. The within-population term `sum_k w_k C_k` was
  being scaled by `m_k^2`, the square of that population's panel sample size,
  because `(m*sumxy - sumx*sumy)` was never divided by `m*(m-1)`.

  Two consequences, both now gone. Populations were effectively re-weighted by
  `w_k * m_k^2` instead of `w_k`; since panel sizes span 52 to 6360 samples,
  user-supplied mixing proportions were largely overridden by panel sample
  sizes. And the between-population term, which is the reason a mixture LD
  matrix differs from a within-population one, was `m_k^2` times smaller than
  the inflated first term and so numerically discarded -- mixture LD behaved
  like a size-weighted within-population average.

  On the package's test fixture (3 populations of 20 samples, weights
  0.5/0.3/0.2), the previous release's numbers differ from the corrected ones
  by up to 0.41 in `distmix()` Z-scores (3 of 101 change sign), 0.20 in
  `computeLD()` correlations (14 of 100 change sign), 0.48 in `qcatmix()`
  chi-squares and 0.37 in `jepegmix()` gene chi-squares. The distortion grows
  with panel sample size, so differences on the 33KG panel are larger.

  Users who need to reproduce results from gauss 1.0.1 or earlier should pin
  that version.

## Known limitations

* The mixing proportions passed as `pop_wgt_df` are not checked to sum to one,
  which Equation (2) requires. Note that `cpw2()` does not rescale its own
  output, so its weights should be normalised before being handed to
  `distmix()` or `jepegmix()`.

# gauss 1.0.1

## Bug fixes

* `jepeg()` and `jepegmix()` no longer terminate the R session when the SNP
  annotation file cannot be opened. A mistyped `annotation_file` path now
  raises an ordinary R error, leaving the session and any unsaved work intact.

* Annotation rows whose `categ` value is not one of `PROTEIN`, `TFBS`,
  `WTH_HAIR`, `WTH_TARGET`, `CIS_EQTL` or `TRANS_EQTL` are now skipped, and
  the number skipped is reported as a warning. Previously such a row was
  assigned to whichever category the preceding row used, which silently
  corrupted the JEPEG weight matrix and the resulting gene-level p-values.
  Results are unchanged for annotation files that use only the six documented
  categories.
