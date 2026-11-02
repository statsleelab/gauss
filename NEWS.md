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

* **Variants that the reference panel does not carry are handled correctly.**
  `MakeSnpVec()` and `MakeSnpVecMix()` looped over every SNP in the map,
  including those present in the user's input but absent from the panel. Such a
  SNP has no genotype record, so its stored file offset is still the default
  `-1`; `bgzf_seek()` failed, no call site checked the return value, and the
  next read returned whichever record the stream happened to be sitting on. The
  variant therefore acquired an unrelated variant's allele frequency, and that
  number then decided — via the minor-allele-frequency cutoff, on evidence
  belonging to a different variant — whether it appeared in the output at all.

  Such variants are now passed through to the output with `type = 2` and an
  `NA` allele frequency, which is what every analysis already assumed: each one
  carries an explicit "if (type == 2) do nothing" branch excluding them from the
  linear algebra. Their own Z-score reaches the output untouched, and the count
  is reported on the console.

  This adds rows to `dist()`, `distmix()`, `qcat()`, `qcatmix()`, `prep_qcat()`
  and `prep_recessive_impute()` output — one per input variant in the reporting
  window that the panel lacks — and those rows carry `NA` in the `af1ref` /
  `af1mix` column. Results for panel variants are unchanged: on a 1 Mb window of
  chromosome 22 against the 33KG panel, one variant of 5869 is recovered this
  way and the other 5868 rows are identical.

* **Every seek into the reference panel is checked.** A file offset that is
  negative, or a seek that fails, now raises an error naming the offset instead
  of silently leaving the stream somewhere else and parsing the wrong record.
  This covers all seven live call sites.

* **`pop_wgt_df` is now validated, and every function that takes ancestry
  proportions rescales them to sum to one before doing any work.** Equation (2)
  assumes `sum_k w_k = 1`; its between-population term is a centred outer
  product, and above 1 it can go negative far enough to make a variant's
  mixture variance negative, at which point the square root in the LD
  normalisation produces `NaN` that spreads through the whole LD matrix. The
  `m_k^2` inflation described above had been masking this, so the two fixes
  belong together.

  This was not hypothetical: the mixing proportions shipped as
  `data(PGC2_SCZ_ANC_Prop)`, used throughout the README and the vignettes,
  summed to 1.061. On the real 33KG panel that produced 402 `NaN` entries
  across every row of a 202-variant LD matrix once the covariance was
  corrected.

  Rescaling is reported on the console rather than raised as a warning:
  `afmix()` and `cpw2()` return unrescaled weights by construction, so this is
  ordinary preprocessing, not an exceptional condition. Only the scale changes
  — each population's share relative to the others is preserved, and passing
  the same proportions at any scale gives identical results.

  Populations named in `pop_wgt_df` that the reference panel does not carry are
  reported by name instead of having their weight silently dropped; a mixture
  in which no population matches the panel, a negative or non-finite weight,
  and weights summing to zero are all now clear errors.

## Data

* `data(PGC2_SCZ_ANC_Prop)` is now rescaled to sum to 1 (it summed to 1.061).
  Only the scale changed: every population's share relative to the others, and
  the set of 21 populations, are exactly as published. The derivation is now
  recorded in `data-raw/make_PGC2_SCZ_ANC_Prop.R`, and the help page documents
  the provenance and the column meanings.

## Build and packaging

* `src/Makevars` is now shipped directly rather than generated by a `configure`
  script whose only job was to copy `src/Makevars.in` into place, and it sets
  `PKG_CXXFLAGS` with `=` instead of `+=`. `+=` is a GNU make extension that the
  package never declared in `SystemRequirements`, so a build under a non-GNU
  make could silently drop the flag.

* `stats` and `utils` are declared in `Imports`, and the symbols `fiqt()` uses
  (`pnorm`, `qnorm`, `p.adjust`) are imported explicitly. They previously
  resolved only because those packages happen to be attached in a default R
  session. `RcppEigen` has moved out of `Imports`: it is a headers-only
  dependency and `LinkingTo` already covers it.

* `fiqt()` spells out `lower.tail = FALSE` where it previously relied on
  partial matching of `low =`. The computation is unchanged.

* `src/zmixrcpp.cpp_` and `src/etc/`, both superseded before the Rcpp rewrite,
  are kept out of the source tarball.

## Documentation

* The `pop_wgt_df` argument is documented properly for `distmix()`,
  `computeLD()`, `qcatmix()`, `jepegmix()`, `simulateLD()` and
  `prep_recessive_impute()`: that it is read positionally, that the weights are
  rescaled, and that `afmix()` returns three columns whose first must be
  dropped before use.

* `?gauss` is a real package help page — what the package does, the four groups
  of entry points, the positional `pop_wgt_df` convention and the DIST
  reference — in place of the `package.skeleton()` placeholder text that had
  been shipping since the first release.

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
