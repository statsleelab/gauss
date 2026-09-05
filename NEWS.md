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
