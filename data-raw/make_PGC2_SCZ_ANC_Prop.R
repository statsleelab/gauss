## Rebuild data/PGC2_SCZ_ANC_Prop.RData with the mixing proportions rescaled to
## sum to one.
##
##   Rscript data-raw/make_PGC2_SCZ_ANC_Prop.R
##
## Provenance. The published values are afmix() output on
## data/PGC2_Chr22_ilmn1M_AF1.txt against the 33KG panel, with afmix()'s first
## column (sup.pop) dropped -- the recipe is preserved as a commented-out
## `wgt.df[,-1]` in vignettes/afmix_example.Rmd. afmix() estimates the weights
## by regressing study allele frequency on the panel populations' allele
## frequencies, clips negative coefficients to zero, rounds to three decimals,
## and returns only the populations with a positive weight (21 of the panel's
## 29). It does not rescale, and clipping the negative coefficients upward
## pushes the total above one, so the published table summed to 1.061.
##
## The mixture covariance the analyses build (DISTMIX Eq. 2) assumes the
## proportions sum to one, so the values are rescaled here rather than left for
## every consumer to fix. Only the scale changes: the proportions relative to
## one another, and therefore every population's share, are untouched.

published <- data.frame(
  pop = c("ACB", "ASW", "BEB", "CCE", "CCS", "CDX", "CEU", "CLM", "CNE",
          "CSE", "FIN", "GBR", "GIH", "IBS", "JPT", "KHV", "MXL", "ORK",
          "PJL", "PUR", "TSI"),
  wgt = c(0.006, 0.036, 0.005, 0.008, 0.004, 0.018, 0.165, 0.025, 0.003,
          0.012, 0.138, 0.165, 0.006, 0.099, 0.011, 0.017, 0.030, 0.166,
          0.016, 0.045, 0.086),
  stringsAsFactors = FALSE
)

stopifnot(nrow(published) == 21L, all(published$wgt > 0))
cat(sprintf("published sum: %.6f\n", sum(published$wgt)))

PGC2_SCZ_ANC_Prop <- published
PGC2_SCZ_ANC_Prop$wgt <- published$wgt / sum(published$wgt)

## Relative proportions must be unchanged.
stopifnot(
  abs(sum(PGC2_SCZ_ANC_Prop$wgt) - 1) < 1e-12,
  all(abs(PGC2_SCZ_ANC_Prop$wgt / published$wgt -
            1 / sum(published$wgt)) < 1e-12),
  identical(PGC2_SCZ_ANC_Prop$pop, published$pop)
)

save(PGC2_SCZ_ANC_Prop, file = "data/PGC2_SCZ_ANC_Prop.RData")
cat(sprintf("rescaled sum:  %.15f\n", sum(PGC2_SCZ_ANC_Prop$wgt)))
print(PGC2_SCZ_ANC_Prop)
