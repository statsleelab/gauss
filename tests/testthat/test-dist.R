## DIST / DISTMIX: Gaussian imputation of Z-scores at panel variants that the
## study did not measure.

fx <- panel()
region <- list(chr = 1L, start = 50000L, end = 150000L, wing = 20000L)

dist_run <- function(...) {
  quietly(gauss::dist(chr = region$chr, start_bp = region$start,
                      end_bp = region$end, wing_size = region$wing,
                      study_pop = "PPA", input_file = fx$z,
                      reference_index_file = fx$index,
                      reference_data_file = fx$data,
                      reference_pop_desc_file = fx$pop_desc, ...))
}

distmix_run <- function(pop_wgt_df = panel_wgt_df(), ...) {
  quietly(gauss::distmix(chr = region$chr, start_bp = region$start,
                         end_bp = region$end, wing_size = region$wing,
                         pop_wgt_df = pop_wgt_df, input_file = fx$z,
                         reference_index_file = fx$index,
                         reference_data_file = fx$data,
                         reference_pop_desc_file = fx$pop_desc, ...))
}

test_that("dist() returns one row per panel variant in the prediction window", {
  res <- dist_run()

  expect_s3_class(res, "data.frame")
  expect_named(res, c("rsid", "chr", "bp", "a1", "a2", "af1ref", "z", "pval",
                      "info", "type"))
  in_window <- panel_bp >= region$start & panel_bp <= region$end
  expect_equal(res$bp, panel_bp[in_window])
  expect_equal(res$rsid, panel_records()$index$rsid[in_window])
  ## Half the panel is measured; the rest is imputed.
  expect_equal(sort(unique(res$type)), c(0L, 1L))
  expect_equal(sum(res$type == 1L), sum(in_window & panel_bp %% 2000 == 1000))
})

test_that("dist() reports the alleles and frequencies of the panel", {
  res <- dist_run()
  index <- panel_records()$index[panel_row(res$bp), ]

  expect_equal(res$a1, index$a1)
  expect_equal(res$a2, index$a2)
  ## af1ref is recomputed from the study population's genotypes, not taken
  ## from the index, and is rounded up to five decimals.
  expected <- vapply(panel_row(res$bp),
                     function(i) mean(panel_dosage(i)[[1]]) / 2, numeric(1))
  expect_true(all(res$af1ref >= expected - 1e-12))
  expect_true(all(res$af1ref - expected <= 1e-5 + 1e-12))
})

test_that("dist() leaves measured Z-scores untouched, in panel orientation", {
  res <- dist_run()
  study <- utils::read.table(fx$z, header = TRUE, stringsAsFactors = FALSE)
  measured <- res[res$type == 1L, ]
  study <- study[match(measured$bp, study$bp), ]

  ## Rows written with the alleles the other way round get their Z-score
  ## negated so that it refers to the panel's allele 1.
  flipped <- study$a1 != measured$a1
  expect_gt(sum(flipped), 0)
  expect_equal(measured$z, ifelse(flipped, -study$z, study$z))
  expect_equal(measured$info, rep(1, nrow(measured)))
})

test_that("dist() imputes finite Z-scores with usable information", {
  res <- dist_run()
  imputed <- res[res$type == 0L, ]

  expect_true(all(is.finite(imputed$z)))
  expect_true(all(imputed$info > 0 & imputed$info <= 1))
  ## Imputation quality should be high in a dense, high-LD region.
  expect_gt(stats::median(imputed$info), 0.5)
})

test_that("the reported p-value is the two-sided normal tail of z", {
  res <- dist_run()
  expect_equal(res$pval, 2 * stats::pnorm(abs(res$z), lower.tail = FALSE))
})

test_that("dist() accepts a superpopulation as the study population", {
  res <- dist_run()
  sup <- quietly(gauss::dist(chr = region$chr, start_bp = region$start,
                             end_bp = region$end, wing_size = region$wing,
                             study_pop = "SPX", input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc))

  expect_equal(sup$bp, res$bp)
  ## SPX pools PPA and PPB, so the frequencies differ from PPA alone.
  expected <- vapply(panel_row(sup$bp),
                     function(i) mean(unlist(panel_dosage(i, 1:2))) / 2,
                     numeric(1))
  expect_true(all(sup$af1ref >= expected - 1e-12))
  expect_true(all(sup$af1ref - expected <= 1e-5 + 1e-12))
  expect_false(isTRUE(all.equal(sup$z, res$z)))
})

test_that("dist() tolerates study variants that are absent from the panel", {
  ## Real GWAS input always contains variants the panel does not carry.
  res <- dist_run()
  with_extra <- quietly(gauss::dist(chr = region$chr, start_bp = region$start,
                                    end_bp = region$end, wing_size = region$wing,
                                    study_pop = "PPA",
                                    input_file = fx$z_unmatched,
                                    reference_index_file = fx$index,
                                    reference_data_file = fx$data,
                                    reference_pop_desc_file = fx$pop_desc))

  expect_equal(with_extra$bp, res$bp)
  expect_equal(with_extra$z, res$z)
})

test_that("dist() results are stable", {
  res <- dist_run()
  expect_snapshot_value(round(res$z, 8), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(round(res$info, 8), style = "deparse", tolerance = 1e-6)
})

test_that("distmix() imputes over a mixture of panel populations", {
  res <- distmix_run()

  expect_s3_class(res, "data.frame")
  expect_named(res, c("rsid", "chr", "bp", "a1", "a2", "af1mix", "z", "pval",
                      "info", "type"))
  in_window <- panel_bp >= region$start & panel_bp <= region$end
  expect_equal(res$bp, panel_bp[in_window])

  ## af1mix is the weighted average of the per-population frequencies stored
  ## in the genotype records.
  expected <- drop(panel_records()$af1[panel_row(res$bp), ] %*% c(0.5, 0.3, 0.2))
  expect_equal(res$af1mix, expected, tolerance = 1e-9)
  expect_true(all(is.finite(res$z)))
  expect_true(all(res$info[res$type == 0L] > 0))
})

test_that("distmix() weights are matched by population name, not by order", {
  ordered <- distmix_run(panel_wgt_df())
  shuffled <- distmix_run(panel_wgt_df(wgt = c(0.2, 0.5, 0.3),
                                       pops = c("PPC", "PPA", "PPB")))

  expect_equal(shuffled$af1mix, ordered$af1mix)
  expect_equal(shuffled$z, ordered$z)
})

test_that("distmix() over a single population still runs", {
  res <- distmix_run(panel_wgt_df(wgt = 1, pops = "PPA"))

  expect_equal(res$af1mix, panel_records()$af1[panel_row(res$bp), 1],
               tolerance = 1e-9)
  expect_true(all(is.finite(res$z)))
})

test_that("distmix() results are stable", {
  res <- distmix_run()
  expect_snapshot_value(round(res$z, 8), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(round(res$info, 8), style = "deparse", tolerance = 1e-6)
})
