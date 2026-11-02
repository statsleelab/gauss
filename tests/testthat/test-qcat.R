## QCAT / QCATMIX: per-variant conditional association tests.

fx <- panel()
region <- list(chr = 1L, start = 50000L, end = 150000L, wing = 20000L)

qcat_run <- function() {
  quietly(gauss::qcat(chr = region$chr, start_bp = region$start,
                      end_bp = region$end, wing_size = region$wing,
                      study_pop = "PPA", input_file = fx$z,
                      reference_index_file = fx$index,
                      reference_data_file = fx$data,
                      reference_pop_desc_file = fx$pop_desc))
}

qcatmix_run <- function(pop_wgt_df = panel_wgt_df()) {
  quietly(gauss::qcatmix(chr = region$chr, start_bp = region$start,
                         end_bp = region$end, wing_size = region$wing,
                         pop_wgt_df = pop_wgt_df, input_file = fx$z,
                         reference_index_file = fx$index,
                         reference_data_file = fx$data,
                         reference_pop_desc_file = fx$pop_desc))
}

test_that("qcat() reports a chi-square test per variant in the window", {
  res <- qcat_run()

  expect_named(res, c("rsid", "chr", "bp", "a1", "a2", "af1ref", "z",
                      "qcat_m", "qcat_t", "qcat_chisq", "qcat_pval", "type"))
  in_window <- panel_bp >= region$start & panel_bp <= region$end
  expect_equal(res$bp, panel_bp[in_window])
  expect_true(all(res$qcat_chisq >= 0))
  expect_true(all(res$qcat_pval > 0 & res$qcat_pval <= 1))
  expect_false(anyNA(res$qcat_pval))
})

test_that("qcatmix() reports the same variants over a population mixture", {
  res <- qcatmix_run()
  ref <- qcat_run()

  ## Same columns, except that the mixture reports af1mix instead of af1ref.
  expect_named(res, c("rsid", "chr", "bp", "a1", "a2", "af1mix", "z",
                      "qcat_m", "qcat_t", "qcat_chisq", "qcat_pval", "type"))
  expect_equal(res$bp, ref$bp)
  expect_true(all(res$qcat_chisq >= 0))
  expect_true(all(res$qcat_pval > 0 & res$qcat_pval <= 1))
  ## A different LD backbone has to give different numbers.
  expect_false(isTRUE(all.equal(res$qcat_chisq, ref$qcat_chisq)))
})

test_that("qcatmix() weights are matched by population name, not by order", {
  ordered <- qcatmix_run()
  shuffled <- qcatmix_run(panel_wgt_df(wgt = c(0.2, 0.5, 0.3),
                                       pops = c("PPC", "PPA", "PPB")))
  expect_equal(shuffled$qcat_chisq, ordered$qcat_chisq)
})

test_that("qcat() results are stable", {
  res <- qcat_run()
  expect_snapshot_value(round(res$qcat_chisq, 8), style = "deparse",
                        tolerance = 1e-6)
})

test_that("qcatmix() results are stable", {
  res <- qcatmix_run()
  expect_snapshot_value(round(res$qcat_chisq, 8), style = "deparse",
                        tolerance = 1e-6)
})

test_that("prep_qcat() returns the pieces qcat() is built from", {
  res <- quietly(gauss::prep_qcat(chr = region$chr, start_bp = region$start,
                                  end_bp = region$end, wing_size = region$wing,
                                  study_pop = "PPA", input_file = fx$z,
                                  reference_index_file = fx$index,
                                  reference_data_file = fx$data,
                                  reference_pop_desc_file = fx$pop_desc))
  expect_type(res, "list")
  expect_gt(length(res), 0)
})
