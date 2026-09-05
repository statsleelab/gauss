## JEPEG / JEPEGMIX: gene-level tests over annotated functional categories.

fx <- panel()

jepeg_run <- function(annotation_file = fx$annot) {
  quietly(gauss::jepeg(study_pop = "PPA", input_file = fx$z,
                       annotation_file = annotation_file,
                       reference_index_file = fx$index,
                       reference_data_file = fx$data,
                       reference_pop_desc_file = fx$pop_desc))
}

jepegmix_run <- function(pop_wgt_df = panel_wgt_df()) {
  quietly(gauss::jepegmix(pop_wgt_df = pop_wgt_df, input_file = fx$z,
                          annotation_file = fx$annot,
                          reference_index_file = fx$index,
                          reference_data_file = fx$data,
                          reference_pop_desc_file = fx$pop_desc))
}

test_that("jepeg() reports one row per annotated gene", {
  res <- jepeg_run()

  expect_named(res, c("geneid", "chisq", "df", "jepeg_pval", "num_snp",
                      "top_categ", "top_categ_pval", "top_snp", "top_snp_pval"))
  expect_equal(res$geneid, sprintf("GENE%02d", 1:8))
  expect_true(all(res$chisq >= 0))
  expect_true(all(res$df >= 1))
  expect_true(all(res$jepeg_pval > 0 & res$jepeg_pval <= 1))
  expect_true(all(res$num_snp > 0))
  expect_true(all(res$top_snp %in% panel_records()$index$rsid))
  ## The output labels, not the annotation file's category names.
  expect_true(all(res$top_categ %in%
                    c("PFS", "TFB", "STR", "TAR", "CIS", "TRN")))
})

test_that("jepegmix() reports the same genes over a population mixture", {
  res <- jepegmix_run()
  ref <- jepeg_run()

  expect_named(res, names(ref))
  expect_equal(res$geneid, ref$geneid)
  expect_true(all(res$jepeg_pval > 0 & res$jepeg_pval <= 1))
  expect_false(isTRUE(all.equal(res$chisq, ref$chisq)))
})

test_that("jepegmix() weights are matched by population name, not by order", {
  ordered <- jepegmix_run()
  shuffled <- jepegmix_run(panel_wgt_df(wgt = c(0.2, 0.5, 0.3),
                                        pops = c("PPC", "PPA", "PPB")))
  expect_equal(shuffled$chisq, ordered$chisq)
})

test_that("an unknown annotation category is reported and the row dropped", {
  ## An unrecognised category used to be written into an uninitialised slot of
  ## the weight matrix; it must now warn and skip the row instead.
  annot <- readLines(fx$annot)
  annot[2:4] <- sub("(\\S+)(\\s+\\S+)$", "NOT_A_CATEGORY\\2", annot[2:4])
  tmp <- temp_lines(annot)

  expect_warning(res <- jepeg_run(tmp))
  expect_s3_class(res, "data.frame")
  expect_true(all(res$jepeg_pval > 0 & res$jepeg_pval <= 1))
})

test_that("jepeg() results are stable", {
  res <- jepeg_run()
  expect_snapshot_value(round(res$chisq, 8), style = "deparse",
                        tolerance = 1e-6)
  expect_snapshot_value(round(res$jepeg_pval, 8), style = "deparse",
                        tolerance = 1e-6)
})

test_that("jepegmix() results are stable", {
  res <- jepegmix_run()
  expect_snapshot_value(round(res$chisq, 8), style = "deparse",
                        tolerance = 1e-6)
})
