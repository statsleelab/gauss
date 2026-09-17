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

dist_unmatched <- function() {
  quietly(gauss::dist(chr = region$chr, start_bp = region$start,
                      end_bp = region$end, wing_size = region$wing,
                      study_pop = "PPA", input_file = fx$z_unmatched,
                      reference_index_file = fx$index,
                      reference_data_file = fx$data,
                      reference_pop_desc_file = fx$pop_desc))
}

test_that("study variants absent from the panel are reported as type 2", {
  ## Real GWAS input always contains variants the panel does not carry. They
  ## have no genotype record, so they carry no panel allele frequency and take
  ## no part in the imputation -- but they are the user's own measurements and
  ## are passed through to the output.
  res <- dist_unmatched()
  extra <- res[res$type == 2L, ]
  study <- utils::read.table(fx$z_unmatched, header = TRUE,
                             stringsAsFactors = FALSE)
  study <- study[grepl("^rsx", study$rsid), ]

  expect_equal(nrow(extra), nrow(study))
  expect_equal(extra$rsid, study$rsid)
  expect_equal(extra$bp, study$bp)
  ## No panel record means no panel frequency: NA, not a number read off some
  ## other variant's record.
  expect_true(all(is.na(extra$af1ref)))
  ## The study's own Z-score reaches the output untouched.
  expect_equal(extra$z, study$z)
  expect_equal(extra$info, rep(1, nrow(extra)))
})

test_that("unmatched variants do not perturb the imputation", {
  ## The decisive invariant: adding variants the panel does not carry must
  ## leave every panel variant's result exactly as it was.
  res <- dist_run()
  with_extra <- dist_unmatched()
  panel_rows <- with_extra[with_extra$type != 2L, ]

  expect_equal(panel_rows$rsid, res$rsid)
  expect_equal(panel_rows$z, res$z)
  expect_equal(panel_rows$info, res$info)
  expect_equal(panel_rows$af1ref, res$af1ref)
})

test_that("the output is deterministic when the panel is missing variants", {
  ## The seek to a non-existent record used to fail unchecked, leaving the
  ## stream wherever it happened to be, so which rows appeared and what
  ## frequency they carried varied between runs.
  expect_equal(dist_unmatched(), dist_unmatched())
})

test_that("an index entry with no genotype record is a clear error", {
  ## Guards the checked seek itself.  A negative file offset can only mean the
  ## index and the genotype file do not belong together; it used to be passed
  ## to bgzf_seek() unchecked, and the caller then parsed whatever record the
  ## stream was sitting on.
  expect_error(
    quietly(gauss::dist(chr = region$chr, start_bp = region$start,
                        end_bp = region$end, wing_size = region$wing,
                        study_pop = "PPA", input_file = fx$z,
                        reference_index_file = fx$index_bad,
                        reference_data_file = fx$data,
                        reference_pop_desc_file = fx$pop_desc)),
    "no record"
  )
})

test_that("distmix() reports unmatched variants with an NA mixture frequency", {
  res <- quietly(gauss::distmix(chr = region$chr, start_bp = region$start,
                                end_bp = region$end, wing_size = region$wing,
                                pop_wgt_df = panel_wgt_df(),
                                input_file = fx$z_unmatched,
                                reference_index_file = fx$index,
                                reference_data_file = fx$data,
                                reference_pop_desc_file = fx$pop_desc))
  extra <- res[res$type == 2L, ]
  ref <- distmix_run()

  expect_equal(nrow(extra), 5L)
  expect_true(all(is.na(extra$af1mix)))
  expect_equal(res$z[res$type != 2L], ref$z)
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
