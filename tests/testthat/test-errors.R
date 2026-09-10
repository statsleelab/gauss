## Error paths.  Every one of these must surface as a catchable R condition --
## the C++ layer used to call exit(), which killed the whole session.

fx <- panel()

dist_with <- function(...) {
  args <- list(chr = 1L, start_bp = 50000L, end_bp = 150000L,
               wing_size = 20000L, study_pop = "PPA", input_file = fx$z,
               reference_index_file = fx$index, reference_data_file = fx$data,
               reference_pop_desc_file = fx$pop_desc)
  args[names(list(...))] <- list(...)
  quietly(do.call(gauss::dist, args))
}

test_that("a missing file is named in the error message", {
  missing <- file.path(tempdir(), "does-not-exist")

  expect_error(dist_with(input_file = missing), "can't open input file")
  expect_error(dist_with(reference_index_file = missing),
               "can't open reference index file")
  expect_error(dist_with(reference_data_file = missing),
               "can't open reference data file")
  expect_error(dist_with(reference_pop_desc_file = missing),
               "can't open reference population description file")
})

test_that("an unknown study population is rejected", {
  expect_error(dist_with(study_pop = "NOPE"), "invalid population name")
})

test_that("too few variants in the window is an error, not a crash", {
  expect_error(dist_with(start_bp = 50000L, end_bp = 52000L, wing_size = 0L),
               "Not enough number of SNPs loaded")
})

test_that("a mixture whose populations are unknown is an error", {
  ## The unknown name is reported first, then the run stops because nothing is
  ## left to weight.
  expect_error(
    suppressWarnings(
      quietly(gauss::distmix(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                             wing_size = 20000L,
                             pop_wgt_df = panel_wgt_df(wgt = 1, pops = "NOPE"),
                             input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc))),
    "none of the populations"
  )
})

test_that("an empty input file is an error, not a crash", {
  empty <- temp_lines("rsid chr bp a1 a2 z")
  expect_error(dist_with(input_file = empty), "Not enough number of SNPs")
})

test_that("a missing annotation file is named in the error message", {
  expect_error(
    quietly(gauss::jepeg(study_pop = "PPA", input_file = fx$z,
                         annotation_file = file.path(tempdir(), "nope.txt"),
                         reference_index_file = fx$index,
                         reference_data_file = fx$data,
                         reference_pop_desc_file = fx$pop_desc)),
    "annotation"
  )
})

test_that("a duplicated study variant is reported", {
  z <- readLines(fx$z)
  dup <- temp_lines(c(z, z[2]))
  ## The duplicate collapses in the SNP map, so the run must either report it
  ## or ignore it -- never produce a different answer for the shared variants.
  res <- tryCatch(dist_with(input_file = dup), error = function(e) e)
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "duplicate")
  } else {
    expect_equal(res$z, dist_with()$z)
  }
})

test_that("mixing proportions are rescaled to sum to one", {
  ## DISTMIX Eq. (2) assumes sum_k w_k = 1.  Above 1 the between-population term
  ## can go negative far enough to make a variant's mixture variance negative,
  ## and the sqrt() in the LD normalisation then returns NaN.  afmix() and
  ## cpw2() return unrescaled weights by construction, so rescaling here is a
  ## normal preprocessing step, not an exception.
  ld <- function(w) {
    quietly(gauss::computeLD(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                             pop_wgt_df = w, input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc))
  }
  ref <- ld(panel_wgt_df())

  ## Same proportions at any scale must give exactly the same answer, with no
  ## NaN from a negative mixture variance.
  for (scale in c(1.061, 100, 0.25)) {
    res <- ld(panel_wgt_df(wgt = c(0.5, 0.3, 0.2) * scale))
    expect_false(anyNA(res$cormat), info = paste("scale", scale))
    expect_equal(res$cormat, ref$cormat, info = paste("scale", scale))
  }
})

test_that("rescaling is reported but is not a warning", {
  expect_no_warning(
    quietly(gauss::computeLD(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                             pop_wgt_df = panel_wgt_df(wgt = c(0.5, 0.3, 0.2) * 1.061),
                             input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc))
  )
})

test_that("populations absent from the panel are reported, not silently dropped", {
  extra <- rbind(panel_wgt_df(), data.frame(pop = "NOPE", wgt = 0))
  expect_warning(
    quietly(gauss::computeLD(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                             pop_wgt_df = extra, input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc)),
    "NOPE"
  )
})

test_that("a mixture with no usable population is a clear error", {
  expect_error(
    suppressWarnings(
      quietly(gauss::computeLD(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                               pop_wgt_df = panel_wgt_df(wgt = 1, pops = "NOPE"),
                               input_file = fx$z,
                               reference_index_file = fx$index,
                               reference_data_file = fx$data,
                               reference_pop_desc_file = fx$pop_desc))),
    "none of the populations"
  )
})

test_that("a negative mixing proportion is rejected", {
  expect_error(
    quietly(gauss::computeLD(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                             pop_wgt_df = panel_wgt_df(wgt = c(1.2, -0.2, 0.0)),
                             input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc)),
    "non-negative"
  )
})
