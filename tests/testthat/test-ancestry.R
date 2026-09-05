## Ancestry-proportion estimation: afmix()/cpw2() from study allele
## frequencies, zmix() from study Z-scores.

fx <- panel()

## The fixture's allele-frequency file is an exact 50/30/20 mixture of the
## three panel populations, so afmix() has a known target to recover.
target <- c(0.5, 0.3, 0.2)

test_that("afmix() recovers the mixing proportions of the study cohort", {
  res <- quietly(gauss::afmix(input_file = fx$af1,
                              reference_index_file = fx$index,
                              reference_data_file = fx$data,
                              reference_pop_desc_file = fx$pop_desc,
                              interval = 1L))

  expect_named(res, c("sup.pop", "pop", "wgt"))
  expect_equal(res$pop, panel_pops)
  expect_equal(res$sup.pop, panel_sup_pops)
  expect_equal(res$wgt, target, tolerance = 1e-3)
  expect_equal(sum(res$wgt), 1, tolerance = 1e-6)
})

test_that("afmix() output needs its first column dropped to feed distmix()", {
  ## pop_wgt_df is read positionally as (population, weight), but afmix()
  ## returns (sup.pop, pop, wgt).  Handing it over unchanged makes distmix()
  ## read superpopulation names as populations and population names as
  ## weights, so the first column has to go first.
  ##
  ## The failing call is deliberately not exercised here: a non-numeric weight
  ## column reaches Rcpp's type cast, which aborts the process outright unless
  ## the package was compiled with NDEBUG.
  res <- quietly(gauss::afmix(input_file = fx$af1,
                              reference_index_file = fx$index,
                              reference_data_file = fx$data,
                              reference_pop_desc_file = fx$pop_desc,
                              interval = 1L))

  expect_equal(names(res)[1:2], c("sup.pop", "pop"))
  expect_type(res[[2]], "character")
  expect_type(res[[3]], "double")

  ## Dropping it gives the canonical two-column shape, which works.
  ok <- quietly(gauss::distmix(chr = 1L, start_bp = 50000L, end_bp = 150000L,
                               wing_size = 20000L, pop_wgt_df = res[, -1],
                               input_file = fx$z,
                               reference_index_file = fx$index,
                               reference_data_file = fx$data,
                               reference_pop_desc_file = fx$pop_desc))
  expect_s3_class(ok, "data.frame")
  expect_equal(nrow(ok), 101L)
})

test_that("cpw2() estimates population weights from allele frequencies", {
  res <- quietly(gauss::cpw2(input_file = fx$af1,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc,
                             interval = 1L))

  expect_named(res, c("pop", "wgt"))
  expect_equal(res$pop, panel_pops)
  expect_true(all(res$wgt >= 0))
  expect_equal(res$wgt, target, tolerance = 0.05)
  ## Unlike afmix(), cpw2() does not rescale its weights to sum to one.
  expect_equal(sum(res$wgt), 1, tolerance = 0.05)
})

test_that("zmix() returns a proper set of population weights", {
  res <- quietly(gauss::zmix(input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc,
                             percentile = 0.5, interval = 2,
                             level = "population"))

  expect_named(res, c("Population", "SuperPopulation", "Weight"))
  expect_equal(res$Population, panel_pops)
  expect_equal(res$SuperPopulation, panel_sup_pops)
  expect_true(all(res$Weight >= 0))
  expect_equal(sum(res$Weight), 1, tolerance = 1e-6)
})

test_that("zmix() can report weights per superpopulation", {
  res <- quietly(gauss::zmix(input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc,
                             percentile = 0.5, interval = 2,
                             level = "superpopulation"))

  expect_named(res, c("SuperPopulation", "Weight"))
  expect_equal(res$SuperPopulation, unique(panel_sup_pops))
  expect_true(all(res$Weight >= 0))
  expect_equal(sum(res$Weight), 1, tolerance = 1e-6)
})

test_that("zmix() rejects an unknown level", {
  expect_error(gauss::zmix(input_file = fx$z,
                           reference_index_file = fx$index,
                           reference_data_file = fx$data,
                           reference_pop_desc_file = fx$pop_desc,
                           level = "family"))
})

test_that("zmix() needs the documented columns in the description file", {
  desc <- readLines(fx$pop_desc)
  desc[1] <- sub("Population_Abbreviation", "pop", desc[1])
  expect_error(quietly(gauss::zmix(input_file = fx$z,
                                   reference_index_file = fx$index,
                                   reference_data_file = fx$data,
                                   reference_pop_desc_file = temp_lines(desc),
                                   percentile = 0.5, interval = 2)),
               "Population_Abbreviation")
})

test_that("prep_zmix5() returns the design matrix behind zmix()", {
  mat <- quietly(gauss::prep_zmix5(input_file = fx$z,
                                   reference_index_file = fx$index,
                                   reference_data_file = fx$data,
                                   reference_pop_desc_file = fx$pop_desc,
                                   percentile = 0.5, interval = 2))
  mat <- as.matrix(mat)

  ## Column 1 is z_i * z_j for a SNP pair; the rest are the per-population LD
  ## correlations for that pair -- not the population IDs and weights the
  ## documentation promises.
  expect_equal(ncol(mat), length(panel_pops) + 1L)
  expect_gt(nrow(mat), 0)
  expect_true(all(abs(mat[, -1]) <= 1 + 1e-12))
})

test_that("ancestry estimates are stable", {
  res <- quietly(gauss::zmix(input_file = fx$z,
                             reference_index_file = fx$index,
                             reference_data_file = fx$data,
                             reference_pop_desc_file = fx$pop_desc,
                             percentile = 0.5, interval = 2,
                             level = "population"))
  expect_snapshot_value(round(res$Weight, 6), style = "deparse",
                        tolerance = 1e-5)
})
