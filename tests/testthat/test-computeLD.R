## computeLD() / simulateLD(): the mixture LD matrix over the measured SNPs.

fx <- panel()
wgt <- c(0.5, 0.3, 0.2)

ld_run <- function(pop_wgt_df = panel_wgt_df(wgt), start = 50000L,
                   end = 150000L) {
  quietly(gauss::computeLD(chr = 1L, start_bp = start, end_bp = end,
                           pop_wgt_df = pop_wgt_df, input_file = fx$z,
                           reference_index_file = fx$index,
                           reference_data_file = fx$data,
                           reference_pop_desc_file = fx$pop_desc))
}

test_that("computeLD() returns a SNP list and a matching correlation matrix", {
  res <- ld_run()

  expect_named(res, c("snplist", "cormat"))
  expect_s3_class(res$snplist, "data.frame")
  expect_named(res$snplist,
               c("rsid", "chr", "bp", "a1", "a2", "af1mix"))
  n <- nrow(res$snplist)
  expect_equal(dim(res$cormat), c(n, n))
  ## Only the measured SNPs inside the window, no wing.
  expect_equal(res$snplist$bp, seq(51000L, 149000L, by = 2000L))
})

test_that("the correlation matrix is symmetric with a unit diagonal", {
  res <- ld_run()

  expect_equal(diag(res$cormat), rep(1, nrow(res$cormat)))
  expect_equal(res$cormat, t(res$cormat))
  expect_false(anyNA(res$cormat))
  expect_true(all(abs(res$cormat) <= 1 + 1e-12))
})

test_that("the correlation matrix is exactly what the kernel produces", {
  ## Cross-check the whole pipeline -- index lookup, seek, genotype parse,
  ## matrix assembly -- against the kernel run directly on the fixture, with
  ## no reference-panel machinery in between.
  res <- ld_run()
  rows <- panel_row(res$snplist$bp)
  geno <- lapply(rows, panel_geno_vec)

  sd_vec <- vapply(geno, function(g) sqrt(.cal_wgt_cov(g, g, wgt)), numeric(1))
  expected <- outer(seq_along(rows), seq_along(rows), Vectorize(function(i, j) {
    if (i == j) return(1)
    .cal_wgt_cov(geno[[i]], geno[[j]], wgt) / (sd_vec[i] * sd_vec[j])
  }))

  expect_equal(res$cormat, expected)
})

test_that("the correlation matrix is Equation (2) of the DISTMIX paper", {
  ## The end-to-end anchor: no C++ in the reference at all.  ref_wgt_cov() is
  ## Eq. (2) transcribed from the paper and ref_wgt_cov_cor() applies the
  ## paper's step 2 (divide by the genotype standard deviations), which is
  ## exactly what computeLD() documents itself as producing.
  res <- ld_run()
  rows <- panel_row(res$snplist$bp)
  dosage <- lapply(rows, panel_dosage)

  expected <- outer(seq_along(rows), seq_along(rows), Vectorize(function(i, j) {
    if (i == j) return(1)
    ref_wgt_cov_cor(dosage[[i]], dosage[[j]], wgt)
  }))

  expect_equal(res$cormat, expected)
})

test_that("computeLD() reports the mixture allele frequency", {
  res <- ld_run()
  expected <- drop(panel_records()$af1[panel_row(res$snplist$bp), ] %*% wgt)
  expect_equal(res$snplist$af1mix, expected, tolerance = 1e-9)
})

test_that("computeLD() results are stable", {
  res <- ld_run()
  expect_snapshot_value(round(res$cormat[1:10, 1:10], 8), style = "deparse",
                        tolerance = 1e-6)
})

test_that("simulateLD() computes LD from a resampled cohort", {
  ## simulateLD() draws individuals from the panel with its own mt19937, so
  ## the numbers are not reproducible from an R seed -- only the shape of the
  ## result is checked here.
  res <- quietly(gauss::simulateLD(chr = 1L, start_bp = 50000L,
                                   end_bp = 150000L,
                                   pop_wgt_df = panel_wgt_df(wgt),
                                   sim_size = 100L, input_file = fx$z,
                                   reference_index_file = fx$index,
                                   reference_data_file = fx$data,
                                   reference_pop_desc_file = fx$pop_desc))
  ld <- ld_run()

  expect_named(res, c("snplist", "cormat"))
  expect_equal(res$snplist$bp, ld$snplist$bp)
  expect_equal(dim(res$cormat), dim(ld$cormat))
  expect_equal(diag(res$cormat), rep(1, nrow(res$cormat)))
  expect_equal(res$cormat, t(res$cormat))
})
