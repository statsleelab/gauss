## The correlation/covariance kernels in src/util.cpp are the arithmetic every
## analysis in the package reduces to.  They stream over ASCII genotype strings,
## so they are checked here against R's cor()/cov() on the same dosages.
##
## The wrappers (.cal_cor_pooled and friends, src/test_kernels.cpp) exist only
## for these tests and are not part of the public API.

variants <- c(11L, 12L, 40L, 41L, 97L, 150L)
wgt <- c(0.5, 0.3, 0.2)

test_that("CalCor over one genotype string is Pearson's r", {
  for (i in variants) {
    for (k in seq_along(panel_pops)) {
      g <- panel_geno_vec(i)[k]
      h <- panel_geno_vec(i + 1L)[k]
      d_g <- panel_dosage(i)[[k]]
      d_h <- panel_dosage(i + 1L)[[k]]
      expect_equal(.cal_cor_one(g, h), stats::cor(d_g, d_h),
                   info = paste("variant", i, "population", panel_pops[k]))
    }
  }
})

test_that("CalCor over a population vector pools the samples", {
  for (i in variants) {
    x <- panel_dosage(i)
    y <- panel_dosage(i + 1L)
    expect_equal(.cal_cor_pooled(panel_geno_vec(i), panel_geno_vec(i + 1L)),
                 ref_pooled_cor(x, y),
                 info = paste("variant", i))
  }
})

test_that("CalCor pooled over a subset uses only the populations given", {
  subset <- c(1L, 3L)
  for (i in variants) {
    expect_equal(
      .cal_cor_pooled(panel_geno_vec(i, subset), panel_geno_vec(i + 1L, subset)),
      ref_pooled_cor(panel_dosage(i, subset), panel_dosage(i + 1L, subset)),
      info = paste("variant", i)
    )
  }
})

test_that("CalCorSup pools the populations of one superpopulation", {
  ## pop_idx is zero-based on the C++ side.
  for (sup in unique(panel_sup_pops)) {
    idx <- which(panel_sup_pops == sup)
    for (i in variants) {
      expect_equal(
        .cal_cor_sup(panel_geno_vec(i), panel_geno_vec(i + 1L), idx - 1L),
        ref_pooled_cor(panel_dosage(i, idx), panel_dosage(i + 1L, idx)),
        info = paste(sup, "variant", i)
      )
    }
  }
})

test_that("CalCorSup on every population equals the pooled correlation", {
  for (i in variants) {
    expect_equal(
      .cal_cor_sup(panel_geno_vec(i), panel_geno_vec(i + 1L),
                   seq_along(panel_pops) - 1L),
      .cal_cor_pooled(panel_geno_vec(i), panel_geno_vec(i + 1L)),
      info = paste("variant", i)
    )
  }
})

test_that("the weighted CalCor is the mixture correlation", {
  for (i in variants) {
    expect_equal(
      .cal_cor_wgt(panel_geno_vec(i), panel_geno_vec(i + 1L), wgt),
      ref_wgt_cor(panel_dosage(i), panel_dosage(i + 1L), wgt),
      info = paste("variant", i)
    )
  }
})

test_that("equal weights on identically sized populations pool the samples", {
  even <- rep(1 / length(panel_pops), length(panel_pops))
  for (i in variants) {
    ## Every population here has 20 samples, so weighting them equally is the
    ## same as pooling -- up to the m vs m-1 difference, which cancels in a
    ## correlation.
    expect_equal(
      .cal_cor_wgt(panel_geno_vec(i), panel_geno_vec(i + 1L), even),
      .cal_cor_pooled(panel_geno_vec(i), panel_geno_vec(i + 1L)),
      info = paste("variant", i)
    )
  }
})

test_that("every kernel returns 1 for a variant against itself", {
  for (i in variants) {
    g <- panel_geno_vec(i)
    expect_equal(.cal_cor_pooled(g, g), 1)
    expect_equal(.cal_cor_wgt(g, g, wgt), 1)
    expect_equal(.cal_cor_sup(g, g, seq_along(panel_pops) - 1L), 1)
  }
})

test_that("CalWgtCov is Equation (2) of the DISTMIX paper", {
  ## The mixture covariance every mixture analysis reduces to -- distmix,
  ## computeLD, qcatmix, jepegmix, prep_recessive_impute.  ref_wgt_cov() is
  ## Eq. (2) transcribed from the paper; see helper-panel.R.
  ##
  ## This assertion is the whole reason the fixture exists: the kernel used to
  ## scale the within-population term by m^2, which overrode the caller's
  ## mixing proportions with panel sample sizes and numerically erased the
  ## between-population term.
  for (i in variants) {
    x <- panel_dosage(i)
    y <- panel_dosage(i + 1L)
    expect_equal(.cal_wgt_cov(panel_geno_vec(i), panel_geno_vec(i + 1L), wgt),
                 ref_wgt_cov(x, y, wgt), info = paste("variant", i))
  }
})

test_that("CalWgtCov is more than the within-population average", {
  ## The between-population term is the reason mixture LD differs from
  ## within-population LD.  Compare the kernel against a within-only weighted
  ## average and require that the difference actually matters: under the m^2
  ## defect the between term survived arithmetically but was suppressed by
  ## 1/m^2 relative to the within term, so this ratio collapsed to ~1e-4.
  rel <- vapply(seq_len(150L), function(i) {
    x <- panel_dosage(i)
    y <- panel_dosage(i + 1L)
    got <- .cal_wgt_cov(panel_geno_vec(i), panel_geno_vec(i + 1L), wgt)
    within_only <- sum(wgt * mapply(stats::cov, x, y))
    abs((got - within_only) / got)
  }, numeric(1))

  ## Median over 150 variant pairs, not a per-pair bound: individual pairs can
  ## legitimately have near-identical frequencies across populations.
  expect_gt(stats::median(rel), 0.02)
})

test_that("a single population reduces to that population's covariance", {
  ## With one population the between term vanishes and Eq. (2) collapses to
  ## C_1, so the kernel must agree with stats::cov() exactly.
  for (i in variants) {
    for (k in seq_along(panel_pops)) {
      g <- panel_geno_vec(i)[k]
      h <- panel_geno_vec(i + 1L)[k]
      expect_equal(.cal_wgt_cov(g, h, 1),
                   stats::cov(panel_dosage(i)[[k]], panel_dosage(i + 1L)[[k]]),
                   info = paste("variant", i, "population", panel_pops[k]))
    }
  }
})

test_that("CalWgtCov is symmetric and its diagonal is positive", {
  for (i in variants) {
    g <- panel_geno_vec(i)
    h <- panel_geno_vec(i + 1L)
    expect_equal(.cal_wgt_cov(g, h, wgt), .cal_wgt_cov(h, g, wgt))
    expect_gt(.cal_wgt_cov(g, g, wgt), 0)
  }
})

test_that("a monomorphic variant makes the kernels return NaN", {
  ## Documents how the kernels behave on degenerate input: the correlation is
  ## an unguarded 0/0, and the NaN then propagates through the LD matrix.  The
  ## fixture deliberately contains no such variant.
  const <- strrep("1", panel_pop_size[1])
  varying <- panel_geno_vec(11L)[1]

  expect_true(is.nan(.cal_cor_one(const, varying)))
  expect_true(is.nan(.cal_cor_pooled(const, varying)))
  expect_true(is.nan(.cal_cor_wgt(const, varying, 1)))
  expect_true(is.nan(.cal_cor_sup(const, varying, 0L)))
})
