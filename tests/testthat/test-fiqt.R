## fiqt(): FDR inverse-quantile transformation (Winner's Curse correction).
## Pure R, no reference panel involved.

test_that("a single Z-score is returned unchanged", {
  ## With one test there is nothing to adjust for, so BH leaves the p-value
  ## alone and the transformation is the identity.
  expect_equal(fiqt(1.5), 1.5)
  expect_equal(fiqt(-2.25), -2.25)
  expect_equal(fiqt(0), 0)
})

test_that("two Z-scores match the value computed by hand", {
  ## z = c(1, 3): p = c(0.3173105, 0.0026998); BH gives
  ## c(0.3173105, 0.0053996); back-transforming yields:
  expect_equal(fiqt(c(1, 3)), c(1, 2.7821750), tolerance = 1e-6)
})

test_that("Z-scores are shrunk towards zero and keep their sign", {
  z <- c(-6, -3.5, -1.2, 0.4, 2, 4.8, 7.1)
  adj <- fiqt(z)

  expect_length(adj, length(z))
  expect_true(all(abs(adj) <= abs(z) + 1e-12))
  expect_equal(sign(adj), sign(z))
})

test_that("the transformation is odd in z", {
  z <- c(1, 2, 3, 4)
  expect_equal(fiqt(-z), -fiqt(z))
})

test_that("stronger signals stay stronger", {
  z <- sort(c(0.3, 1.1, 2.4, 3.9, 5.6))
  expect_true(all(diff(fiqt(z)) >= 0))
})

test_that("Z-scores beyond the min.p bound are passed through", {
  ## qnorm(1e-300 / 2, lower = FALSE) is about 37, so |z| = 50 cannot be
  ## adjusted and is returned as is.
  expect_equal(fiqt(c(50, 1)), c(50, 1))
})

test_that("min.p is honoured", {
  ## A tighter bound brings the pass-through threshold down to about 5.7.
  expect_equal(fiqt(c(8, 1), min.p = 1e-8)[1], 8)
  expect_lt(fiqt(c(8, 1), min.p = 1e-300)[1], 8)
})

test_that("missing values stay missing", {
  adj <- fiqt(c(1, NA, 3))

  expect_length(adj, 3L)
  expect_true(is.na(adj[2]))
  expect_true(all(is.finite(adj[c(1, 3)])))
})

test_that("an empty input gives an empty result", {
  adj <- fiqt(numeric(0))

  expect_type(adj, "double")
  expect_length(adj, 0L)
})
