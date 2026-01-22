# Tests for cf_calibration function

test_that("cf_calibration returns correct structure", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.9 * x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0
  )

  expect_s3_class(result, "cf_calibration")
  expect_s3_class(result, "cf_performance")

  expect_true("predicted" %in% names(result))
  expect_true("observed" %in% names(result))
  expect_true("ici" %in% names(result))
  expect_true("e50" %in% names(result))
  expect_true("emax" %in% names(result))
})


test_that("cf_calibration metrics are non-negative", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x)
  )

  expect_true(result$ici >= 0)
  expect_true(result$e50 >= 0)
  expect_true(result$e90 >= 0)
  expect_true(result$emax >= 0)
})


test_that("cf_calibration binned smoother works", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    smoother = "binned",
    n_bins = 5
  )

  expect_equal(result$smoother, "binned")
  expect_equal(length(result$predicted), 5)
})


test_that("cf_calibration requires binary outcomes", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rnorm(n)  # Continuous
  pred <- runif(n)

  expect_error(
    cf_calibration(predictions = pred, outcomes = y, treatment = a,
                   covariates = data.frame(x = x)),
    "binary"
  )
})
