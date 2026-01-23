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


test_that("cf_calibration CL estimator works", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.9 * x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "cl"
  )

  expect_s3_class(result, "cf_calibration")
  expect_equal(result$estimator, "cl")
  expect_equal(result$n_obs, n)  # CL uses all observations
  expect_true(result$ici >= 0)
})


test_that("cf_calibration DR estimator works", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.9 * x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr"
  )

  expect_s3_class(result, "cf_calibration")
  expect_equal(result$estimator, "dr")
  expect_equal(result$n_obs, n)  # DR uses all observations
  expect_true(result$ici >= 0)
  expect_true(!is.null(result$propensity_model))
  expect_true(!is.null(result$outcome_model))
})


test_that("cf_calibration DR is default estimator", {
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

  expect_equal(result$estimator, "dr")
})


test_that("cf_calibration all three estimators produce reasonable results", {
  set.seed(123)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covs <- data.frame(x = x)

  result_ipw <- cf_calibration(pred, y, a, covs, treatment_level = 0, estimator = "ipw")
  result_cl <- cf_calibration(pred, y, a, covs, treatment_level = 0, estimator = "cl")
  result_dr <- cf_calibration(pred, y, a, covs, treatment_level = 0, estimator = "dr")

  # All should produce valid ICI values
  expect_true(result_ipw$ici >= 0 && result_ipw$ici <= 1)
  expect_true(result_cl$ici >= 0 && result_cl$ici <= 1)
  expect_true(result_dr$ici >= 0 && result_dr$ici <= 1)

  # Results should be somewhat similar (within a reasonable range)
  # This is a loose check since estimators can vary
  ici_values <- c(result_ipw$ici, result_cl$ici, result_dr$ici)
  expect_true(max(ici_values) - min(ici_values) < 0.2)
})
