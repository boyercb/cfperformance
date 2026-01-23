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


# ==============================================================================
# Bootstrap Standard Error Tests
# ==============================================================================

test_that("cf_calibration bootstrap returns SE and CI", {
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
    se_method = "bootstrap",
    n_boot = 50  # Small for testing
  )

  expect_s3_class(result, "cf_calibration")
  expect_equal(result$se_method, "bootstrap")

  # Check SE structure
  expect_true(!is.null(result$se))
  expect_true("ici" %in% names(result$se))
  expect_true("e50" %in% names(result$se))
  expect_true("emax" %in% names(result$se))

  # Check CI structure
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))

  # SEs should be positive

  expect_true(result$se$ici > 0)
  expect_true(result$se$e50 > 0)
})


test_that("cf_calibration bootstrap returns curves for CI bands", {
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
    se_method = "bootstrap",
    n_boot = 50
  )

  expect_true(!is.null(result$boot_curves))
  expect_true(is.data.frame(result$boot_curves))
  expect_true("predicted" %in% names(result$boot_curves))
  expect_true("ci_lower" %in% names(result$boot_curves))
  expect_true("ci_upper" %in% names(result$boot_curves))

  # CI bounds should be reasonable
  expect_true(all(result$boot_curves$ci_lower <= result$boot_curves$ci_upper, na.rm = TRUE))
})


test_that("cf_calibration without bootstrap has NULL SE fields", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + x)

  result <- cf_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  # Use [["se"]] to avoid partial matching with se_method
  expect_null(result[["se"]])
  expect_null(result[["ci_lower"]])
  expect_null(result[["ci_upper"]])
  expect_null(result[["boot_curves"]])
})
