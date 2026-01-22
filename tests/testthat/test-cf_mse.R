# Tests for cf_mse function

test_that("cf_mse returns correct structure", {
  # Generate test data
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr",
    se_method = "none"
  )

  # Check class

expect_s3_class(result, "cf_mse")
  expect_s3_class(result, "cf_performance")

  # Check required elements
  expect_true("estimate" %in% names(result))
  expect_true("estimator" %in% names(result))
  expect_true("treatment_level" %in% names(result))
  expect_true("n_obs" %in% names(result))
  expect_true("naive_estimate" %in% names(result))

  # Check values are reasonable
  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_equal(result$estimator, "dr")
  expect_equal(result$n_obs, n)
})


test_that("cf_mse naive estimator works", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)  # Random treatment
  y <- rbinom(n, 1, plogis(-1 + x))  # No treatment effect
  pred <- plogis(-1 + x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    estimator = "naive",
    se_method = "none"
  )

  # Naive estimate should equal mean squared error
  expected_mse <- mean((y - pred)^2)
  expect_equal(result$estimate, expected_mse)
})


test_that("cf_mse different estimators produce different results with confounding", {
  set.seed(123)
  n <- 500
  x <- rnorm(n)
  # Strong confounding
  a <- rbinom(n, 1, plogis(-1 + 2 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 1 * a))
  pred <- plogis(-1 + 0.8 * x)

  naive_result <- cf_mse(
    predictions = pred, outcomes = y, treatment = a,
    covariates = data.frame(x = x), estimator = "naive", se_method = "none"
  )

  dr_result <- cf_mse(
    predictions = pred, outcomes = y, treatment = a,
    covariates = data.frame(x = x), estimator = "dr", se_method = "none"
  )

  # With confounding, estimates should differ
  expect_false(isTRUE(all.equal(naive_result$estimate, dr_result$estimate,
                                 tolerance = 0.001)))
})


test_that("cf_mse input validation works", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.3)
  pred <- runif(n)

  # Wrong length predictions
  expect_error(
    cf_mse(predictions = pred[1:50], outcomes = y, treatment = a,
           covariates = data.frame(x = x), se_method = "none"),
    "same length"
  )

  # Non-binary treatment
  expect_error(
    cf_mse(predictions = pred, outcomes = y, treatment = c(a[-1], 2),
           covariates = data.frame(x = x), se_method = "none"),
    "binary"
  )
})


test_that("cf_mse bootstrap SE works", {
  skip_on_cran()  # Skip slow test on CRAN

  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50  # Small number for testing
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$ci_lower < result$ci_upper)
})
