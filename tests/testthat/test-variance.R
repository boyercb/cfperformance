# Tests for variance estimation functions

test_that("influence function SE works for naive MSE estimator", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "naive",
    se_method = "influence"
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})

test_that("influence function SE works for DR MSE estimator", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "influence"
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})

test_that("bootstrap SE works for MSE estimators", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50  # Small for testing
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})

test_that("cross-fitting works for DR MSE", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result_cf <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "influence",
    cross_fit = TRUE,
    n_folds = 5
  )

  result_nocf <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "influence",
    cross_fit = FALSE
  )

  # Both should produce valid estimates
  expect_true(!is.na(result_cf$estimate))
  expect_true(!is.na(result_nocf$estimate))
  expect_true(!is.na(result_cf$se))
  expect_true(!is.na(result_nocf$se))

  # Estimates should be similar (both are consistent)
  expect_equal(result_cf$estimate, result_nocf$estimate, tolerance = 0.1)
})

test_that("influence function SE works for AUC estimators", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "naive",
    se_method = "influence"
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
})

test_that("bootstrap SE works for AUC estimators", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50
  )

  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
})

test_that("confint method returns correct format", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "influence"
  )

  ci <- confint(result)

  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_equal(ci[1, 1], result$ci_lower)
  expect_equal(ci[1, 2], result$ci_upper)
})

test_that("coef method returns estimate", {
  set.seed(123)
  n <- 200

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "none"
  )

  est <- coef(result)
  expect_equal(est, result$estimate)
})

test_that("SE methods are comparable", {
  skip_on_cran()
  skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "true",
              "Skipping slow comparison test")

  set.seed(123)
  n <- 500

  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.3 * a))
  pred <- plogis(-1 + 0.8 * x)
  covariates <- data.frame(x = x)

  result_boot <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 200
  )

  result_inf <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covariates,
    treatment_level = 0,
    estimator = "dr",
    se_method = "influence"
  )

  # SEs should be of similar magnitude
  # (allowing for bootstrap variability)
  ratio <- result_boot$se / result_inf$se
  expect_true(ratio > 0.5 && ratio < 2,
              info = sprintf("SE ratio: %.2f", ratio))
})
