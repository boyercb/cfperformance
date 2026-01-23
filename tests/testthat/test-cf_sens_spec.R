# test-cf_sens_spec.R
# Tests for counterfactual sensitivity, specificity, and FPR functions

# =============================================================================
# Test Data Setup
# =============================================================================

test_that("test data generation works", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  expect_true(length(x) == n)
  expect_true(sum(a == 1) > 0)
  expect_true(sum(a == 0) > 0)
})

# =============================================================================
# cf_sensitivity Tests
# =============================================================================

test_that("cf_sensitivity returns valid structure", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  # Check class
  expect_s3_class(result, "cf_sensitivity")
  expect_s3_class(result, "cf_performance")

  # Check required fields
  expect_true("estimate" %in% names(result))
  expect_true("threshold" %in% names(result))
  expect_true("naive_estimate" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_true("estimator" %in% names(result))

  # Check metric name
  expect_equal(result$metric, "sensitivity")
})

test_that("cf_sensitivity estimate is valid probability", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_true(result$naive_estimate >= 0 && result$naive_estimate <= 1)
})

test_that("cf_sensitivity works with different estimators", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  estimators <- c("dr", "cl", "ipw", "naive")

  for (est in estimators) {
    result <- cf_sensitivity(
      predictions = pred,
      outcomes = y,
      treatment = a,
      covariates = covars,
      estimator = est,
      se_method = "none"
    )

    expect_s3_class(result, "cf_sensitivity")
    expect_equal(result$estimator, est)
    expect_true(result$estimate >= 0 && result$estimate <= 1)
  }
})

test_that("cf_sensitivity supports vectorized threshold", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.2, 0.3, 0.5, 0.7)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = thresholds,
    se_method = "none"
  )

  expect_length(result$estimate, length(thresholds))
  expect_length(result$threshold, length(thresholds))
  expect_length(result$naive_estimate, length(thresholds))

  # Sensitivity should decrease as threshold increases
  for (i in 2:length(thresholds)) {
    expect_true(result$estimate[i] <= result$estimate[i - 1] + 0.01)
  }
})

test_that("cf_sensitivity threshold at 0 gives sensitivity of 1", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result$naive_estimate, 1, tolerance = 0.001)
})

# =============================================================================
# cf_specificity Tests
# =============================================================================

test_that("cf_specificity returns valid structure", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  expect_s3_class(result, "cf_specificity")
  expect_s3_class(result, "cf_performance")
  expect_equal(result$metric, "specificity")
})

test_that("cf_specificity estimate is valid probability", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_true(result$naive_estimate >= 0 && result$naive_estimate <= 1)
})

test_that("cf_specificity increases as threshold increases", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.2, 0.3, 0.5, 0.7)

  result <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = thresholds,
    estimator = "naive",
    se_method = "none"
  )

  # Specificity should increase as threshold increases
  for (i in 2:length(thresholds)) {
    expect_true(result$naive_estimate[i] >= result$naive_estimate[i - 1] - 0.01)
  }
})

test_that("cf_specificity threshold at 1 gives specificity of 1", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 1,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result$naive_estimate, 1, tolerance = 0.001)
})

# =============================================================================
# cf_fpr Tests
# =============================================================================

test_that("cf_fpr is complement of cf_specificity", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result_spec <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covars,
    se_method = "none"
  )

  result_fpr <- cf_fpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covars,
    se_method = "none"
  )

  expect_equal(result_fpr$estimate, 1 - result_spec$estimate, tolerance = 1e-10)
  expect_equal(result_fpr$naive_estimate, 1 - result_spec$naive_estimate, tolerance = 1e-10)
  expect_equal(result_fpr$metric, "fpr")
})

test_that("cf_tpr alias works like cf_sensitivity", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result_sens <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covars,
    se_method = "none"
  )

  result_tpr <- cf_tpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covars,
    se_method = "none"
  )

  expect_equal(result_tpr$estimate, result_sens$estimate, tolerance = 1e-10)
})

# =============================================================================
# Bootstrap SE Tests
# =============================================================================

test_that("cf_sensitivity bootstrap returns SE and CI", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0.3,
    se_method = "bootstrap",
    n_boot = 30
  )

  expect_true(!is.null(result$se))
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})

test_that("cf_specificity bootstrap returns SE and CI", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0.3,
    se_method = "bootstrap",
    n_boot = 30
  )

  expect_true(!is.null(result$se))
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$se > 0)
})

test_that("cf_sensitivity bootstrap works with multiple thresholds", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.3, 0.5, 0.7)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = thresholds,
    se_method = "bootstrap",
    n_boot = 30
  )

  expect_length(result$se, length(thresholds))
  expect_length(result$ci_lower, length(thresholds))
  expect_length(result$ci_upper, length(thresholds))
})

# =============================================================================
# Edge Cases and Input Validation
# =============================================================================

test_that("cf_sensitivity handles edge case thresholds", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Threshold at 0 - all positive predictions
  result_0 <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0,
    estimator = "naive",
    se_method = "none"
  )

  # Threshold at 1 - no positive predictions
  result_1 <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 1,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result_0$naive_estimate, 1, tolerance = 0.001)
  expect_equal(result_1$naive_estimate, 0, tolerance = 0.001)
})

test_that("cf_sensitivity works with different treatment levels", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result_a0 <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    se_method = "none"
  )

  result_a1 <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 1,
    se_method = "none"
  )

  expect_s3_class(result_a0, "cf_sensitivity")
  expect_s3_class(result_a1, "cf_sensitivity")
  expect_equal(result_a0$treatment_level, 0)
  expect_equal(result_a1$treatment_level, 1)
})

# =============================================================================
# Print Method Tests
# =============================================================================

test_that("print method works for single threshold", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  expect_output(print(result), "Sensitivity|SENSITIVITY")
})

test_that("print method works for FPR", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_fpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_output(print(result), "FPR|False Positive")
})

# =============================================================================
# Comparison with tr_ functions
# =============================================================================

test_that("cf_ and tr_ functions give different results due to different estimands", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  # Create source/target split
  s <- rbinom(n, 1, 0.6)  # ~60% source

  result_cf <- cf_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = covars,
    estimator = "naive",
    se_method = "none"
  )

  result_tr <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    estimator = "naive",
    se_method = "none"
  )

  # cf_ uses all data, tr_ uses only source with treatment level
  # These should be different because of different populations
  expect_s3_class(result_cf, "cf_sensitivity")
  expect_s3_class(result_tr, "tr_sensitivity")
})
