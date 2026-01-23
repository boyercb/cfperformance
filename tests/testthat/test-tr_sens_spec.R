# test-tr_sensitivity.R
# Tests for transportable sensitivity, specificity, and FPR functions

# =============================================================================
# Test Data Setup
# =============================================================================

test_that("test data generation works", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  expect_true(length(x) == n)
  expect_true(sum(s == 1) > 0)
  expect_true(sum(s == 0) > 0)
})

# =============================================================================
# tr_sensitivity Tests
# =============================================================================

test_that("tr_sensitivity returns valid structure", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  # Check class

expect_s3_class(result, "tr_sensitivity")
  expect_s3_class(result, "tr_performance")

  # Check required fields
  expect_true("estimate" %in% names(result))
  expect_true("threshold" %in% names(result))
  expect_true("naive_estimate" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_true("estimator" %in% names(result))
  expect_true("analysis" %in% names(result))

  # Check metric name
  expect_equal(result$metric, "sensitivity")
})

test_that("tr_sensitivity estimate is valid probability", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_true(result$naive_estimate >= 0 && result$naive_estimate <= 1)
})

test_that("tr_sensitivity works with different estimators", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  estimators <- c("dr", "om", "ipw", "naive")

  for (est in estimators) {
    result <- tr_sensitivity(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = covars,
      estimator = est,
      se_method = "none"
    )

    expect_s3_class(result, "tr_sensitivity")
    expect_equal(result$estimator, est)
    expect_true(result$estimate >= 0 && result$estimate <= 1)
  }
})

test_that("tr_sensitivity supports vectorized threshold", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.2, 0.3, 0.5, 0.7)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = thresholds,
    se_method = "none"
  )

  expect_length(result$estimate, length(thresholds))
  expect_length(result$threshold, length(thresholds))
  expect_length(result$naive_estimate, length(thresholds))

  # Sensitivity should decrease as threshold increases
  for (i in 2:length(thresholds)) {
    expect_true(result$estimate[i] <= result$estimate[i - 1] + 0.01)  # small tolerance
  }
})

test_that("tr_sensitivity threshold at 0 gives sensitivity of 1",
{
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0,  # All predictions > 0
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result$naive_estimate, 1, tolerance = 0.001)
})

# =============================================================================
# tr_specificity Tests
# =============================================================================

test_that("tr_specificity returns valid structure", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  expect_s3_class(result, "tr_specificity")
  expect_s3_class(result, "tr_performance")
  expect_equal(result$metric, "specificity")
})

test_that("tr_specificity estimate is valid probability", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_true(result$naive_estimate >= 0 && result$naive_estimate <= 1)
})

test_that("tr_specificity increases as threshold increases", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.2, 0.3, 0.5, 0.7)

  result <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
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

test_that("tr_specificity threshold at 1 gives specificity of 1", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 1,  # All predictions <= 1
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result$naive_estimate, 1, tolerance = 0.001)
})

# =============================================================================
# tr_fpr Tests
# =============================================================================

test_that("tr_fpr is complement of tr_specificity", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result_spec <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    se_method = "none"
  )

  result_fpr <- tr_fpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    se_method = "none"
  )

  expect_equal(result_fpr$estimate, 1 - result_spec$estimate, tolerance = 1e-10)
  expect_equal(result_fpr$naive_estimate, 1 - result_spec$naive_estimate, tolerance = 1e-10)
  expect_equal(result_fpr$metric, "fpr")
})

test_that("tr_tpr alias works like tr_sensitivity", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result_sens <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    se_method = "none"
  )

  result_tpr <- tr_tpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    se_method = "none"
  )

  expect_equal(result_tpr$estimate, result_sens$estimate, tolerance = 1e-10)
})

# =============================================================================
# Bootstrap SE Tests
# =============================================================================

test_that("tr_sensitivity bootstrap returns SE and CI", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0.3,
    se_method = "bootstrap",
    n_boot = 30  # Small for speed
  )

  expect_true(!is.null(result$se))
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})

test_that("tr_specificity bootstrap returns SE and CI", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_specificity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
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

test_that("tr_sensitivity bootstrap works with multiple thresholds", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  thresholds <- c(0.3, 0.5, 0.7)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
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

test_that("tr_sensitivity handles edge case thresholds", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Threshold at 0 - all positive predictions
  result_0 <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0,
    estimator = "naive",
    se_method = "none"
  )

  # Threshold at 1 - no positive predictions
  result_1 <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 1,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result_0$naive_estimate, 1, tolerance = 0.001)
  expect_equal(result_1$naive_estimate, 0, tolerance = 0.001)
})

test_that("tr_sensitivity works with different treatment levels", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result_a0 <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    se_method = "none"
  )

  result_a1 <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 1,
    se_method = "none"
  )

  expect_s3_class(result_a0, "tr_sensitivity")
  expect_s3_class(result_a1, "tr_sensitivity")
  expect_equal(result_a0$treatment_level, 0)
  expect_equal(result_a1$treatment_level, 1)
})

# =============================================================================
# Analysis Type Tests
# =============================================================================

test_that("tr_sensitivity works with joint analysis", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    analysis = "joint",
    se_method = "none"
  )

  expect_s3_class(result, "tr_sensitivity")
  expect_equal(result$analysis, "joint")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("transport and joint analysis give different estimates", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))  # non-uniform selection
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result_transport <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    analysis = "transport",
    se_method = "none"
  )

  result_joint <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    analysis = "joint",
    se_method = "none"
  )

  # With differential selection, estimates should differ
  expect_true(result_transport$estimate != result_joint$estimate)
})

# =============================================================================
# Print Method Tests
# =============================================================================

test_that("print method works for single threshold", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = 0.5,
    se_method = "none"
  )

  expect_output(print(result), "Transportable Sensitivity")
  expect_output(print(result), "Threshold:")
  expect_output(print(result), "Estimate:")
})

test_that("print method works for multiple thresholds", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_sensitivity(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    threshold = c(0.3, 0.5, 0.7),
    se_method = "none"
  )

  expect_output(print(result), "Results by threshold")
})

test_that("print method works for FPR", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_fpr(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  expect_output(print(result), "False Positive Rate")
})
