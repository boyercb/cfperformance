# test-roc.R
# Tests for ROC curve functions

# =============================================================================
# Test Data Setup
# =============================================================================

test_that("test data generation works", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  expect_true(length(x) == n)
  expect_true(sum(s == 1) > 0)
  expect_true(sum(s == 0) > 0)
})

# =============================================================================
# tr_roc Tests
# =============================================================================

test_that("tr_roc returns valid structure", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  expect_s3_class(roc, "tr_roc")
  expect_s3_class(roc, "roc_curve")
  expect_true("thresholds" %in% names(roc))
  expect_true("sensitivity" %in% names(roc))
  expect_true("fpr" %in% names(roc))
  expect_true("auc" %in% names(roc))
})

test_that("tr_roc computes valid AUC", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 21
  )

  # AUC should be between 0 and 1
  expect_true(roc$auc >= 0 && roc$auc <= 1)
  expect_true(roc$naive_auc >= 0 && roc$naive_auc <= 1)
})

test_that("tr_roc works with custom thresholds", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  custom_thresholds <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    thresholds = custom_thresholds
  )

  expect_equal(roc$thresholds, custom_thresholds)
  expect_length(roc$sensitivity, length(custom_thresholds))
})

test_that("tr_roc works with different estimators", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  for (est in c("dr", "om", "ipw", "naive")) {
    roc <- tr_roc(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = data.frame(x = x),
      estimator = est,
      n_thresholds = 11
    )

    expect_s3_class(roc, "tr_roc")
    expect_equal(roc$estimator, est)
    expect_true(roc$auc >= 0 && roc$auc <= 1)
  }
})

test_that("tr_roc can exclude naive curve", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 11,
    include_naive = FALSE
  )

  expect_null(roc$naive_sensitivity)
  expect_null(roc$naive_fpr)
  expect_null(roc$naive_auc)
})

# =============================================================================
# cf_roc Tests
# =============================================================================

test_that("cf_roc returns valid structure", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  expect_s3_class(roc, "cf_roc")
  expect_s3_class(roc, "roc_curve")
  expect_true("thresholds" %in% names(roc))
  expect_true("sensitivity" %in% names(roc))
  expect_true("fpr" %in% names(roc))
  expect_true("auc" %in% names(roc))
})

test_that("cf_roc computes valid AUC", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    n_thresholds = 21
  )

  expect_true(roc$auc >= 0 && roc$auc <= 1)
  expect_true(roc$naive_auc >= 0 && roc$naive_auc <= 1)
})

test_that("cf_roc works with different estimators", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  for (est in c("dr", "cl", "ipw", "naive")) {
    roc <- cf_roc(
      predictions = pred,
      outcomes = y,
      treatment = a,
      covariates = data.frame(x = x),
      estimator = est,
      n_thresholds = 11
    )

    expect_s3_class(roc, "cf_roc")
    expect_equal(roc$estimator, est)
    expect_true(roc$auc >= 0 && roc$auc <= 1)
  }
})

# =============================================================================
# AUC Computation Tests
# =============================================================================

test_that("AUC of random classifier is approximately 0.5", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.3)  # Random outcome
  pred <- runif(n)  # Random predictions

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    estimator = "naive",
    n_thresholds = 101
  )

  # Should be close to 0.5 for random predictions
  expect_true(abs(roc$naive_auc - 0.5) < 0.15)
})

test_that("AUC of perfect classifier is approximately 1", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + 2 * x))  # Strongly x-dependent
  pred <- plogis(-1 + 2 * x)  # Perfect predictions

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    estimator = "naive",
    n_thresholds = 101
  )

  # Should be high for good predictions
  expect_true(roc$naive_auc > 0.7)
})

# =============================================================================
# as.data.frame Tests
# =============================================================================

test_that("as.data.frame.tr_roc works correctly", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  df <- as.data.frame(roc)

  expect_s3_class(df, "data.frame")
  expect_true("threshold" %in% names(df))
  expect_true("fpr" %in% names(df))
  expect_true("sensitivity" %in% names(df))
  expect_true("type" %in% names(df))

  # Should have both adjusted and naive rows
  expect_true("adjusted" %in% df$type)
  expect_true("naive" %in% df$type)
})

test_that("as.data.frame.cf_roc works correctly", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  df <- as.data.frame(roc)

  expect_s3_class(df, "data.frame")
  expect_true(nrow(df) == 22)  # 11 adjusted + 11 naive
})

# =============================================================================
# Print Method Tests
# =============================================================================

test_that("print method works for tr_roc", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  expect_output(print(roc), "Transportable ROC")
  expect_output(print(roc), "AUC:")
})

test_that("print method works for cf_roc", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  expect_output(print(roc), "Counterfactual ROC")
  expect_output(print(roc), "AUC:")
})

# =============================================================================
# Plot Method Tests
# =============================================================================

test_that("plot method works for tr_roc", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- tr_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  # Should not error
  expect_no_error(plot(roc))
})

test_that("plot method works for cf_roc", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    n_thresholds = 11
  )

  expect_no_error(plot(roc))
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("ROC handles extreme thresholds", {
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  roc <- cf_roc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    thresholds = c(0, 0.5, 1),
    estimator = "naive"
  )

  # At threshold 0, sensitivity = 1 (all predicted positive)
  expect_equal(roc$naive_sensitivity[1], 1)

  # At threshold 1, sensitivity = 0 (none predicted positive)
  expect_equal(roc$naive_sensitivity[3], 0)
})
