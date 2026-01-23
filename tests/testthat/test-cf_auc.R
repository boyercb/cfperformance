# Tests for cf_auc function

test_that("cf_auc returns correct structure", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr",
    se_method = "none"
  )

  expect_s3_class(result, "cf_auc")
  expect_s3_class(result, "cf_performance")

  expect_true("estimate" %in% names(result))
  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_equal(result$metric, "auc")
})


test_that("cf_auc naive equals standard AUC", {
  set.seed(42)
  n <- 500
  y <- rbinom(n, 1, 0.3)
  pred <- runif(n)
  a <- rbinom(n, 1, 0.5)
  x <- rnorm(n)

  result <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    estimator = "naive",
    se_method = "none"
  )

  # Calculate standard AUC
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(c(pred[y == 1], pred[y == 0]))
  expected_auc <- (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)

  expect_equal(result$estimate, expected_auc)
})


test_that("cf_auc requires binary outcomes", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  y <- rnorm(n)  # Continuous outcomes
  pred <- runif(n)

  expect_error(
    cf_auc(predictions = pred, outcomes = y, treatment = a,
           covariates = data.frame(x = x), se_method = "none"),
    "binary"
  )
})


test_that("cf_auc different estimators work", {
  set.seed(123)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  estimators <- c("naive", "cl", "ipw", "dr")
  results <- lapply(estimators, function(est) {
    cf_auc(
      predictions = pred, outcomes = y, treatment = a,
      covariates = data.frame(x = x), estimator = est, se_method = "none"
    )
  })

  # All should return valid AUC values
  for (r in results) {
    expect_true(r$estimate >= 0 && r$estimate <= 1)
  }
})


test_that("cf_auc cross_fit with ML learners works", {
  skip_if_not_installed("ranger")

  set.seed(42)
  n <- 200
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr",
    propensity_model = ml_learner("ranger", num.trees = 50),
    outcome_model = ml_learner("ranger", num.trees = 50),
    cross_fit = TRUE,
    n_folds = 3,
    se_method = "influence"
  )

  expect_s3_class(result, "cf_auc")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_true(!is.null(result$se))
  expect_true(result$se > 0)
})


test_that("cf_auc cross_fit only works with DR estimator",
{
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Cross-fitting with DR should work
  result_dr <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr",
    cross_fit = TRUE,
    n_folds = 3,
    se_method = "none"
  )
  expect_s3_class(result_dr, "cf_auc")

  # Cross-fitting with other estimators falls back to non-cross-fit
  result_ipw <- cf_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "ipw",
    cross_fit = TRUE,
    n_folds = 3,
    se_method = "none"
  )
  expect_s3_class(result_ipw, "cf_auc")
})
