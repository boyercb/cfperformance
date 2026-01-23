# Tests for tr_mse() - Transportable MSE estimation

test_that("tr_mse works with basic inputs", {
  set.seed(123)
  n <- 500

  # Generate data
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))  # Source indicator
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),  # Randomized in source
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))  # Confounded in target
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Test naive estimator
  result_naive <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "transport",
    estimator = "naive",
    se_method = "none"
  )

  expect_s3_class(result_naive, "tr_mse")
  expect_s3_class(result_naive, "tr_performance")
  expect_true(!is.na(result_naive$estimate))
  expect_equal(result_naive$estimator, "naive")
  expect_equal(result_naive$analysis, "transport")
})


test_that("tr_mse works with all estimators for transport analysis", {
  set.seed(456)
  n <- 500

  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  covariates <- data.frame(x = x)

  for (est in c("naive", "om", "ipw", "dr")) {
    result <- tr_mse(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = covariates,
      treatment_level = 0,
      analysis = "transport",
      estimator = est,
      se_method = "none"
    )

    expect_s3_class(result, "tr_mse")
    expect_true(!is.na(result$estimate), info = paste("Estimator:", est))
    expect_equal(result$estimator, est)
  }
})


test_that("tr_mse works with all estimators for joint analysis", {
  set.seed(789)
  n <- 500

  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  covariates <- data.frame(x = x)

  for (est in c("naive", "om", "ipw", "dr")) {
    result <- tr_mse(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = covariates,
      treatment_level = 0,
      analysis = "joint",
      estimator = est,
      se_method = "none"
    )

    expect_s3_class(result, "tr_mse")
    expect_true(!is.na(result$estimate), info = paste("Estimator:", est))
    expect_equal(result$analysis, "joint")
  }
})


test_that("tr_mse bootstrap SE works", {
  set.seed(111)
  n <- 300

  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50  # Small for speed
  )

  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
  expect_true(result$ci_lower < result$estimate)
  expect_true(result$ci_upper > result$estimate)
})


test_that("tr_mse validates inputs correctly", {
  set.seed(222)
  n <- 100
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)
  pred <- runif(n)
  covariates <- data.frame(x = x)

  # Mismatched lengths
  expect_error(
    tr_mse(pred[1:50], y, a, s, covariates, se_method = "none"),
    "same length"
  )

  # Invalid treatment
  expect_error(
    tr_mse(pred, y, rep(2, n), s, covariates, se_method = "none"),
    "binary"
  )

  # Invalid source
  expect_error(
    tr_mse(pred, y, a, rep(2, n), covariates, se_method = "none"),
    "binary"
  )

  # No target observations
  expect_error(
    tr_mse(pred, y, a, rep(1, n), covariates, se_method = "none"),
    "No target"
  )

  # No source observations
  expect_error(
    tr_mse(pred, y, a, rep(0, n), covariates, se_method = "none"),
    "No source"
  )
})


test_that("tr_mse print and summary methods work", {
  set.seed(333)
  n <- 200

  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)
  pred <- runif(n)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  # Test print
  expect_output(print(result), "Transportable MSE")
  expect_output(print(result), "transport")
  expect_output(print(result), "dr")

  # Test summary
  expect_output(summary(result), "Summary")
  expect_output(summary(result), "Transportable")

  # Test coef
  expect_equal(coef(result), result$estimate)
})


test_that("tr_mse treatment_level = 1 works", {
  set.seed(444)
  n <- 300

  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + x + 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 1,
    estimator = "dr",
    se_method = "none"
  )

  expect_equal(result$treatment_level, 1)
  expect_true(!is.na(result$estimate))
})


test_that("tr_mse returns correct structure", {
  set.seed(555)
  n <- 200

  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)
  pred <- runif(n)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    se_method = "none"
  )

  # Check all expected elements are present

  expect_true("estimate" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("estimator" %in% names(result))
  expect_true("analysis" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_true("treatment_level" %in% names(result))
  expect_true("n_target" %in% names(result))
  expect_true("n_source" %in% names(result))
  expect_true("naive_estimate" %in% names(result))
  expect_true("call" %in% names(result))

  expect_equal(result$metric, "mse")
  expect_equal(result$n_target + result$n_source, n)
})


# =============================================================================
# Tests for influence function SE
# =============================================================================

test_that("tr_mse influence SE works for DR transport", {
  skip_on_cran()

  set.seed(1234)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "influence"
  )

  expect_s3_class(result, "tr_mse")
  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
})

test_that("tr_mse influence SE works for DR joint", {
  skip_on_cran()

  set.seed(2345)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "joint",
    estimator = "dr",
    se_method = "influence"
  )

  expect_s3_class(result, "tr_mse")
  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
})

test_that("tr_mse influence SE works for all estimators", {
  skip_on_cran()

  set.seed(3456)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  for (est in c("om", "ipw", "dr")) {
    result <- tr_mse(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = data.frame(x = x),
      treatment_level = 0,
      analysis = "transport",
      estimator = est,
      se_method = "influence"
    )

    expect_true(!is.na(result$se), info = paste("Estimator:", est))
    expect_true(result$se > 0, info = paste("Estimator:", est))
  }
})


# =============================================================================
# Tests for stratified bootstrap
# =============================================================================

test_that("tr_mse stratified bootstrap works", {
  skip_on_cran()

  set.seed(4567)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Stratified bootstrap (default)
  result_strat <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50,
    stratified_boot = TRUE
  )

  expect_s3_class(result_strat, "tr_mse")
  expect_true(!is.na(result_strat$se))
  expect_true(result_strat$se > 0)
})

test_that("tr_mse non-stratified bootstrap works", {
  skip_on_cran()

  set.seed(5678)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)

  # Non-stratified bootstrap
  result_simple <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50,
    stratified_boot = FALSE
  )

  expect_s3_class(result_simple, "tr_mse")
  expect_true(!is.na(result_simple$se))
  expect_true(result_simple$se > 0)
})
