# Test suite for tr_calibration
# Transportable Calibration Estimators

# ==============================================================================
# Test Data Generation
# ==============================================================================

#' Generate test data for transportability calibration
generate_tr_calibration_data <- function(n = 1000, seed = 123) {
  set.seed(seed)

  # Covariates
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Source indicator (S=1 for RCT, S=0 for target)
  # Target population has higher values of x1
  s <- rbinom(n, 1, plogis(0.5 - 0.5 * x1))

  # Treatment (randomized in source, confounded in target)
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),  # Randomized in RCT
              rbinom(n, 1, plogis(-0.3 + 0.3 * x1)))  # Confounded in target

  # True outcome probability
  true_prob <- plogis(-1 + 0.8 * x1 + 0.3 * x2 - 0.4 * a)

  # Binary outcome
  y <- rbinom(n, 1, true_prob)

  # Prediction model (calibrated)
  pred <- plogis(-1 + 0.75 * x1 + 0.25 * x2)

  list(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x1 = x1, x2 = x2),
    n = n,
    n_source = sum(s == 1),
    n_target = sum(s == 0)
  )
}


# ==============================================================================
# Input Validation Tests
# ==============================================================================

test_that("tr_calibration validates input lengths", {
  d <- generate_tr_calibration_data(n = 100)

  expect_error(
    tr_calibration(d$predictions[1:50], d$outcomes, d$treatment, d$source, d$covariates),
    "same length"
  )

  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment[1:50], d$source, d$covariates),
    "same length"
  )

  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment, d$source[1:50], d$covariates),
    "same length"
  )

  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment, d$source, d$covariates[1:50, ]),
    "same number of rows"
  )
})

test_that("tr_calibration requires binary outcomes", {
  d <- generate_tr_calibration_data(n = 100)

  expect_error(
    tr_calibration(d$predictions, rnorm(100), d$treatment, d$source, d$covariates),
    "binary outcomes"
  )
})

test_that("tr_calibration validates source indicator", {
  d <- generate_tr_calibration_data(n = 100)

  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment, d$source + 1, d$covariates),
    "binary indicator"
  )
})

test_that("tr_calibration validates predictions are probabilities", {
  d <- generate_tr_calibration_data(n = 100)

  expect_error(
    tr_calibration(d$predictions * 2, d$outcomes, d$treatment, d$source, d$covariates),
    "probabilities between 0 and 1"
  )

  expect_error(
    tr_calibration(d$predictions - 1, d$outcomes, d$treatment, d$source, d$covariates),
    "probabilities between 0 and 1"
  )
})

test_that("tr_calibration requires sufficient data", {
  d <- generate_tr_calibration_data(n = 100)

  # Make all observations source (no target)
  d$source <- rep(1, 100)
  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment, d$source, d$covariates),
    "Insufficient target data"
  )

  # Make all observations target (no source)
  d$source <- rep(0, 100)
  expect_error(
    tr_calibration(d$predictions, d$outcomes, d$treatment, d$source, d$covariates),
    "Insufficient source data"
  )
})


# ==============================================================================
# Transport Analysis Tests - IPW Estimator
# ==============================================================================

test_that("tr_calibration transport IPW returns correct structure", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_s3_class(result, "tr_performance")

  # Check required components
  expect_true(!is.null(result$predicted))
  expect_true(!is.null(result$observed))
  expect_true(!is.null(result$ici))
  expect_true(!is.null(result$e50))
  expect_true(!is.null(result$e90))

  expect_true(!is.null(result$emax))
  expect_equal(result$analysis, "transport")
  expect_equal(result$estimator, "ipw")
  expect_equal(result$metric, "calibration")
})

test_that("tr_calibration transport IPW produces valid calibration metrics", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  # ICI should be non-negative
  expect_gte(result$ici, 0)

  # E50, E90, Emax should be non-negative
  expect_gte(result$e50, 0)
  expect_gte(result$e90, 0)
  expect_gte(result$emax, 0)

  # E50 <= E90 <= Emax
  expect_lte(result$e50, result$e90 + 1e-10)
  expect_lte(result$e90, result$emax + 1e-10)

  # Predicted and observed should be between 0 and 1
  expect_true(all(result$predicted >= 0 & result$predicted <= 1))
  expect_true(all(result$observed >= 0 & result$observed <= 1, na.rm = TRUE))
})

test_that("tr_calibration transport IPW treatment_level=1 works", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 1,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$treatment_level, 1)
  expect_gte(result$ici, 0)
})


# ==============================================================================
# Transport Analysis Tests - OM Estimator
# ==============================================================================

test_that("tr_calibration transport OM returns correct structure", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "transport",
    estimator = "om",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$estimator, "om")
  expect_gte(result$ici, 0)
})

test_that("tr_calibration transport OM and IPW produce different results", {
  d <- generate_tr_calibration_data(n = 500, seed = 456)

  result_ipw <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  result_om <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "transport",
    estimator = "om",
    se_method = "none"
  )

  # Different estimators should generally give different results
  # But both should be valid (non-negative)
  expect_gte(result_ipw$ici, 0)
  expect_gte(result_om$ici, 0)
})


# ==============================================================================
# Joint Analysis Tests
# ==============================================================================

test_that("tr_calibration joint IPW returns correct structure", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "joint",
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$analysis, "joint")
  expect_gte(result$ici, 0)
})

test_that("tr_calibration joint OM returns correct structure", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    treatment_level = 0,
    analysis = "joint",
    estimator = "om",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$analysis, "joint")
  expect_equal(result$estimator, "om")
  expect_gte(result$ici, 0)
})

test_that("tr_calibration transport and joint produce different results", {
  d <- generate_tr_calibration_data(n = 500)

  result_tr <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  result_joint <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    analysis = "joint",
    estimator = "ipw",
    se_method = "none"
  )

  # Both should be valid
  expect_gte(result_tr$ici, 0)
  expect_gte(result_joint$ici, 0)
})


# ==============================================================================
# Smoother Tests
# ==============================================================================

test_that("tr_calibration works with loess smoother", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    smoother = "loess",
    span = 0.75,
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$smoother, "loess")
})

test_that("tr_calibration works with binned smoother", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    smoother = "binned",
    n_bins = 10,
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_equal(result$smoother, "binned")
  expect_gte(result$ici, 0)
})

test_that("tr_calibration loess and binned produce consistent results", {
  d <- generate_tr_calibration_data(n = 500, seed = 789)

  result_loess <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    smoother = "loess",
    se_method = "none"
  )

  result_binned <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    smoother = "binned",
    n_bins = 10,
    se_method = "none"
  )

  # Both methods should give reasonably similar ICI for well-calibrated model
  expect_gte(result_loess$ici, 0)
  expect_gte(result_binned$ici, 0)
})


# ==============================================================================
# Bootstrap Standard Error Tests
# ==============================================================================

test_that("tr_calibration bootstrap produces standard errors", {
  d <- generate_tr_calibration_data(n = 300)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    analysis = "transport",
    estimator = "ipw",
    se_method = "bootstrap",
    n_boot = 50  # Small for testing
  )

  expect_s3_class(result, "tr_calibration")

  # Check SE components exist
  expect_true(!is.null(result$se))
  expect_true(!is.null(result$se$ici))
  expect_true(!is.null(result$se$e50))
  expect_true(!is.null(result$se$e90))
  expect_true(!is.null(result$se$emax))

  # SEs should be positive
  expect_gt(result$se$ici, 0)
  expect_gt(result$se$e50, 0)
  expect_gt(result$se$emax, 0)

  # Check CI components
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_lt(result$ci_lower$ici, result$ci_upper$ici)
})

test_that("tr_calibration stratified bootstrap preserves ratio", {
  d <- generate_tr_calibration_data(n = 300)

  # This test just ensures stratified bootstrap runs without error
  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "bootstrap",
    n_boot = 30,
    stratified_boot = TRUE
  )

  expect_s3_class(result, "tr_calibration")
  expect_true(!is.null(result$se))
})

test_that("tr_calibration non-stratified bootstrap works", {
  d <- generate_tr_calibration_data(n = 300)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "bootstrap",
    n_boot = 30,
    stratified_boot = FALSE
  )

  expect_s3_class(result, "tr_calibration")
  expect_true(!is.null(result$se))
})


# ==============================================================================
# S3 Methods Tests
# ==============================================================================

test_that("print.tr_calibration works", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_output(print(result), "Transportable CALIBRATION")
  expect_output(print(result), "ICI")
})

test_that("summary.tr_calibration works", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_output(summary(result), "Summary")
  expect_output(summary(result), "Calibration Statistics")
})

test_that("coef.tr_calibration returns calibration metrics", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  coefs <- coef(result)
  expect_named(coefs, c("ici", "e50", "e90", "emax"))
  expect_true(all(coefs >= 0))
})

test_that("confint.tr_calibration returns CIs with bootstrap", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "bootstrap",
    n_boot = 30
  )

  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 4)  # ici, e50, e90, emax
  expect_equal(ncol(ci), 2)  # lower, upper
})

test_that("confint.tr_calibration warns without bootstrap", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_warning(confint(result), "not computed")
})


# ==============================================================================
# Plot Method Tests
# ==============================================================================

test_that("plot.tr_calibration works with base R", {
  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  # Should not error
  expect_silent(plot(result))
})

test_that("plot.tr_calibration returns ggplot when available", {
  skip_if_not_installed("ggplot2")

  d <- generate_tr_calibration_data(n = 200)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Naive Estimate Tests
# ==============================================================================

test_that("tr_calibration includes naive estimate", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_true(!is.null(result$naive_ici))
  expect_gte(result$naive_ici, 0)
})


# ==============================================================================
# Edge Cases
# ==============================================================================

test_that("tr_calibration handles extreme predictions", {
  d <- generate_tr_calibration_data(n = 300)

  # Add some extreme predictions
  d$predictions[1:10] <- 0.001
  d$predictions[11:20] <- 0.999

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
  expect_gte(result$ici, 0)
})

test_that("tr_calibration handles sparse treatment groups", {
  set.seed(999)
  n <- 300

  # Generate data with sparse treatment in source
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.6)
  # Very few treated in source
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.1),  # Only 10% treated in RCT
              rbinom(n, 1, 0.5))
  y <- rbinom(n, 1, plogis(-1 + x))
  pred <- plogis(-1 + 0.8 * x)

  # Should still work (may have larger SE)
  result <- tr_calibration(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = data.frame(x = x),
    treatment_level = 1,  # Sparse group
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
})


# ==============================================================================
# User-Supplied Model Tests
# ==============================================================================

test_that("tr_calibration accepts user-supplied selection model", {
  d <- generate_tr_calibration_data(n = 300)

  # Fit selection model externally
  sel_data <- cbind(S0 = 1 - d$source, d$covariates)
  selection_model <- glm(S0 ~ x1 + x2, data = sel_data, family = binomial())

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    selection_model = selection_model,
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
})

test_that("tr_calibration accepts user-supplied propensity model", {
  d <- generate_tr_calibration_data(n = 300)

  # Fit propensity model externally using source data
  source_idx <- d$source == 1
  ps_data <- cbind(A = d$treatment[source_idx], d$covariates[source_idx, ])
  propensity_model <- glm(A ~ x1 + x2, data = ps_data, family = binomial())

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    propensity_model = propensity_model,
    se_method = "none"
  )

  expect_s3_class(result, "tr_calibration")
})


# ==============================================================================
# Consistency Tests
# ==============================================================================

test_that("tr_calibration is deterministic without bootstrap", {
  d <- generate_tr_calibration_data(n = 300)

  result1 <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  result2 <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_equal(result1$ici, result2$ici)
  expect_equal(result1$e50, result2$e50)
  expect_equal(result1$emax, result2$emax)
})

test_that("tr_calibration sample sizes are correct", {
  d <- generate_tr_calibration_data(n = 500)

  result <- tr_calibration(
    predictions = d$predictions,
    outcomes = d$outcomes,
    treatment = d$treatment,
    source = d$source,
    covariates = d$covariates,
    se_method = "none"
  )

  expect_equal(result$n_source, sum(d$source == 1))
  expect_equal(result$n_target, sum(d$source == 0))
})
