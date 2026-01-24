# Benchmark Tests
#
# This file verifies that cfperformance estimators match established
# implementations when they should agree. These tests help validate
# correctness by comparing against known-good implementations.
#
# Key insight: The IPW estimator uses observation-level weights.
# When we have a simple random sample (no confounding), IPW weights = 1
# and our estimators should match standard unweighted versions.
# More importantly, when we manually compute IPW weights, our estimators
# should match weighted versions from other packages.

# =============================================================================
# BENCHMARK 1: cf_auc IPW vs WeightedROC
# =============================================================================
# The IPW estimator computes:
#   AUC_IPW = sum_{i,j} w_i w_j I(Y_i > Y_j) I(f_i > f_j) /
#             sum_{i,j} w_i w_j I(Y_i > Y_j)
# where w_i = I(A_i = a) / P(A = a | X_i)
#
# WeightedROC computes a weighted AUC using observation weights.
# These should match when we provide the same weights.

test_that("cf_auc IPW matches WeightedROC for weighted AUC", {
  skip_if_not_installed("WeightedROC")

  set.seed(42)
  n <- 200

  # Generate data
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  ps_true <- plogis(0.5 + 0.3 * X$x1 - 0.2 * X$x2)
  A <- rbinom(n, 1, ps_true)
  Y <- rbinom(n, 1, plogis(-0.5 + 0.5 * X$x1 + 0.3 * X$x2 - 0.5 * A))
  pred <- runif(n, 0, 1)  # Random predictions for testing

  # Fit propensity model
  ps_model <- glm(A ~ x1 + x2, data = cbind(X, A = A), family = binomial())
  ps_fitted <- predict(ps_model, type = "response")

  # Test for treatment_level = 1
  # IPW weights for A=1: I(A=1) / ps
  ipw_weights_1 <- ifelse(A == 1, 1 / ps_fitted, 0)

  # cf_auc with IPW
  cf_result_1 <- cf_auc(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 1,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"  # No trimming for exact comparison
  )

  # WeightedROC - filter to treated only and use weights
  idx_1 <- A == 1
  if (sum(Y[idx_1] == 1) > 0 && sum(Y[idx_1] == 0) > 0) {
    roc_1 <- WeightedROC::WeightedROC(
      guess = pred[idx_1],
      label = Y[idx_1],
      weight = ipw_weights_1[idx_1]
    )
    wroc_auc_1 <- WeightedROC::WeightedAUC(roc_1)

    # These should be very close (numerical precision)
    expect_equal(cf_result_1$estimate, wroc_auc_1, tolerance = 1e-6)
  }

  # Test for treatment_level = 0
  # IPW weights for A=0: I(A=0) / (1-ps)
  ipw_weights_0 <- ifelse(A == 0, 1 / (1 - ps_fitted), 0)

  cf_result_0 <- cf_auc(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 0,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"
  )

  idx_0 <- A == 0
  if (sum(Y[idx_0] == 1) > 0 && sum(Y[idx_0] == 0) > 0) {
    roc_0 <- WeightedROC::WeightedROC(
      guess = pred[idx_0],
      label = Y[idx_0],
      weight = ipw_weights_0[idx_0]
    )
    wroc_auc_0 <- WeightedROC::WeightedAUC(roc_0)

    expect_equal(cf_result_0$estimate, wroc_auc_0, tolerance = 1e-6)
  }
})


test_that("cf_auc naive matches pROC::auc for unweighted data", {
  skip_if_not_installed("pROC")

  set.seed(123)
  n <- 200

  Y <- rbinom(n, 1, 0.4)
  pred <- runif(n)

  # cf_auc naive
  cf_result <- cf_auc(
    predictions = pred,
    outcomes = Y,
    treatment = rep(1, n),  # Everyone treated
    covariates = data.frame(x = rnorm(n)),
    treatment_level = 1,
    estimator = "naive",
    se_method = "none"
  )

  # pROC::auc
  proc_result <- pROC::auc(Y, pred, direction = "<", quiet = TRUE)

  expect_equal(cf_result$estimate, as.numeric(proc_result), tolerance = 1e-6)
})


# =============================================================================
# BENCHMARK 2: cf_mse IPW vs weighted.mean
# =============================================================================
# The IPW estimator for MSE computes:
#   MSE_IPW = sum_i w_i (Y_i - f_i)^2 / sum_i w_i
# where w_i = I(A_i = a) / P(A = a | X_i)
#
# This is just a weighted mean of squared errors.

test_that("cf_mse IPW matches manual Horvitz-Thompson calculation", {
  set.seed(42)
  n <- 200

  # Generate data
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  ps_true <- plogis(0.5 + 0.3 * X$x1 - 0.2 * X$x2)
  A <- rbinom(n, 1, ps_true)
  Y <- rbinom(n, 1, plogis(-0.5 + 0.5 * X$x1 + 0.3 * X$x2 - 0.5 * A))
  pred <- runif(n)

  # Fit propensity model
  ps_model <- glm(A ~ x1 + x2, data = cbind(X, A = A), family = binomial())
  ps_fitted <- predict(ps_model, type = "response")

  # Test for treatment_level = 1
  # Note: cfperformance uses Horvitz-Thompson (HT) style IPW, not normalized
  # HT: mean(I(A=a)/ps * L) estimates E[L(Y(a), f)]
  I_a_1 <- as.numeric(A == 1)
  ipw_weights_1 <- I_a_1 / ps_fitted

  cf_result_1 <- cf_mse(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 1,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"
  )

  # Manual HT-style MSE
  loss_1 <- (Y - pred)^2
  manual_mse_1 <- mean(ipw_weights_1 * loss_1)

  expect_equal(cf_result_1$estimate, manual_mse_1, tolerance = 1e-10)

  # Test for treatment_level = 0
  I_a_0 <- as.numeric(A == 0)
  ipw_weights_0 <- I_a_0 / (1 - ps_fitted)

  cf_result_0 <- cf_mse(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 0,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"
  )

  manual_mse_0 <- mean(ipw_weights_0 * loss_1)

  expect_equal(cf_result_0$estimate, manual_mse_0, tolerance = 1e-10)
})


test_that("cf_mse naive matches mean squared error", {
  set.seed(123)
  n <- 200

  Y <- rbinom(n, 1, 0.4)
  pred <- runif(n)

  cf_result <- cf_mse(
    predictions = pred,
    outcomes = Y,
    treatment = rep(1, n),
    covariates = data.frame(x = rnorm(n)),
    treatment_level = 1,
    estimator = "naive",
    se_method = "none"
  )

  manual_mse <- mean((Y - pred)^2)

  expect_equal(cf_result$estimate, manual_mse, tolerance = 1e-10)
})


# =============================================================================
# BENCHMARK 3: cf_sensitivity/specificity IPW vs weighted proportions
# =============================================================================
# IPW sensitivity at threshold t:
#   Se_IPW(t) = sum_i w_i I(f_i >= t) I(Y_i = 1) / sum_i w_i I(Y_i = 1)
#
# This is a weighted proportion.

test_that("cf_sensitivity IPW matches weighted proportion", {
  set.seed(42)
  n <- 300

  X <- data.frame(x1 = rnorm(n))
  ps_true <- plogis(0.5 + 0.3 * X$x1)
  A <- rbinom(n, 1, ps_true)
  Y <- rbinom(n, 1, plogis(0.5 * X$x1 - 0.5 * A))
  pred <- plogis(0.3 * X$x1)  # Predictions between 0 and 1

  ps_model <- glm(A ~ x1, data = cbind(X, A = A), family = binomial())
  ps_fitted <- predict(ps_model, type = "response")

  threshold <- 0.5

  # Test treatment_level = 1
  ipw_weights_1 <- ifelse(A == 1, 1 / ps_fitted, 0)

  cf_sens_1 <- cf_sensitivity(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 1,
    threshold = threshold,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"
  )

  # Manual weighted sensitivity
  idx_1 <- A == 1
  cases_1 <- Y[idx_1] == 1
  predicted_pos_1 <- pred[idx_1] >= threshold

  # Weighted sensitivity = weighted sum of TP / weighted sum of cases
  manual_sens_1 <- sum(ipw_weights_1[idx_1] * cases_1 * predicted_pos_1) /
                   sum(ipw_weights_1[idx_1] * cases_1)

  expect_equal(cf_sens_1$estimate, manual_sens_1, tolerance = 1e-10)

  # Test specificity for treatment_level = 0
  ipw_weights_0 <- ifelse(A == 0, 1 / (1 - ps_fitted), 0)

  cf_spec_0 <- cf_specificity(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = X,
    treatment_level = 0,
    threshold = threshold,
    estimator = "ipw",
    propensity_model = ps_model,
    se_method = "none",
    ps_trim = "none"
  )

  # Manual weighted specificity
  idx_0 <- A == 0
  non_cases_0 <- Y[idx_0] == 0
  predicted_neg_0 <- pred[idx_0] < threshold

  manual_spec_0 <- sum(ipw_weights_0[idx_0] * non_cases_0 * predicted_neg_0) /
                   sum(ipw_weights_0[idx_0] * non_cases_0)

  expect_equal(cf_spec_0$estimate, manual_spec_0, tolerance = 1e-10)
})


# =============================================================================
# BENCHMARK 4: No confounding case - should match naive
# =============================================================================
# When treatment is randomized (no confounding), IPW should equal naive.

test_that("IPW equals naive when treatment is randomized", {
  set.seed(42)
  n <- 500

  # No confounding: treatment is random
  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, 0.5)  # Random treatment, no dependence on X
  Y <- rbinom(n, 1, plogis(0.5 * X$x1))  # Outcome depends on X, not A
  pred <- runif(n)

  # Fit propensity model (should estimate ~0.5 for everyone)
  ps_model <- glm(A ~ x1, data = cbind(X, A = A), family = binomial())

  # MSE comparison
  mse_naive <- cf_mse(pred, Y, A, X, 1, "naive", se_method = "none")
  mse_ipw <- cf_mse(pred, Y, A, X, 1, "ipw", propensity_model = ps_model,
                    se_method = "none", ps_trim = "none")

  # Should be close (not exact due to finite sample)
  expect_equal(mse_naive$estimate, mse_ipw$estimate, tolerance = 0.05)

  # AUC comparison
  auc_naive <- cf_auc(pred, Y, A, X, 1, "naive", se_method = "none")
  auc_ipw <- cf_auc(pred, Y, A, X, 1, "ipw", propensity_model = ps_model,
                    se_method = "none", ps_trim = "none")

  expect_equal(auc_naive$estimate, auc_ipw$estimate, tolerance = 0.05)
})


# =============================================================================
# BENCHMARK 5: Transportability - tr_mse IPW manual verification
# =============================================================================

test_that("tr_mse IPW matches manual formula from implementation", {
  set.seed(42)
  n <- 300

  # Generate combined source + target population
  x <- rnorm(n)

  # Source indicator - source has lower x values on average
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))  # Source = 1

  # Treatment - randomized in source, confounded in target
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))

  # Outcomes
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))

  # Predictions
  pred <- plogis(-1 + 0.8 * x)

  covariates <- data.frame(x = x)

  # Fit selection model: P(S=1|X) - note: selection model predicts P(S=0|X)
  # but we fit it as P(S=0|X) directly
  selection_model <- glm(I(1 - s) ~ x, data = data.frame(s = s, x = x), family = binomial())

  # Fit propensity model in source: P(A=1|X, S=1)
  source_idx <- s == 1
  propensity_model <- glm(a ~ x, data = data.frame(a = a, x = x)[source_idx, ],
                          family = binomial())

  treatment_level <- 0
  n0 <- sum(s == 0)  # Target sample size

  # Get model predictions for source observations
  p_s0_source <- predict(selection_model, newdata = data.frame(x = x[source_idx]),
                         type = "response")
  p_s1_source <- 1 - p_s0_source

  ps_treatment_source <- predict(propensity_model, type = "response")
  if (treatment_level == 0) {
    ps_treatment_source <- 1 - ps_treatment_source
  }

  # IPW formula from .tr_mse_transport:
  # psi = (1/n0) * sum_{S=1, A=a} [P(S=0|X) / (P(A=a|X,S=1) * P(S=1|X))] * L
  I_a_source <- a[source_idx] == treatment_level
  weights <- (p_s0_source / (ps_treatment_source * p_s1_source))[I_a_source]
  loss_source <- (y[source_idx] - pred[source_idx])^2

  # tr_mse
  tr_result <- tr_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covariates,
    treatment_level = treatment_level,
    analysis = "transport",
    estimator = "ipw",
    selection_model = selection_model,
    propensity_model = propensity_model,
    se_method = "none",
    ps_trim = "none"
  )

  # Manual calculation following exact formula
  manual_mse <- sum(weights * loss_source[I_a_source]) / n0

  # These should match exactly
  expect_equal(tr_result$estimate, manual_mse, tolerance = 1e-10)
})


# =============================================================================
# BENCHMARK 6: Consistency check - predictions at outcome should have AUC = 1
# =============================================================================

test_that("Perfect predictions yield AUC = 1", {
  set.seed(42)
  n <- 200

  X <- data.frame(x1 = rnorm(n))
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, 0.5)

  # Perfect predictions (with small noise to break ties)
  pred <- Y + rnorm(n, 0, 0.001)

  ps_model <- glm(A ~ x1, data = cbind(X, A = A), family = binomial())

  # All estimators should give AUC ≈ 1
  auc_naive <- cf_auc(pred, Y, A, X, 1, "naive", se_method = "none")
  auc_ipw <- cf_auc(pred, Y, A, X, 1, "ipw", propensity_model = ps_model,
                    se_method = "none", ps_trim = "none")

  expect_equal(auc_naive$estimate, 1.0, tolerance = 0.01)
  expect_equal(auc_ipw$estimate, 1.0, tolerance = 0.01)
})


test_that("Random predictions yield AUC ≈ 0.5",
  {
  set.seed(42)
  n <- 1000  # Large n for stable estimate

  Y <- rbinom(n, 1, 0.5)
  pred <- runif(n)  # Random predictions, independent of Y

  result <- cf_auc(
    predictions = pred,
    outcomes = Y,
    treatment = rep(1, n),
    covariates = data.frame(x = rnorm(n)),
    treatment_level = 1,
    estimator = "naive",
    se_method = "none"
  )

  # Should be close to 0.5
  expect_equal(result$estimate, 0.5, tolerance = 0.05)
})


# =============================================================================
# BENCHMARK 7: Known formula verification
# =============================================================================
# For a specific simple case, verify against hand calculation

test_that("Manual AUC calculation matches for simple case", {
  # Simple case: 4 observations
  Y <- c(1, 1, 0, 0)
  pred <- c(0.9, 0.7, 0.4, 0.2)

  # AUC = proportion of concordant pairs

# Pairs (Yi=1, Yj=0):
  # (0.9, 0.4): concordant
  # (0.9, 0.2): concordant
  # (0.7, 0.4): concordant
  # (0.7, 0.2): concordant
  # All 4 pairs concordant, AUC = 4/4 = 1

  result <- cf_auc(
    predictions = pred,
    outcomes = Y,
    treatment = rep(1, 4),
    covariates = data.frame(x = 1:4),
    treatment_level = 1,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result$estimate, 1.0)

  # Now with some discordance
  pred2 <- c(0.9, 0.3, 0.6, 0.2)  # Second case has lower prediction than third
  # Pairs:
  # (0.9, 0.6): concordant
  # (0.9, 0.2): concordant
  # (0.3, 0.6): discordant (case has lower pred than non-case)
  # (0.3, 0.2): concordant
  # 3 concordant, 1 discordant -> AUC = 3/4 = 0.75

  result2 <- cf_auc(
    predictions = pred2,
    outcomes = Y,
    treatment = rep(1, 4),
    covariates = data.frame(x = 1:4),
    treatment_level = 1,
    estimator = "naive",
    se_method = "none"
  )

  expect_equal(result2$estimate, 0.75)
})


# =============================================================================
# Additional benchmark ideas for future expansion:
# =============================================================================
#
# 1. Compare calibration curves to rms::val.prob or CalibrationCurves package
#
# 2. Compare bootstrap SEs to theoretical SEs for simple cases where
#    analytical formulas exist (e.g., binomial proportions)
#
# 3. For DR estimator: verify double robustness property by showing
#    consistency when one (but not both) nuisance model is misspecified
#    (this is partly done in test-simulation-study.R)
#
# 4. Compare cross-fitting variance to honest variance bounds
#
# 5. Verify that influence function SEs match bootstrap SEs asymptotically
#    (would need large n simulation)
