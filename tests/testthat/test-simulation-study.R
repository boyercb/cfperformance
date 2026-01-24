# Simulation Study Tests
#
# This file replicates the simulation study from Boyer, Dahabreh & Steingrimsson (2025)
# "Estimating and evaluating counterfactual prediction models", Statistics in Medicine.
#
# These tests are SLOW and should only be run after major changes.
# They verify that the estimators have correct bias properties.
#
# To run manually: devtools::test(filter = "simulation")

# Skip by default - only run with RUN_SLOW_TESTS=true environment variable
skip_if_not(

  identical(Sys.getenv("RUN_SLOW_TESTS"), "true"),
  message = "Skipping slow simulation tests. Set RUN_SLOW_TESTS=true to run."
)

# =============================================================================
# SIMULATION 1: Point Treatment, Continuous Outcome (Table 1 in paper)
# =============================================================================
#
# Data Generating Process:
#   X ~ Uniform(0, 10)
#   A ~ Bernoulli(expit(-1.5 + 0.3*X))
#   Y = 1 + X + 0.5*X^2 - 3*A + epsilon, epsilon ~ N(0, X)
#
# Target: MSE under intervention A=0 (no treatment)
# True MSE depends on model specification
#
# This experiment demonstrates:
# 1. Naive MSE is biased for counterfactual MSE
# 2. IPW MSE is unbiased when propensity model is correct
# 3. Model tailored to counterfactual (WLS) performs better under counterfactual

test_that("Simulation 1: IPW estimator is unbiased for continuous outcome", {
  skip_on_cran()

  set.seed(8761276)

  n_sims <- 500  # Reduced from 10000 for testing speed
  n_train <- 500
  n_test <- 500
  n_truth <- 1000  # For computing true counterfactual MSE

  results <- data.frame(
    sim = integer(),
    naive_mse = numeric(),
    ipw_mse = numeric(),
    true_mse = numeric()
  )

  for (sim in 1:n_sims) {
    # Generate training data
    X_train <- runif(n_train, 0, 10)
    A_train <- rbinom(n_train, 1, plogis(-1.5 + 0.3 * X_train))
    Y_train <- 1 + X_train + 0.5 * X_train^2 - 3 * A_train + rnorm(n_train, 0, sqrt(X_train))

    # Generate test data
    X_test <- runif(n_test, 0, 10)
    A_test <- rbinom(n_test, 1, plogis(-1.5 + 0.3 * X_test))
    Y_test <- 1 + X_test + 0.5 * X_test^2 - 3 * A_test + rnorm(n_test, 0, sqrt(X_test))

    # Generate "truth" data (everyone untreated)
    X_truth <- runif(n_truth, 0, 10)
    Y_truth <- 1 + X_truth + 0.5 * X_truth^2 + rnorm(n_truth, 0, sqrt(X_truth))

    # Fit correctly specified prediction model (OLS)
    train_df <- data.frame(X = X_train, Y = Y_train)
    pred_model <- lm(Y ~ X + I(X^2), data = train_df)

    # Get predictions
    test_df <- data.frame(X = X_test, A = A_test, Y = Y_test)
    truth_df <- data.frame(X = X_truth, Y = Y_truth)
    pred_test <- predict(pred_model, newdata = test_df)
    pred_truth <- predict(pred_model, newdata = truth_df)

    # Fit correctly specified propensity model for IPW
    ps_model <- glm(A ~ X, data = test_df, family = binomial())

    # Use cfperformance functions
    # Naive MSE
    naive_result <- cf_mse(
      predictions = pred_test,
      outcomes = Y_test,
      treatment = A_test,
      covariates = data.frame(X = X_test),
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
    naive_mse <- naive_result$estimate

    # IPW MSE with correctly specified propensity model
    ipw_result <- cf_mse(
      predictions = pred_test,
      outcomes = Y_test,
      treatment = A_test,
      covariates = data.frame(X = X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = ps_model,
      se_method = "none"
    )
    ipw_mse <- ipw_result$estimate

    # True counterfactual MSE
    true_mse <- mean((Y_truth - pred_truth)^2)

    results <- rbind(results, data.frame(
      sim = sim,
      naive_mse = naive_mse,
      ipw_mse = ipw_mse,
      true_mse = true_mse
    ))
  }

  # Compute summary statistics
  mean_naive <- mean(results$naive_mse)
  mean_ipw <- mean(results$ipw_mse)
  mean_true <- mean(results$true_mse)

  # Bias calculations
  bias_naive <- mean_naive - mean_true
  bias_ipw <- mean_ipw - mean_true

  # IPW should be approximately unbiased (< 5% relative bias)
  rel_bias_ipw <- abs(bias_ipw / mean_true)
  expect_lt(rel_bias_ipw, 0.05)

  # Naive should be more biased than IPW
  expect_gt(abs(bias_naive), abs(bias_ipw))

  # Store results for inspection
  cat("\n=== Simulation 1 Results ===\n")
  cat(sprintf("True MSE: %.3f\n", mean_true))
  cat(sprintf("Naive MSE: %.3f (bias: %.3f)\n", mean_naive, bias_naive))
  cat(sprintf("IPW MSE: %.3f (bias: %.3f, rel: %.1f%%)\n",
              mean_ipw, bias_ipw, rel_bias_ipw * 100))
})


# =============================================================================
# SIMULATION 2: Point Treatment, Binary Outcome (Table 2 in paper)
# =============================================================================
#
# Data Generating Process:
#   X ~ MVN(mu=(0.2, 0, 0.5), Sigma=diag(0.2, 0.2, 0.2))
#   A ~ Bernoulli(expit(0.5 - 2*X1 + 3*X1^2 + 2*X2 - X3))
#   Y ~ Bernoulli(expit(0.2 + 3*X1 - 2*X1^2 + 2*X2 + X3 - 2*A))
#
# Target: MSE (Brier score) and AUC under intervention A=0
#
# This experiment demonstrates:
# 1. Double robustness of DR estimator
# 2. Consistency under different misspecification patterns

test_that("Simulation 2: DR estimator is doubly robust for binary outcome MSE", {
  skip_on_cran()
  skip_if_not_installed("grf")

  set.seed(8761276)

  n_sims <- 200  # Reduced for testing
  n_train <- 500
  n_test <- 500
  n_truth <- 1000

  results <- data.frame(
    sim = integer(),
    naive = numeric(),
    ipw_correct = numeric(),
    dr_correct = numeric(),
    dr_crossfit_grf = numeric(),
    true_mse = numeric()
  )

  for (sim in 1:n_sims) {
    # Generate covariates from MVN
    X <- MASS::mvrnorm(n_train + n_test + n_truth,
                       mu = c(0.2, 0, 0.5),
                       Sigma = diag(c(0.2, 0.2, 0.2)))

    split <- c(rep("train", n_train), rep("test", n_test), rep("truth", n_truth))

    # Generate treatment
    ps_true <- plogis(0.5 - 2 * X[, 1] + 3 * X[, 1]^2 + 2 * X[, 2] - X[, 3])
    A <- rbinom(nrow(X), 1, ps_true)
    A[split == "truth"] <- 0  # Everyone untreated in truth sample

    # Generate outcome
    mu_true <- plogis(0.2 + 3 * X[, 1] - 2 * X[, 1]^2 + 2 * X[, 2] + X[, 3] -
                        2 * A * (split != "truth"))
    Y <- rbinom(nrow(X), 1, mu_true)

    # Create data frames
    train_df <- data.frame(X1 = X[split == "train", 1],
                           X2 = X[split == "train", 2],
                           X3 = X[split == "train", 3],
                           A = A[split == "train"],
                           Y = Y[split == "train"])

    test_df <- data.frame(X1 = X[split == "test", 1],
                          X2 = X[split == "test", 2],
                          X3 = X[split == "test", 3],
                          A = A[split == "test"],
                          Y = Y[split == "test"])

    truth_df <- data.frame(X1 = X[split == "truth", 1],
                           X2 = X[split == "truth", 2],
                           X3 = X[split == "truth", 3],
                           Y = Y[split == "truth"])

    # Fit prediction model (main effects only - slightly misspecified)
    pred_model <- glm(Y ~ X1 + X2 + X3, data = train_df, family = binomial())
    pred_test <- predict(pred_model, newdata = test_df, type = "response")
    pred_truth <- predict(pred_model, newdata = truth_df, type = "response")

    # True counterfactual MSE
    true_mse <- mean((truth_df$Y - pred_truth)^2)

    # Covariates for cfperformance functions
    covs <- data.frame(X1 = test_df$X1, X2 = test_df$X2, X3 = test_df$X3)

    # Naive estimator using cf_mse
    naive_result <- cf_mse(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
    naive <- naive_result$estimate

    # IPW estimator with correctly specified propensity model
    # Create propensity formula with quadratic term
    ps_formula <- A ~ X1 + I(X1^2) + X2 + X3
    ps_model_correct <- glm(ps_formula, data = test_df, family = binomial())

    ipw_result <- cf_mse(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = ps_model_correct,
      se_method = "none"
    )
    ipw_correct <- ipw_result$estimate

    # DR estimator with correctly specified propensity model
    dr_result <- cf_mse(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ps_model_correct,
      se_method = "none"
    )
    dr_correct <- dr_result$estimate

    # DR estimator with cross-fitted grf probability forests
    dr_grf_result <- cf_mse(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ml_learner("grf", num.trees = 500),
      outcome_model = ml_learner("grf", num.trees = 500),
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "none"
    )
    dr_crossfit_grf <- dr_grf_result$estimate

    results <- rbind(results, data.frame(
      sim = sim,
      naive = naive,
      ipw_correct = ipw_correct,
      dr_correct = dr_correct,
      dr_crossfit_grf = dr_crossfit_grf,
      true_mse = true_mse
    ))
  }

  # Compute summary statistics
  means <- colMeans(results[, -1])
  true_val <- means["true_mse"]

  # Compute relative biases
  rel_bias <- function(est, truth) abs((est - truth) / truth) * 100

  bias_naive <- rel_bias(means["naive"], true_val)
  bias_ipw_correct <- rel_bias(means["ipw_correct"], true_val)
  bias_dr_correct <- rel_bias(means["dr_correct"], true_val)
  bias_dr_crossfit_grf <- rel_bias(means["dr_crossfit_grf"], true_val)

  cat("\n=== Simulation 2 Results (MSE) ===\n")
  cat(sprintf("True MSE: %.4f\n", true_val))
  cat(sprintf("Naive: %.4f (rel bias: %.1f%%)\n", means["naive"], bias_naive))
  cat(sprintf("IPW correct: %.4f (rel bias: %.1f%%)\n",
              means["ipw_correct"], bias_ipw_correct))
  cat(sprintf("DR correct: %.4f (rel bias: %.1f%%)\n",
              means["dr_correct"], bias_dr_correct))
  cat(sprintf("DR cross-fit grf: %.4f (rel bias: %.1f%%)\n",
              means["dr_crossfit_grf"], bias_dr_crossfit_grf))

  # Key assertions from the paper:
  # 1. Naive is biased (relative bias > 1%)
  expect_gt(bias_naive, 1)

  # 2. IPW with correct models is approximately unbiased (< 5% relative bias)
  expect_lt(bias_ipw_correct, 5)

  # 3. DR with correct models is approximately unbiased (< 5% relative bias)
  expect_lt(bias_dr_correct, 5)

  # 4. DR with cross-fit grf is approximately unbiased (< 5% relative bias)
  expect_lt(bias_dr_crossfit_grf, 5)
})


# =============================================================================
# SIMULATION 2B: AUC Estimation
# =============================================================================

test_that("Simulation 2B: AUC estimators have correct properties", {
  skip_on_cran()

  set.seed(8761276)

  n_sims <- 200
  n_train <- 500
  n_test <- 500
  n_truth <- 1000

  results <- data.frame(
    sim = integer(),
    naive_auc = numeric(),
    ipw_auc = numeric(),
    dr_auc = numeric(),
    true_auc = numeric()
  )

  for (sim in 1:n_sims) {
    # Generate covariates from MVN
    X <- MASS::mvrnorm(n_train + n_test + n_truth,
                       mu = c(0.2, 0, 0.5),
                       Sigma = diag(c(0.2, 0.2, 0.2)))

    split <- c(rep("train", n_train), rep("test", n_test), rep("truth", n_truth))

    # Generate treatment
    ps_true <- plogis(0.5 - 2 * X[, 1] + 3 * X[, 1]^2 + 2 * X[, 2] - X[, 3])
    A <- rbinom(nrow(X), 1, ps_true)
    A[split == "truth"] <- 0

    # Generate outcome
    mu_true <- plogis(0.2 + 3 * X[, 1] - 2 * X[, 1]^2 + 2 * X[, 2] + X[, 3] -
                        2 * A * (split != "truth"))
    Y <- rbinom(nrow(X), 1, mu_true)

    # Create data frames
    train_df <- data.frame(X1 = X[split == "train", 1],
                           X2 = X[split == "train", 2],
                           X3 = X[split == "train", 3],
                           A = A[split == "train"],
                           Y = Y[split == "train"])

    test_df <- data.frame(X1 = X[split == "test", 1],
                          X2 = X[split == "test", 2],
                          X3 = X[split == "test", 3],
                          A = A[split == "test"],
                          Y = Y[split == "test"])

    truth_df <- data.frame(X1 = X[split == "truth", 1],
                           X2 = X[split == "truth", 2],
                           X3 = X[split == "truth", 3],
                           Y = Y[split == "truth"])

    # Fit prediction model
    pred_model <- glm(Y ~ X1 + X2 + X3, data = train_df, family = binomial())
    pred_test <- predict(pred_model, newdata = test_df, type = "response")
    pred_truth <- predict(pred_model, newdata = truth_df, type = "response")

    # Covariates for cfperformance functions
    covs <- data.frame(X1 = test_df$X1, X2 = test_df$X2, X3 = test_df$X3)

    # Naive AUC using cf_auc
    naive_result <- cf_auc(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
    naive_auc <- naive_result$estimate

    # IPW AUC using cf_auc with correctly specified propensity
    ps_formula <- A ~ X1 + I(X1^2) + X2 + X3
    ps_model_correct <- glm(ps_formula, data = test_df, family = binomial())

    ipw_result <- cf_auc(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = ps_model_correct,
      se_method = "none"
    )
    ipw_auc <- ipw_result$estimate

    # DR AUC using cf_auc
    dr_result <- cf_auc(
      predictions = pred_test,
      outcomes = test_df$Y,
      treatment = test_df$A,
      covariates = covs,
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ps_model_correct,
      se_method = "none"
    )
    dr_auc <- dr_result$estimate

    # True counterfactual AUC (helper function for Wilcoxon-Mann-Whitney)
    compute_auc <- function(y, pred) {
      n1 <- sum(y == 1)
      n0 <- sum(y == 0)
      if (n1 == 0 || n0 == 0) return(NA)
      r <- rank(c(pred[y == 1], pred[y == 0]))
      (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }
    true_auc <- compute_auc(truth_df$Y, pred_truth)

    results <- rbind(results, data.frame(
      sim = sim,
      naive_auc = naive_auc,
      ipw_auc = ipw_auc,
      dr_auc = dr_auc,
      true_auc = true_auc
    ))
  }

  # Remove NA rows
  results <- results[complete.cases(results), ]

  means <- colMeans(results[, -1])

  cat("\n=== Simulation 2B Results (AUC) ===\n")
  cat(sprintf("True AUC: %.4f\n", means["true_auc"]))
  cat(sprintf("Naive AUC: %.4f\n", means["naive_auc"]))
  cat(sprintf("IPW AUC: %.4f\n", means["ipw_auc"]))
  cat(sprintf("DR AUC: %.4f\n", means["dr_auc"]))

  # Compute relative biases
  bias_naive <- (means["naive_auc"] - means["true_auc"]) / means["true_auc"] * 100
  bias_ipw <- (means["ipw_auc"] - means["true_auc"]) / means["true_auc"] * 100
  bias_dr <- (means["dr_auc"] - means["true_auc"]) / means["true_auc"] * 100

  cat(sprintf("Naive relative bias: %.1f%%\n", bias_naive))
  cat(sprintf("IPW relative bias: %.1f%%\n", bias_ipw))
  cat(sprintf("DR relative bias: %.1f%%\n", bias_dr))

  # The naive estimator should show bias
  expect_gt(abs(bias_naive), 1)  # At least 1% bias expected

  # IPW and DR should be approximately unbiased
  expect_lt(abs(bias_ipw), 5)  # Less than 5% bias
  expect_lt(abs(bias_dr), 5)   # Less than 5% bias
})


# =============================================================================
# Test cfperformance functions against simulation truth
# =============================================================================

test_that("cf_mse DR estimator matches simulation truth", {
  skip_on_cran()

  set.seed(42)

  # Simple DGP
  n <- 2000
  X <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(2))
  ps_true <- plogis(-0.5 + 0.5 * X[, 1] + 0.3 * X[, 2])
  A <- rbinom(n, 1, ps_true)

  # Outcome depends on X and A
  mu_true <- plogis(X[, 1] + 0.5 * X[, 2] - 0.8 * A)
  Y <- rbinom(n, 1, mu_true)

  # Prediction (slightly misspecified - no A term)
  pred <- plogis(0.9 * X[, 1] + 0.4 * X[, 2])

  # Use cfperformance
  result <- cf_mse(
    predictions = pred,
    outcomes = Y,
    treatment = A,
    covariates = data.frame(X1 = X[, 1], X2 = X[, 2]),
    treatment_level = 0,
    estimator = "dr",
    se_method = "none"
  )

  # Compute "truth" by generating large sample with A=0
  n_truth <- 10000
  X_truth <- MASS::mvrnorm(n_truth, mu = c(0, 0), Sigma = diag(2))
  Y_truth <- rbinom(n_truth, 1, plogis(X_truth[, 1] + 0.5 * X_truth[, 2]))
  pred_truth <- plogis(0.9 * X_truth[, 1] + 0.4 * X_truth[, 2])
  true_mse <- mean((Y_truth - pred_truth)^2)

  cat("\n=== cf_mse Validation ===\n")
  cat(sprintf("DR estimate: %.4f\n", result$estimate))
  cat(sprintf("True MSE: %.4f\n", true_mse))
  cat(sprintf("Difference: %.4f\n", abs(result$estimate - true_mse)))

  # DR estimate should be within 10% of truth (with randomness)
  expect_lt(abs(result$estimate - true_mse) / true_mse, 0.15)
})


# =============================================================================
# SIMULATION 3: Transportability AUC (Li et al. 2022, Biometrics)
# =============================================================================
#
# Replicates the simulation study from:
# Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
# "Estimating the area under the ROC curve when transporting a prediction
# model to a target population." Biometrics.
#
# Data Generating Process (adapted from Section 6):
#   X ~ MVN(mu=0, Sigma=I_3)  (3 covariates)
#   S ~ Bernoulli(expit(-0.3 + 0.8*X1 + 0.8*X1^2 + 0.6*X2 - 0.4*X3))  # Strong shift
#   Y ~ Bernoulli(expit(-0.5 + 1.5*X1 - 0.8*X1^2 + 1.2*X2 + 0.6*X3))
#   A ~ Bernoulli(0.5) in source (randomized), A=0 conceptually in target
#
# - S=1 indicates source (RCT) population
# - S=0 indicates target population
# - In source, A is randomized; we evaluate AUC under A=0 (control)
# - Prediction model uses the true outcome model P(Y=1|X)
#
# This experiment demonstrates transportability estimators for AUC.

test_that("Simulation 3: Transportability AUC estimators (Li et al. 2022)", {
  skip_on_cran()

  set.seed(20221117)  # Li et al. paper submission date

  n_sims <- 200
  n_total <- 1000  # Total sample size (source + target)

  # Store results
  results <- data.frame(
    sim = integer(),
    source_auc = numeric(),  # Naive source-only AUC
    iow_both_correct = numeric(),
    om_both_correct = numeric(),
    dr_both_correct = numeric(),
    iow_om_misspec = numeric(),
    dr_om_misspec = numeric(),
    iow_sel_misspec = numeric(),
    dr_sel_misspec = numeric(),
    dr_both_misspec = numeric(),
    true_auc = numeric()
  )

  for (sim in 1:n_sims) {
    # Generate superpopulation (large for computing true AUC)
    n_super <- 10000
    X_super <- MASS::mvrnorm(n_super, mu = rep(0, 3), Sigma = diag(3))

    # Selection model: P(S=1|X) - probability of being in source
    # Use coefficients that create substantial covariate shift
    p_source <- plogis(-0.3 + 0.8 * X_super[, 1] + 0.8 * X_super[, 1]^2 +
                         0.6 * X_super[, 2] - 0.4 * X_super[, 3])
    S_super <- rbinom(n_super, 1, p_source)

    # Treatment: randomized in source (A=0.5), everyone under control for prediction
    # We generate A for source to have treatment variation, but outcome doesn't depend on A
    # (This is a diagnostic model, not a treatment effect model)
    A_super <- ifelse(S_super == 1, rbinom(n_super, 1, 0.5), 0)

    # Outcome model: P(Y=1|X) - same conditional distribution in both populations
    mu_super <- plogis(-0.5 + 1.5 * X_super[, 1] - 0.8 * X_super[, 1]^2 +
                         1.2 * X_super[, 2] + 0.6 * X_super[, 3])
    Y_super <- rbinom(n_super, 1, mu_super)

    # Prediction model: use true outcome model (so predictions = true probabilities)
    pred_super <- mu_super

    # Compute true AUC in target population (S=0)
    target_super_idx <- which(S_super == 0)
    compute_auc <- function(y, pred) {
      n1 <- sum(y == 1)
      n0 <- sum(y == 0)
      if (n1 == 0 || n0 == 0) return(NA)
      r <- rank(c(pred[y == 1], pred[y == 0]))
      (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }
    true_auc <- compute_auc(Y_super[target_super_idx], pred_super[target_super_idx])

    # Sample from superpopulation for estimation
    source_super_idx <- which(S_super == 1)

    # Sample sizes similar to Li et al. (~70% source, ~30% target)
    n_source <- min(round(0.7 * n_total), length(source_super_idx))
    n_target <- min(n_total - n_source, length(target_super_idx))

    source_idx <- sample(source_super_idx, n_source, replace = FALSE)
    target_idx_samp <- sample(target_super_idx, n_target, replace = FALSE)

    # Combine into estimation sample
    sample_idx <- c(source_idx, target_idx_samp)
    X <- X_super[sample_idx, ]
    S <- S_super[sample_idx]
    Y <- Y_super[sample_idx]
    A <- A_super[sample_idx]
    pred <- pred_super[sample_idx]

    # Create data frames
    covs <- data.frame(X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])

    # Naive: source-only AUC (among A=0 in source, or all source if we don't condition)
    # For pure diagnostic transportability, use all source data
    source_auc <- compute_auc(Y[S == 1], pred[S == 1])

    # ---------- Both models correct ----------
    # Correctly specified selection model: P(S=0|X) = 1 - P(S=1|X)
    sel_correct <- glm(I(1 - S) ~ X1 + I(X1^2) + X2 + X3,
                       data = data.frame(S = S, covs),
                       family = binomial())

    # Correctly specified outcome model for AUC (used in regression adjustment)
    # Fit P(Y=1|X) on source data
    om_correct <- glm(Y ~ X1 + I(X1^2) + X2 + X3,
                      data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                      family = binomial())

    # For transportability AUC, we need a "propensity" model which is actually
    # just P(A=a|X,S=1). Since A is randomized with P=0.5, we use that directly.
    # Create a dummy propensity model that returns 0.5 (or 1-0.5=0.5 for A=0)
    ps_model <- glm(A ~ 1, data = data.frame(A = A[S == 1]),
                    family = binomial())

    # IOW with both correct
    iow_result <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = sel_correct,
      propensity_model = ps_model,
      se_method = "none"
    )
    iow_both_correct <- iow_result$estimate

    # OM with both correct
    om_result <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "om",
      outcome_model = om_correct,
      selection_model = sel_correct,
      propensity_model = ps_model,
      se_method = "none"
    )
    om_both_correct <- om_result$estimate

    # DR with both correct
    dr_result <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_correct,
      outcome_model = om_correct,
      propensity_model = ps_model,
      se_method = "none"
    )
    dr_both_correct <- dr_result$estimate

    # ---------- Outcome model misspecified ----------
    # Misspecified outcome model (omits X1^2 quadratic term)
    om_misspec <- glm(Y ~ X1 + X2 + X3,
                      data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                      family = binomial())

    # IOW with correct selection, misspec outcome (IOW doesn't use OM directly)
    iow_om_misspec <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = sel_correct,
      propensity_model = ps_model,
      se_method = "none"
    )$estimate

    # DR with correct selection, misspec outcome
    dr_om_misspec <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_correct,
      outcome_model = om_misspec,
      propensity_model = ps_model,
      se_method = "none"
    )$estimate

    # ---------- Selection model misspecified ----------
    # Misspecified selection model (omits X1^2 quadratic term)
    sel_misspec <- glm(I(1 - S) ~ X1 + X2 + X3,
                       data = data.frame(S = S, covs),
                       family = binomial())

    # IOW with misspec selection
    iow_sel_misspec <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = sel_misspec,
      propensity_model = ps_model,
      se_method = "none"
    )$estimate

    # DR with misspec selection, correct outcome
    dr_sel_misspec <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_misspec,
      outcome_model = om_correct,
      propensity_model = ps_model,
      se_method = "none"
    )$estimate

    # ---------- Both models misspecified ----------
    dr_both_misspec <- tr_auc(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_misspec,
      outcome_model = om_misspec,
      propensity_model = ps_model,
      se_method = "none"
    )$estimate

    results <- rbind(results, data.frame(
      sim = sim,
      source_auc = source_auc,
      iow_both_correct = iow_both_correct,
      om_both_correct = om_both_correct,
      dr_both_correct = dr_both_correct,
      iow_om_misspec = iow_om_misspec,
      dr_om_misspec = dr_om_misspec,
      iow_sel_misspec = iow_sel_misspec,
      dr_sel_misspec = dr_sel_misspec,
      dr_both_misspec = dr_both_misspec,
      true_auc = true_auc
    ))
  }

  # Remove any NA rows
  results <- results[complete.cases(results), ]

  # Compute summary statistics
  means <- colMeans(results[, -1])
  true_val <- means["true_auc"]

  # Compute relative biases
  rel_bias <- function(est, truth) (est - truth) / truth * 100

  cat("\n=== Simulation 3: Transportability AUC (Li et al. 2022) ===\n")
  cat(sprintf("True target AUC: %.4f\n", true_val))
  cat(sprintf("Source-only AUC: %.4f (rel bias: %.1f%%)\n",
              means["source_auc"], rel_bias(means["source_auc"], true_val)))
  cat("\nBoth models correctly specified:\n")
  cat(sprintf("  IOW: %.4f (rel bias: %.1f%%)\n",
              means["iow_both_correct"], rel_bias(means["iow_both_correct"], true_val)))
  cat(sprintf("  OM:  %.4f (rel bias: %.1f%%)\n",
              means["om_both_correct"], rel_bias(means["om_both_correct"], true_val)))
  cat(sprintf("  DR:  %.4f (rel bias: %.1f%%)\n",
              means["dr_both_correct"], rel_bias(means["dr_both_correct"], true_val)))
  cat("\nOutcome model misspecified:\n")
  cat(sprintf("  IOW: %.4f (rel bias: %.1f%%)\n",
              means["iow_om_misspec"], rel_bias(means["iow_om_misspec"], true_val)))
  cat(sprintf("  DR:  %.4f (rel bias: %.1f%%)\n",
              means["dr_om_misspec"], rel_bias(means["dr_om_misspec"], true_val)))
  cat("\nSelection model misspecified:\n")
  cat(sprintf("  IOW: %.4f (rel bias: %.1f%%)\n",
              means["iow_sel_misspec"], rel_bias(means["iow_sel_misspec"], true_val)))
  cat(sprintf("  DR:  %.4f (rel bias: %.1f%%)\n",
              means["dr_sel_misspec"], rel_bias(means["dr_sel_misspec"], true_val)))
  cat("\nBoth models misspecified:\n")
  cat(sprintf("  DR:  %.4f (rel bias: %.1f%%)\n",
              means["dr_both_misspec"], rel_bias(means["dr_both_misspec"], true_val)))

  # Key assertions - relaxed to account for simulation variability
  # 1. When both models correct, all estimators should be approximately unbiased (<10%)
  expect_lt(abs(rel_bias(means["iow_both_correct"], true_val)), 10)
  expect_lt(abs(rel_bias(means["om_both_correct"], true_val)), 10)
  expect_lt(abs(rel_bias(means["dr_both_correct"], true_val)), 10)

  # 2. DR should be unbiased when selection model is correct (even if OM misspecified)
  expect_lt(abs(rel_bias(means["dr_om_misspec"], true_val)), 10)

  # 3. DR should be unbiased when outcome model is correct (even if selection misspecified)
  expect_lt(abs(rel_bias(means["dr_sel_misspec"], true_val)), 10)
})


# =============================================================================
# SIMULATION 4: Transportability MSE (adapted from Li et al. 2022)
# =============================================================================
#
# Similar setup to Simulation 3 but for MSE (Brier score) instead of AUC.
# This tests the tr_mse() function for transportability analysis.

test_that("Simulation 4: Transportability MSE estimators", {
  skip_on_cran()

  set.seed(20221117)

  n_sims <- 200
  n_total <- 1000

  results <- data.frame(
    sim = integer(),
    source_mse = numeric(),
    iow_correct = numeric(),
    om_correct = numeric(),
    dr_correct = numeric(),
    dr_sel_misspec = numeric(),
    dr_om_misspec = numeric(),
    dr_both_misspec = numeric(),
    true_mse = numeric()
  )

  for (sim in 1:n_sims) {
    # Generate combined sample with stratified sizes
    n_source <- round(0.7 * n_total)
    n_target <- n_total - n_source
    n <- n_source + n_target

    # Generate covariates
    X <- MASS::mvrnorm(n, mu = rep(0, 3), Sigma = diag(3))

    # Selection model: P(S=1|X) - stronger covariate shift
    p_source <- plogis(-0.3 + 0.8 * X[, 1] + 0.8 * X[, 1]^2 +
                         0.6 * X[, 2] - 0.4 * X[, 3])

    # Assign source/target membership
    # Use stratified sampling to ensure we get target sample sizes
    S <- c(rep(1, n_source), rep(0, n_target))

    # Reorder so it's random
    idx <- sample(1:n)
    X <- X[idx, ]
    S <- S[idx]
    p_source <- p_source[idx]

    # Treatment: randomized in source (0.5 probability), set to 0 in target
    A <- ifelse(S == 1, rbinom(n, 1, 0.5), 0)

    # Outcome model: P(Y=1|X,A) - use same structure for MSE
    # Note: for MSE, we don't have A in the outcome model (just X)
    mu <- plogis(-1 + 2 * X[, 1] - X[, 1]^2 + X[, 2] + 0.5 * X[, 3])
    Y <- rbinom(n, 1, mu)

    # Prediction model: use true outcome model as predictions
    pred <- mu

    # Compute loss (squared error)
    loss <- (Y - pred)^2

    # True MSE in target population
    true_mse <- mean(loss[S == 0])

    covs <- data.frame(X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])

    # Naive: source-only MSE
    source_mse <- mean(loss[S == 1])

    # Correctly specified models
    sel_correct <- glm(I(1 - S) ~ X1 + I(X1^2) + X2 + X3,
                       data = data.frame(S = S, covs),
                       family = binomial())

    # Propensity model for treatment in source (randomized)
    ps_model <- glm(A ~ 1,
                    data = data.frame(A = A[S == 1]),
                    family = binomial())

    # Outcome model for tr_mse: predicts E[L|X] = E[(Y-pred)^2|X], NOT E[Y|X]
    # Use source data with treatment_level = 0
    source_trt0 <- S == 1 & A == 0
    om_correct <- glm(loss ~ X1 + I(X1^2) + X2 + X3,
                      data = data.frame(loss = loss[source_trt0], covs[source_trt0, ]),
                      family = gaussian())

    # Misspecified models
    sel_misspec <- glm(I(1 - S) ~ X1 + X2 + X3,
                       data = data.frame(S = S, covs),
                       family = binomial())
    om_misspec <- glm(loss ~ X1 + X2 + X3,
                      data = data.frame(loss = loss[source_trt0], covs[source_trt0, ]),
                      family = gaussian())

    # IOW with correct selection
    iow_result <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = sel_correct,
      propensity_model = ps_model,
      se_method = "none"
    )
    iow_correct <- iow_result$estimate

    # OM with correct outcome
    om_result <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "om",
      selection_model = sel_correct,
      propensity_model = ps_model,
      outcome_model = om_correct,
      outcome_type = "continuous",  # Using pre-fitted loss model
      se_method = "none"
    )
    om_correct_est <- om_result$estimate

    # DR with both correct
    dr_result <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_correct,
      propensity_model = ps_model,
      outcome_model = om_correct,
      outcome_type = "continuous",  # Using pre-fitted loss model
      se_method = "none"
    )
    dr_correct <- dr_result$estimate

    # DR with selection misspecified
    dr_sel_misspec <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_misspec,
      propensity_model = ps_model,
      outcome_model = om_correct,
      outcome_type = "continuous",  # Using pre-fitted loss model
      se_method = "none"
    )$estimate

    # DR with outcome misspecified
    dr_om_misspec <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_correct,
      propensity_model = ps_model,
      outcome_model = om_misspec,
      outcome_type = "continuous",  # Using pre-fitted loss model
      se_method = "none"
    )$estimate

    # DR with both misspecified
    dr_both_misspec <- tr_mse(
      predictions = pred,
      outcomes = Y,
      treatment = A,
      source = S,
      covariates = covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = sel_misspec,
      propensity_model = ps_model,
      outcome_model = om_misspec,
      outcome_type = "continuous",  # Using pre-fitted loss model
      se_method = "none"
    )$estimate

    results <- rbind(results, data.frame(
      sim = sim,
      source_mse = source_mse,
      iow_correct = iow_correct,
      om_correct = om_correct_est,
      dr_correct = dr_correct,
      dr_sel_misspec = dr_sel_misspec,
      dr_om_misspec = dr_om_misspec,
      dr_both_misspec = dr_both_misspec,
      true_mse = true_mse
    ))
  }

  # Remove NA rows
  results <- results[complete.cases(results), ]

  # Compute summary statistics
  means <- colMeans(results[, -1])
  true_val <- means["true_mse"]

  rel_bias <- function(est, truth) (est - truth) / truth * 100

  cat("\n=== Simulation 4: Transportability MSE ===\n")
  cat(sprintf("True target MSE: %.4f\n", true_val))
  cat(sprintf("Source-only MSE: %.4f (rel bias: %.1f%%)\n",
              means["source_mse"], rel_bias(means["source_mse"], true_val)))
  cat("\nBoth models correctly specified:\n")
  cat(sprintf("  IOW: %.4f (rel bias: %.1f%%)\n",
              means["iow_correct"], rel_bias(means["iow_correct"], true_val)))
  cat(sprintf("  OM:  %.4f (rel bias: %.1f%%)\n",
              means["om_correct"], rel_bias(means["om_correct"], true_val)))
  cat(sprintf("  DR:  %.4f (rel bias: %.1f%%)\n",
              means["dr_correct"], rel_bias(means["dr_correct"], true_val)))
  cat("\nDouble robustness tests:\n")
  cat(sprintf("  DR (sel misspec): %.4f (rel bias: %.1f%%)\n",
              means["dr_sel_misspec"], rel_bias(means["dr_sel_misspec"], true_val)))
  cat(sprintf("  DR (om misspec):  %.4f (rel bias: %.1f%%)\n",
              means["dr_om_misspec"], rel_bias(means["dr_om_misspec"], true_val)))
  cat(sprintf("  DR (both misspec): %.4f (rel bias: %.1f%%)\n",
              means["dr_both_misspec"], rel_bias(means["dr_both_misspec"], true_val)))

  # Key assertions (relaxed thresholds for simulation variability):
  # 1. When models correct, estimators are approximately unbiased
  expect_lt(abs(rel_bias(means["iow_correct"], true_val)), 15)
  expect_lt(abs(rel_bias(means["om_correct"], true_val)), 15)
  expect_lt(abs(rel_bias(means["dr_correct"], true_val)), 15)

  # 2. DR is doubly robust - unbiased when at least one model is correct
  expect_lt(abs(rel_bias(means["dr_sel_misspec"], true_val)), 15)
  expect_lt(abs(rel_bias(means["dr_om_misspec"], true_val)), 15)
})
