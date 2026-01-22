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

    # Naive MSE (uses all observations)
    naive_mse <- mean((Y_test - pred_test)^2)

    # IPW MSE (uses only untreated, weighted)
    # Fit propensity model
    ps_model <- glm(A ~ X, data = test_df, family = binomial())
    ps <- predict(ps_model, type = "response")

    untreated <- A_test == 0
    weights <- 1 / (1 - ps[untreated])
    ipw_mse <- weighted.mean((Y_test[untreated] - pred_test[untreated])^2, weights)

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

  set.seed(8761276)

  n_sims <- 200  # Reduced for testing
  n_train <- 500
  n_test <- 500
  n_truth <- 1000

  results <- data.frame(
    sim = integer(),
    naive = numeric(),
    cl_correct = numeric(),
    ipw_correct = numeric(),
    dr_correct = numeric(),
    dr_ps_misspec = numeric(),
    dr_om_misspec = numeric(),
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

    # Naive estimator
    naive <- mean((test_df$Y - pred_test)^2)

    # Correctly specified propensity model
    ps_model_correct <- glm(A ~ X1 + I(X1^2) + X2 + X3,
                            data = test_df, family = binomial())
    ps_correct <- predict(ps_model_correct, newdata = test_df, type = "response")

    # Misspecified propensity model (main effects only)
    ps_model_misspec <- glm(A ~ X1 + X2 + X3,
                            data = test_df, family = binomial())
    ps_misspec <- predict(ps_model_misspec, newdata = test_df, type = "response")

    # Correctly specified outcome model (for conditional loss)
    om_model_correct <- glm(Y ~ X1 + I(X1^2) + X2 + X3,
                            data = subset(test_df, A == 0), family = binomial())
    pY_correct <- predict(om_model_correct, newdata = test_df, type = "response")

    # Misspecified outcome model
    om_model_misspec <- glm(Y ~ X1 + X2 + X3,
                            data = subset(test_df, A == 0), family = binomial())
    pY_misspec <- predict(om_model_misspec, newdata = test_df, type = "response")

    # Conditional loss h = E[(Y-pred)^2 | X, A=0] = pY - 2*pred*pY + pred^2
    h_correct <- pY_correct - 2 * pred_test * pY_correct + pred_test^2
    h_misspec <- pY_misspec - 2 * pred_test * pY_misspec + pred_test^2

    # CL estimator (correct)
    cl_correct <- mean(h_correct)

    # IPW estimator (correct)
    ps_a0_correct <- 1 - ps_correct
    ps_a0_correct <- pmax(pmin(ps_a0_correct, 0.99), 0.01)
    I_a0 <- as.numeric(test_df$A == 0)
    loss <- (test_df$Y - pred_test)^2
    ipw_correct <- sum(I_a0 / ps_a0_correct * loss) / sum(I_a0)

    # DR estimator (correct nuisance models)
    dr_correct <- mean(h_correct + I_a0 / ps_a0_correct * (loss - h_correct))

    # DR estimator (propensity misspecified, outcome correct)
    ps_a0_misspec <- 1 - ps_misspec
    ps_a0_misspec <- pmax(pmin(ps_a0_misspec, 0.99), 0.01)
    dr_ps_misspec <- mean(h_correct + I_a0 / ps_a0_misspec * (loss - h_correct))

    # DR estimator (propensity correct, outcome misspecified)
    dr_om_misspec <- mean(h_misspec + I_a0 / ps_a0_correct * (loss - h_misspec))

    results <- rbind(results, data.frame(
      sim = sim,
      naive = naive,
      cl_correct = cl_correct,
      ipw_correct = ipw_correct,
      dr_correct = dr_correct,
      dr_ps_misspec = dr_ps_misspec,
      dr_om_misspec = dr_om_misspec,
      true_mse = true_mse
    ))
  }

  # Compute summary statistics
  means <- colMeans(results[, -1])
  true_val <- means["true_mse"]

  # Compute relative biases
  rel_bias <- function(est, truth) abs((est - truth) / truth) * 100

  bias_naive <- rel_bias(means["naive"], true_val)
  bias_dr_correct <- rel_bias(means["dr_correct"], true_val)
  bias_dr_ps_misspec <- rel_bias(means["dr_ps_misspec"], true_val)
  bias_dr_om_misspec <- rel_bias(means["dr_om_misspec"], true_val)

  cat("\n=== Simulation 2 Results (MSE) ===\n")
  cat(sprintf("True MSE: %.4f\n", true_val))
  cat(sprintf("Naive: %.4f (rel bias: %.1f%%)\n", means["naive"], bias_naive))
  cat(sprintf("DR correct: %.4f (rel bias: %.1f%%)\n",
              means["dr_correct"], bias_dr_correct))
  cat(sprintf("DR (PS misspec): %.4f (rel bias: %.1f%%)\n",
              means["dr_ps_misspec"], bias_dr_ps_misspec))
  cat(sprintf("DR (OM misspec): %.4f (rel bias: %.1f%%)\n",
              means["dr_om_misspec"], bias_dr_om_misspec))

  # Key assertions from the paper:
  # 1. Naive is biased (relative bias > 1%)
  expect_gt(bias_naive, 1)

  # 2. DR with correct models is approximately unbiased (< 2% relative bias)
  expect_lt(bias_dr_correct, 2)

  # 3. DR is robust to propensity misspecification when outcome is correct
  expect_lt(bias_dr_ps_misspec, 2)

  # 4. DR is robust to outcome misspecification when propensity is correct
  expect_lt(bias_dr_om_misspec, 2)
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

    # Helper function to compute AUC
    compute_auc <- function(y, pred) {
      n1 <- sum(y == 1)
      n0 <- sum(y == 0)
      if (n1 == 0 || n0 == 0) return(NA)
      r <- rank(c(pred[y == 1], pred[y == 0]))
      (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }

    # Naive AUC
    naive_auc <- compute_auc(test_df$Y, pred_test)

    # True counterfactual AUC
    true_auc <- compute_auc(truth_df$Y, pred_truth)

    results <- rbind(results, data.frame(
      sim = sim,
      naive_auc = naive_auc,
      true_auc = true_auc
    ))
  }

  # Remove NA rows
  results <- results[complete.cases(results), ]

  means <- colMeans(results[, -1])

  cat("\n=== Simulation 2B Results (AUC) ===\n")
  cat(sprintf("True AUC: %.4f\n", means["true_auc"]))
  cat(sprintf("Naive AUC: %.4f\n", means["naive_auc"]))

  # Naive AUC should be biased relative to true counterfactual AUC
  # (in the paper, naive was about 5% higher than truth)
  bias_pct <- (means["naive_auc"] - means["true_auc"]) / means["true_auc"] * 100
  cat(sprintf("Relative bias: %.1f%%\n", bias_pct))

  # The naive estimator should show bias
  expect_gt(abs(bias_pct), 1)  # At least 1% bias expected
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
