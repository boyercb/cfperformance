#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive Simulation Study for cfperformance Package
# =============================================================================
#
# This script replicates and extends simulation studies from:
#
# 1. Boyer, Dahabreh & Steingrimsson (2025). "Estimating and evaluating
#    counterfactual prediction models." Statistics in Medicine.
#    doi:10.1002/sim.70287
#
# 2. Li, Gatsonis, Dahabreh & Steingrimsson (2022). "Estimating the area
#    under the ROC curve when transporting a prediction model to a target
#    population." Biometrics. doi:10.1111/biom.13796
#
# Run with: Rscript inst/simulation/run_simulation_study.R
# Or source in R: source("inst/simulation/run_simulation_study.R")
#
# Results are saved to inst/simulation/results/
# =============================================================================

# Load required packages
library(cfperformance)
library(MASS)
library(parallel)

# Configuration
CONFIG <- list(
  # Simulation parameters
  n_sims = 1000,          # Number of simulations
  n_train = 1000,         # Training sample size
  n_test = 1000,          # Test sample size
  n_truth = 10000,        # Large sample for computing true values

  # Bootstrap parameters for SE estimation
  n_boot = 200,           # Bootstrap replications
  conf_level = 0.95,      # Confidence level


  # Parallel processing
  use_parallel = TRUE,

  n_cores = max(1, detectCores() - 2),

  # Output
  output_dir = "inst/simulation/results",
  seed = 20250124
)

# Create output directory if needed
if (!dir.exists(CONFIG$output_dir)) {
  dir.create(CONFIG$output_dir, recursive = TRUE)
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Compute simulation metrics
#' @param estimates Vector of point estimates
#' @param ses Vector of standard errors (optional)
#' @param ci_lowers Vector of CI lower bounds (optional)
#' @param ci_uppers Vector of CI upper bounds (optional)
#' @param true_value True parameter value
#' @param conf_level Confidence level for coverage
compute_metrics <- function(estimates, ses = NULL, ci_lowers = NULL,
                            ci_uppers = NULL, true_value, conf_level = 0.95) {

  # Remove NAs
  valid <- !is.na(estimates)
  if (!is.null(ses)) valid <- valid & !is.na(ses)
  if (!is.null(ci_lowers)) valid <- valid & !is.na(ci_lowers)
  if (!is.null(ci_uppers)) valid <- valid & !is.na(ci_uppers)

  n_valid <- sum(valid)
  estimates <- estimates[valid]

  metrics <- list(
    n_sims = n_valid,
    true_value = true_value,
    mean_estimate = mean(estimates),
    bias = mean(estimates) - true_value,
    rel_bias_pct = (mean(estimates) - true_value) / abs(true_value) * 100,
    empirical_se = sd(estimates),
    rmse = sqrt(mean((estimates - true_value)^2))
  )

  if (!is.null(ses)) {
    ses <- ses[valid]
    metrics$mean_se = mean(ses)
    metrics$se_ratio = mean(ses) / sd(estimates)  # Should be ~1
  }

  if (!is.null(ci_lowers) && !is.null(ci_uppers)) {
    ci_lowers <- ci_lowers[valid]
    ci_uppers <- ci_uppers[valid]
    metrics$coverage = mean(ci_lowers <= true_value & ci_uppers >= true_value)
    metrics$mean_ci_width = mean(ci_uppers - ci_lowers)
  }

  return(metrics)
}


#' Format metrics as a data frame row
format_metrics_row <- function(metrics, estimator_name) {
  data.frame(
    Estimator = estimator_name,
    N_sims = metrics$n_sims,
    True = round(metrics$true_value, 4),
    Mean = round(metrics$mean_estimate, 4),
    Bias = round(metrics$bias, 4),
    Rel_Bias_Pct = round(metrics$rel_bias_pct, 2),
    Emp_SE = round(metrics$empirical_se, 4),
    Mean_SE = if (!is.null(metrics$mean_se)) round(metrics$mean_se, 4) else NA,
    SE_Ratio = if (!is.null(metrics$se_ratio)) round(metrics$se_ratio, 2) else NA,
    RMSE = round(metrics$rmse, 4),
    Coverage = if (!is.null(metrics$coverage)) round(metrics$coverage, 3) else NA,
    CI_Width = if (!is.null(metrics$mean_ci_width)) round(metrics$mean_ci_width, 4) else NA,
    stringsAsFactors = FALSE
  )
}


#' Compute AUC using Wilcoxon-Mann-Whitney
compute_auc <- function(y, pred) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  r <- rank(c(pred[y == 1], pred[y == 0]))
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


#' Run a single simulation replicate
#' @param sim Simulation number
#' @param dgp Data generating process function
#' @param estimators List of estimator functions to run
#' @param ... Additional arguments passed to DGP
run_single_sim <- function(sim, dgp, estimators, ...) {
  # Generate data
  data <- dgp(...)

  # Run each estimator
  results <- lapply(names(estimators), function(est_name) {
    tryCatch({
      result <- estimators[[est_name]](data)
      data.frame(
        sim = sim,
        estimator = est_name,
        estimate = result$estimate,
        se = if (!is.null(result$se)) result$se else NA,
        ci_lower = if (!is.null(result$ci_lower)) result$ci_lower else NA,
        ci_upper = if (!is.null(result$ci_upper)) result$ci_upper else NA,
        true_value = data$true_value,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        sim = sim,
        estimator = est_name,
        estimate = NA,
        se = NA,
        ci_lower = NA,
        ci_upper = NA,
        true_value = data$true_value,
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, results)
}


# =============================================================================
# SIMULATION 1: Counterfactual MSE - Continuous Outcome
# (Boyer, Dahabreh & Steingrimsson 2025, Table 1)
# =============================================================================
#
# DGP:
#   X ~ Uniform(0, 10)
#   A ~ Bernoulli(expit(-1.5 + 0.3*X))
#   Y = 1 + X + 0.5*X^2 - 3*A + epsilon, epsilon ~ N(0, sqrt(X))
#
# Target: MSE under intervention A=0

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 1: Counterfactual MSE - Continuous Outcome\n")
cat("(Boyer, Dahabreh & Steingrimsson 2025)\n")
cat(strrep("=", 70), "\n\n")

dgp_continuous <- function(n_train, n_test, n_truth) {
  # Training data
  X_train <- runif(n_train, 0, 10)
  A_train <- rbinom(n_train, 1, plogis(-1.5 + 0.3 * X_train))
  Y_train <- 1 + X_train + 0.5 * X_train^2 - 3 * A_train +
             rnorm(n_train, 0, sqrt(X_train))

  # Test data
  X_test <- runif(n_test, 0, 10)
  A_test <- rbinom(n_test, 1, plogis(-1.5 + 0.3 * X_test))
  Y_test <- 1 + X_test + 0.5 * X_test^2 - 3 * A_test +
            rnorm(n_test, 0, sqrt(X_test))

  # Truth data (everyone untreated)
  X_truth <- runif(n_truth, 0, 10)
  Y_truth <- 1 + X_truth + 0.5 * X_truth^2 + rnorm(n_truth, 0, sqrt(X_truth))

  # Fit prediction model on training data
  train_df <- data.frame(X = X_train, Y = Y_train)
  pred_model <- lm(Y ~ X + I(X^2), data = train_df)

  # Get predictions
  pred_test <- predict(pred_model, newdata = data.frame(X = X_test))
  pred_truth <- predict(pred_model, newdata = data.frame(X = X_truth))

  # Fit propensity model on test data
  ps_model <- glm(A_test ~ X_test, family = binomial())

  # True MSE
  true_value <- mean((Y_truth - pred_truth)^2)

  list(
    X_test = X_test,
    A_test = A_test,
    Y_test = Y_test,
    pred_test = pred_test,
    ps_model = ps_model,
    true_value = true_value
  )
}

estimators_sim1 <- list(
  Naive = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
  },

  IPW = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      conf_level = CONFIG$conf_level
    )
  },

  DR = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      conf_level = CONFIG$conf_level
    )
  },

  IPW_Influence = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 1
set.seed(CONFIG$seed)
cat("Running", CONFIG$n_sims, "simulations...\n")

if (CONFIG$use_parallel && CONFIG$n_cores > 1) {
  cl <- makeCluster(CONFIG$n_cores)
  clusterExport(cl, c("dgp_continuous", "estimators_sim1", "CONFIG",
                       "run_single_sim"))
  clusterEvalQ(cl, library(cfperformance))

  results_sim1_list <- parLapply(cl, 1:CONFIG$n_sims, function(sim) {
    run_single_sim(sim, dgp_continuous, estimators_sim1,
                   n_train = CONFIG$n_train,
                   n_test = CONFIG$n_test,
                   n_truth = CONFIG$n_truth)
  })
  stopCluster(cl)
  results_sim1 <- do.call(rbind, results_sim1_list)
} else {
  results_sim1_list <- lapply(1:CONFIG$n_sims, function(sim) {
    if (sim %% 100 == 0) cat("  Simulation", sim, "of", CONFIG$n_sims, "\n")
    run_single_sim(sim, dgp_continuous, estimators_sim1,
                   n_train = CONFIG$n_train,
                   n_test = CONFIG$n_test,
                   n_truth = CONFIG$n_truth)
  })
  results_sim1 <- do.call(rbind, results_sim1_list)
}

# Compute metrics for each estimator
true_val_sim1 <- mean(results_sim1$true_value, na.rm = TRUE)
table_sim1 <- do.call(rbind, lapply(unique(results_sim1$estimator), function(est) {
  idx <- results_sim1$estimator == est
  m <- compute_metrics(
    estimates = results_sim1$estimate[idx],
    ses = results_sim1$se[idx],
    ci_lowers = results_sim1$ci_lower[idx],
    ci_uppers = results_sim1$ci_upper[idx],
    true_value = true_val_sim1,
    conf_level = CONFIG$conf_level
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 1 Results:\n")
print(table_sim1, row.names = FALSE)


# =============================================================================
# SIMULATION 2: Counterfactual MSE & AUC - Binary Outcome
# (Boyer, Dahabreh & Steingrimsson 2025, Table 2)
# =============================================================================
#
# DGP:
#   X ~ MVN(mu=(0.2, 0, 0.5), Sigma=diag(0.2, 0.2, 0.2))
#   A ~ Bernoulli(expit(0.5 - 2*X1 + 3*X1^2 + 2*X2 - X3))
#   Y ~ Bernoulli(expit(0.2 + 3*X1 - 2*X1^2 + 2*X2 + X3 - 2*A))
#
# Target: MSE and AUC under intervention A=0

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 2: Counterfactual MSE & AUC - Binary Outcome\n")
cat("(Boyer, Dahabreh & Steingrimsson 2025)\n")
cat(strrep("=", 70), "\n\n")

dgp_binary <- function(n_train, n_test, n_truth) {
  n_total <- n_train + n_test + n_truth

  # Generate covariates
  X <- mvrnorm(n_total, mu = c(0.2, 0, 0.5), Sigma = diag(c(0.2, 0.2, 0.2)))
  split <- c(rep("train", n_train), rep("test", n_test), rep("truth", n_truth))

  # Treatment model
  ps_true <- plogis(0.5 - 2 * X[, 1] + 3 * X[, 1]^2 + 2 * X[, 2] - X[, 3])
  A <- rbinom(n_total, 1, ps_true)
  A[split == "truth"] <- 0

  # Outcome model (depends on A in train/test, not in truth)
  mu_true <- plogis(0.2 + 3 * X[, 1] - 2 * X[, 1]^2 + 2 * X[, 2] + X[, 3] -
                      2 * A * (split != "truth"))
  Y <- rbinom(n_total, 1, mu_true)

  # Split data
  train_idx <- split == "train"
  test_idx <- split == "test"
  truth_idx <- split == "truth"

  # Fit prediction model on training data (main effects only - slightly misspecified)
  train_df <- data.frame(X1 = X[train_idx, 1], X2 = X[train_idx, 2],
                         X3 = X[train_idx, 3], Y = Y[train_idx])
  pred_model <- glm(Y ~ X1 + X2 + X3, data = train_df, family = binomial())

  # Predictions
  test_df <- data.frame(X1 = X[test_idx, 1], X2 = X[test_idx, 2],
                        X3 = X[test_idx, 3])
  truth_df <- data.frame(X1 = X[truth_idx, 1], X2 = X[truth_idx, 2],
                         X3 = X[truth_idx, 3])
  pred_test <- predict(pred_model, newdata = test_df, type = "response")
  pred_truth <- predict(pred_model, newdata = truth_df, type = "response")

  # Fit propensity model with correct specification on test data
  ps_df <- data.frame(A = A[test_idx], test_df)
  ps_model_correct <- glm(A ~ X1 + I(X1^2) + X2 + X3, data = ps_df, family = binomial())
  ps_model_misspec <- glm(A ~ X1 + X2 + X3, data = ps_df, family = binomial())

  # True values
  true_mse <- mean((Y[truth_idx] - pred_truth)^2)
  true_auc <- compute_auc(Y[truth_idx], pred_truth)

  list(
    X_test = X[test_idx, ],
    A_test = A[test_idx],
    Y_test = Y[test_idx],
    pred_test = pred_test,
    ps_model_correct = ps_model_correct,
    ps_model_misspec = ps_model_misspec,
    true_mse = true_mse,
    true_auc = true_auc,
    true_value = true_mse  # For MSE simulations
  )
}

estimators_sim2_mse <- list(
  Naive = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
  },

  IPW_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  },

  IPW_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_misspec,
      se_method = "none"
    )
  },

  DR_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  }
)

# Run simulation 2 - MSE
set.seed(CONFIG$seed + 1)
cat("Running", CONFIG$n_sims, "MSE simulations...\n")

results_sim2_mse_list <- lapply(1:CONFIG$n_sims, function(sim) {
  if (sim %% 100 == 0) cat("  Simulation", sim, "of", CONFIG$n_sims, "\n")
  run_single_sim(sim, dgp_binary, estimators_sim2_mse,
                 n_train = CONFIG$n_train,
                 n_test = CONFIG$n_test,
                 n_truth = CONFIG$n_truth)
})
results_sim2_mse <- do.call(rbind, results_sim2_mse_list)

# Compute metrics
true_val_sim2_mse <- mean(results_sim2_mse$true_value, na.rm = TRUE)
table_sim2_mse <- do.call(rbind, lapply(unique(results_sim2_mse$estimator), function(est) {
  idx <- results_sim2_mse$estimator == est
  m <- compute_metrics(
    estimates = results_sim2_mse$estimate[idx],
    ses = results_sim2_mse$se[idx],
    ci_lowers = results_sim2_mse$ci_lower[idx],
    ci_uppers = results_sim2_mse$ci_upper[idx],
    true_value = true_val_sim2_mse
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 2 MSE Results:\n")
print(table_sim2_mse, row.names = FALSE)


# =============================================================================
# SIMULATION 2B: Counterfactual AUC - Binary Outcome
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 2B: Counterfactual AUC - Binary Outcome\n")
cat(strrep("=", 70), "\n\n")

# Modify DGP to return true_auc as true_value
dgp_binary_auc <- function(n_train, n_test, n_truth) {
  data <- dgp_binary(n_train, n_test, n_truth)
  data$true_value <- data$true_auc
  data
}

estimators_sim2_auc <- list(
  Naive = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "naive",
      se_method = "none"
    )
  },

  IPW_Correct = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  },

  DR_Correct = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  }
)

# Run simulation 2B - AUC
set.seed(CONFIG$seed + 2)
cat("Running", CONFIG$n_sims, "AUC simulations...\n")

results_sim2_auc_list <- lapply(1:CONFIG$n_sims, function(sim) {
  if (sim %% 100 == 0) cat("  Simulation", sim, "of", CONFIG$n_sims, "\n")
  run_single_sim(sim, dgp_binary_auc, estimators_sim2_auc,
                 n_train = CONFIG$n_train,
                 n_test = CONFIG$n_test,
                 n_truth = CONFIG$n_truth)
})
results_sim2_auc <- do.call(rbind, results_sim2_auc_list)

# Compute metrics
true_val_sim2_auc <- mean(results_sim2_auc$true_value, na.rm = TRUE)
table_sim2_auc <- do.call(rbind, lapply(unique(results_sim2_auc$estimator), function(est) {
  idx <- results_sim2_auc$estimator == est
  m <- compute_metrics(
    estimates = results_sim2_auc$estimate[idx],
    ses = results_sim2_auc$se[idx],
    ci_lowers = results_sim2_auc$ci_lower[idx],
    ci_uppers = results_sim2_auc$ci_upper[idx],
    true_value = true_val_sim2_auc
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 2B AUC Results:\n")
print(table_sim2_auc, row.names = FALSE)


# =============================================================================
# SIMULATION 3: Transportability AUC (Li et al. 2022)
# =============================================================================
#
# DGP:
#   X ~ MVN(0, I_3)
#   S ~ Bernoulli(expit(-0.3 + 0.8*X1 + 0.8*X1^2 + 0.6*X2 - 0.4*X3))
#   A ~ Bernoulli(0.5) in source (randomized)
#   Y ~ Bernoulli(expit(-0.5 + 1.5*X1 - 0.8*X1^2 + 1.2*X2 + 0.6*X3))
#
# Target: AUC in target population (S=0)

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 3: Transportability AUC\n")
cat("(Li, Gatsonis, Dahabreh & Steingrimsson 2022)\n")
cat(strrep("=", 70), "\n\n")

dgp_transport_auc <- function(n_total, n_truth = 10000) {
  # Generate superpopulation for true AUC
  X_super <- mvrnorm(n_truth, mu = rep(0, 3), Sigma = diag(3))
  p_source_super <- plogis(-0.3 + 0.8 * X_super[, 1] + 0.8 * X_super[, 1]^2 +
                            0.6 * X_super[, 2] - 0.4 * X_super[, 3])
  S_super <- rbinom(n_truth, 1, p_source_super)
  mu_super <- plogis(-0.5 + 1.5 * X_super[, 1] - 0.8 * X_super[, 1]^2 +
                      1.2 * X_super[, 2] + 0.6 * X_super[, 3])
  Y_super <- rbinom(n_truth, 1, mu_super)

  # True AUC in target (S=0)
  target_idx <- S_super == 0
  true_auc <- compute_auc(Y_super[target_idx], mu_super[target_idx])

  # Generate estimation sample
  X <- mvrnorm(n_total, mu = rep(0, 3), Sigma = diag(3))
  p_source <- plogis(-0.3 + 0.8 * X[, 1] + 0.8 * X[, 1]^2 +
                      0.6 * X[, 2] - 0.4 * X[, 3])
  S <- rbinom(n_total, 1, p_source)

  # Treatment: randomized in source
  A <- ifelse(S == 1, rbinom(n_total, 1, 0.5), 0)

  # Outcome
  mu <- plogis(-0.5 + 1.5 * X[, 1] - 0.8 * X[, 1]^2 +
                1.2 * X[, 2] + 0.6 * X[, 3])
  Y <- rbinom(n_total, 1, mu)

  # Predictions (use true model)
  pred <- mu

  covs <- data.frame(X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])

  # Correctly specified models
  sel_correct <- glm(I(1 - S) ~ X1 + I(X1^2) + X2 + X3,
                     data = data.frame(S = S, covs), family = binomial())
  sel_misspec <- glm(I(1 - S) ~ X1 + X2 + X3,
                     data = data.frame(S = S, covs), family = binomial())

  om_correct <- glm(Y ~ X1 + I(X1^2) + X2 + X3,
                    data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                    family = binomial())
  om_misspec <- glm(Y ~ X1 + X2 + X3,
                    data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                    family = binomial())

  ps_model <- glm(A ~ 1, data = data.frame(A = A[S == 1]), family = binomial())

  list(
    X = X,
    S = S,
    A = A,
    Y = Y,
    pred = pred,
    covs = covs,
    sel_correct = sel_correct,
    sel_misspec = sel_misspec,
    om_correct = om_correct,
    om_misspec = om_misspec,
    ps_model = ps_model,
    true_value = true_auc,
    source_auc = compute_auc(Y[S == 1], pred[S == 1])
  )
}

estimators_sim3 <- list(
  Source_Only = function(data) {
    list(estimate = data$source_auc, se = NULL, ci_lower = NULL, ci_upper = NULL)
  },

  IOW_Both_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = data$sel_correct,
      propensity_model = data$ps_model,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  },

  DR_Both_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_correct,
      outcome_model = data$om_correct,
      propensity_model = data$ps_model,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot
    )
  },

  DR_Sel_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_misspec,
      outcome_model = data$om_correct,
      propensity_model = data$ps_model,
      se_method = "none"
    )
  },

  DR_OM_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_correct,
      outcome_model = data$om_misspec,
      propensity_model = data$ps_model,
      se_method = "none"
    )
  },

  DR_Both_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_misspec,
      outcome_model = data$om_misspec,
      propensity_model = data$ps_model,
      se_method = "none"
    )
  }
)

# Run simulation 3
set.seed(CONFIG$seed + 3)
cat("Running", CONFIG$n_sims, "transportability AUC simulations...\n")

results_sim3_list <- lapply(1:CONFIG$n_sims, function(sim) {
  if (sim %% 100 == 0) cat("  Simulation", sim, "of", CONFIG$n_sims, "\n")
  run_single_sim(sim, dgp_transport_auc, estimators_sim3,
                 n_total = CONFIG$n_test * 2)
})
results_sim3 <- do.call(rbind, results_sim3_list)

# Compute metrics
true_val_sim3 <- mean(results_sim3$true_value, na.rm = TRUE)
table_sim3 <- do.call(rbind, lapply(unique(results_sim3$estimator), function(est) {
  idx <- results_sim3$estimator == est
  m <- compute_metrics(
    estimates = results_sim3$estimate[idx],
    ses = results_sim3$se[idx],
    ci_lowers = results_sim3$ci_lower[idx],
    ci_uppers = results_sim3$ci_upper[idx],
    true_value = true_val_sim3
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 3 Transportability AUC Results:\n")
print(table_sim3, row.names = FALSE)


# =============================================================================
# Save Results
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Saving Results\n")
cat(strrep("=", 70), "\n\n")

# Save individual result tables
write.csv(table_sim1, file.path(CONFIG$output_dir, "sim1_continuous_mse.csv"),
          row.names = FALSE)
write.csv(table_sim2_mse, file.path(CONFIG$output_dir, "sim2_binary_mse.csv"),
          row.names = FALSE)
write.csv(table_sim2_auc, file.path(CONFIG$output_dir, "sim2_binary_auc.csv"),
          row.names = FALSE)
write.csv(table_sim3, file.path(CONFIG$output_dir, "sim3_transport_auc.csv"),
          row.names = FALSE)

# Save raw results for further analysis
saveRDS(list(
  config = CONFIG,
  sim1 = results_sim1,
  sim2_mse = results_sim2_mse,
  sim2_auc = results_sim2_auc,
  sim3 = results_sim3
), file.path(CONFIG$output_dir, "simulation_results_raw.rds"))

# Create summary markdown
summary_md <- paste0(
  "# cfperformance Simulation Study Results\n\n",
  "Date: ", Sys.Date(), "\n\n",
  "Configuration:\n",
  "- N simulations: ", CONFIG$n_sims, "\n",
  "- N train: ", CONFIG$n_train, "\n",
  "- N test: ", CONFIG$n_test, "\n",
  "- N truth: ", CONFIG$n_truth, "\n",
  "- Bootstrap replications: ", CONFIG$n_boot, "\n",
  "- Confidence level: ", CONFIG$conf_level, "\n\n",

  "## Simulation 1: Counterfactual MSE - Continuous Outcome\n\n",
  "DGP from Boyer, Dahabreh & Steingrimsson (2025)\n\n",
  "```\n",
  paste(capture.output(print(table_sim1, row.names = FALSE)), collapse = "\n"),
  "\n```\n\n",

  "## Simulation 2: Counterfactual MSE - Binary Outcome\n\n",
  "```\n",
  paste(capture.output(print(table_sim2_mse, row.names = FALSE)), collapse = "\n"),
  "\n```\n\n",

  "## Simulation 2B: Counterfactual AUC - Binary Outcome\n\n",
  "```\n",
  paste(capture.output(print(table_sim2_auc, row.names = FALSE)), collapse = "\n"),
  "\n```\n\n",

  "## Simulation 3: Transportability AUC\n\n",
  "DGP from Li, Gatsonis, Dahabreh & Steingrimsson (2022)\n\n",
  "```\n",
  paste(capture.output(print(table_sim3, row.names = FALSE)), collapse = "\n"),
  "\n```\n\n",

  "## Interpretation\n\n",
  "- **Bias**: Difference between mean estimate and true value\n",
  "- **Rel_Bias_Pct**: Relative bias as percentage of true value\n",
  "- **Emp_SE**: Empirical standard error (SD of estimates)\n",
  "- **Mean_SE**: Mean of estimated standard errors\n",
  "- **SE_Ratio**: Mean_SE / Emp_SE (should be ~1 if SEs are calibrated)\n",
  "- **RMSE**: Root mean squared error\n",
  "- **Coverage**: Proportion of CIs containing true value (target: ", CONFIG$conf_level, ")\n",
  "- **CI_Width**: Mean confidence interval width\n"
)

writeLines(summary_md, file.path(CONFIG$output_dir, "simulation_summary.md"))

cat("Results saved to:", CONFIG$output_dir, "\n")
cat("- sim1_continuous_mse.csv\n")
cat("- sim2_binary_mse.csv\n")
cat("- sim2_binary_auc.csv\n")
cat("- sim3_transport_auc.csv\n")
cat("- simulation_results_raw.rds\n")
cat("- simulation_summary.md\n")

cat("\n", strrep("=", 70), "\n")
cat("Simulation Study Complete\n")
cat(strrep("=", 70), "\n")
