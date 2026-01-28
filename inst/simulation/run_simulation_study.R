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
# Use devtools::load_all() if running from source, otherwise library()
if (file.exists("DESCRIPTION")) {
  devtools::load_all(quiet = TRUE)
} else if (file.exists("../../DESCRIPTION")) {
  devtools::load_all("../..", quiet = TRUE)
} else {
  library(cfperformance)
}
library(MASS)
library(parallel)

# Check if grf is available for cross-fitting simulations
if (!requireNamespace("grf", quietly = TRUE)) {
  message("Install 'grf' package for GRF cross-fitting: install.packages('grf')")
  USE_GRF <- FALSE
} else {
  library(grf)
  USE_GRF <- TRUE
}

# Check if running interactively for progress bars
if (!requireNamespace("progress", quietly = TRUE)) {
  message("Install 'progress' package for progress bars: install.packages('progress')")
  USE_PROGRESS <- FALSE
} else {
  library(progress)
  USE_PROGRESS <- TRUE
}

# Configuration
CONFIG <- list(
  # Simulation parameters
  n_sims = 1000,           # Number of simulations (use 1000 for full study)

  n_train = 1000,          # Training sample size
  n_test = 1000,           # Test sample size
  n_truth = 10000,         # Large sample for computing true values (per-sim)
  n_truth_fixed = 1000000, # Very large sample for fixed true value (reduces MC error)

  # Bootstrap parameters for SE estimation
  n_boot = 100,           # Bootstrap replications (use 200 for full study)
  conf_level = 0.95,      # Confidence level
  boot_ci_type = "percentile",  # Bootstrap CI type: "percentile", "normal", or "basic"


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
#' @param true_value True parameter value (scalar for fixed, vector for simulation-specific)
#' @param conf_level Confidence level for coverage
compute_metrics <- function(estimates, ses = NULL, ci_lowers = NULL,
                            ci_uppers = NULL, true_value, conf_level = 0.95) {

  # For basic metrics (bias, RMSE), only filter on valid estimates
  valid_est <- !is.na(estimates)
  n_valid <- sum(valid_est)
  
  if (n_valid == 0) {
    return(list(
      n_sims = 0,
      true_value = NA,
      mean_estimate = NA,
      bias = NA,
      rel_bias_pct = NA,
      empirical_se = NA,
      rmse = NA
    ))
  }
  
  estimates_valid <- estimates[valid_est]
  
  # Handle simulation-specific true values
  if (length(true_value) > 1) {
    true_value_valid <- true_value[valid_est]
    mean_true <- mean(true_value_valid)
  } else {
    true_value_valid <- true_value
    mean_true <- true_value
  }

  metrics <- list(
    n_sims = n_valid,
    true_value = mean_true,
    mean_estimate = mean(estimates_valid),
    bias = mean(estimates_valid) - mean_true,
    rel_bias_pct = (mean(estimates_valid) - mean_true) / abs(mean_true) * 100,
    empirical_se = sd(estimates_valid),
    rmse = sqrt(mean((estimates_valid - mean_true)^2))
  )

  # For SE-based metrics, require valid SEs
  if (!is.null(ses)) {
    valid_se <- valid_est & !is.na(ses)
    if (sum(valid_se) > 0) {
      ses_valid <- ses[valid_se]
      metrics$mean_se = mean(ses_valid)
      # SE ratio uses empirical SE from those same observations
      metrics$se_ratio = mean(ses_valid) / sd(estimates[valid_se])
    }
  }

  # For coverage, require valid CIs
  if (!is.null(ci_lowers) && !is.null(ci_uppers)) {
    valid_ci <- valid_est & !is.na(ci_lowers) & !is.na(ci_uppers)
    if (sum(valid_ci) > 0) {
      ci_lowers_valid <- ci_lowers[valid_ci]
      ci_uppers_valid <- ci_uppers[valid_ci]
      # Use simulation-specific true values for coverage if available
      if (length(true_value) > 1) {
        true_value_ci <- true_value[valid_ci]
      } else {
        true_value_ci <- true_value
      }
      metrics$coverage = mean(ci_lowers_valid <= true_value_ci & ci_uppers_valid >= true_value_ci)
      metrics$mean_ci_width = mean(ci_uppers_valid - ci_lowers_valid)
    }
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
  # Use as.numeric to avoid integer overflow with large samples
  (sum(r[1:n1]) - as.numeric(n1) * (n1 + 1) / 2) / (as.numeric(n1) * n0)
}


#' Create a progress bar
#' @param n Total number of iterations
#' @param name Name to display
create_progress_bar <- function(n, name = "Simulation") {
  if (USE_PROGRESS) {
    progress_bar$new(
      format = paste0("  ", name, " [:bar] :current/:total (:percent) eta: :eta"),
      total = n,
      clear = FALSE,
      width = 60
    )
  } else {
    NULL
  }
}

#' Update progress bar or print status
update_progress <- function(pb, sim, n_sims) {
  if (!is.null(pb)) {
    pb$tick()
  } else if (sim %% 50 == 0) {
    cat(sprintf("  Simulation %d of %d (%.0f%%)\n", sim, n_sims, sim/n_sims*100))
  }
}


#' Run simulations (parallel or sequential)
#' @param n_sims Number of simulations
#' @param dgp Data generating process function
#' @param estimators List of estimator functions
#' @param use_parallel Whether to use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param sim_name Name for progress display
#' @param ... Additional arguments passed to DGP
run_simulations <- function(n_sims, dgp, estimators, use_parallel, n_cores,
                            sim_name = "Simulation", ...) {
  if (use_parallel && n_cores > 1) {
    # Parallel execution with mclapply
    cat(sprintf("  Using %d cores for parallel execution\n", n_cores))
    results_list <- mclapply(1:n_sims, function(sim) {
      run_single_sim(sim, dgp, estimators, ...)
    }, mc.cores = n_cores, mc.set.seed = TRUE)
  } else {
    # Sequential execution with progress bar
    pb <- create_progress_bar(n_sims, sim_name)
    results_list <- lapply(1:n_sims, function(sim) {
      update_progress(pb, sim, n_sims)
      run_single_sim(sim, dgp, estimators, ...)
    })
  }
  do.call(rbind, results_list)
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
#
# NOTE: Training data is generated ONCE and the prediction model is fixed
# across all simulations. This ensures the bootstrap SE (which conditions on
# a fixed prediction model) is estimating the correct source of variability.

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 1: Counterfactual MSE - Continuous Outcome\n")
cat("(Boyer, Dahabreh & Steingrimsson 2025)\n")
cat(strrep("=", 70), "\n\n")

# Generate FIXED training data and prediction model ONCE before simulations
cat("Generating fixed training data and prediction model...\n")
set.seed(CONFIG$seed)  # Ensure reproducibility

# Training data (fixed across all simulations)
X_train_fixed <- runif(CONFIG$n_train, 0, 10)
A_train_fixed <- rbinom(CONFIG$n_train, 1, plogis(-1.5 + 0.3 * X_train_fixed))
Y_train_fixed <- 1 + X_train_fixed + 0.5 * X_train_fixed^2 - 3 * A_train_fixed + rnorm(CONFIG$n_train)

# Fit prediction model ONCE on fixed training data
train_df_fixed <- data.frame(X = X_train_fixed, Y = Y_train_fixed)
pred_model_fixed <- lm(Y ~ X + I(X^2), data = train_df_fixed)

# Compute true MSE using very large truth sample (reduces MC error to ~0.003)
cat(sprintf("  Computing true MSE with n=%d samples...\n", CONFIG$n_truth_fixed))
X_truth_fixed <- runif(CONFIG$n_truth_fixed, 0, 10)
Y_truth_fixed <- 1 + X_truth_fixed + 0.5 * X_truth_fixed^2 + rnorm(CONFIG$n_truth_fixed)
pred_truth_fixed <- predict(pred_model_fixed, newdata = data.frame(X = X_truth_fixed))
true_mse_fixed <- mean((Y_truth_fixed - pred_truth_fixed)^2)

# Estimate MC error in truth calculation
mc_se_truth <- sd((Y_truth_fixed - pred_truth_fixed)^2) / sqrt(CONFIG$n_truth_fixed)
cat(sprintf("  True MSE (fixed): %.4f (MC SE: %.4f)\n", true_mse_fixed, mc_se_truth))

dgp_continuous <- function(n_train, n_test, n_truth, pred_model, true_value) {
  # Test data (varies each simulation - this is what bootstrap estimates SE for)
  X_test <- runif(n_test, 0, 10)
  A_test <- rbinom(n_test, 1, plogis(-1.5 + 0.3 * X_test))
  Y_test <- 1 + X_test + 0.5 * X_test^2 - 3 * A_test + rnorm(n_test)

  # Get predictions using FIXED prediction model
  pred_test <- predict(pred_model, newdata = data.frame(X = X_test))

  # Fit propensity models on test data (these are refit each simulation)
  ps_model_correct <- glm(A_test ~ X_test, family = binomial())  # Correct: linear in X
  ps_model_misspec <- glm(A_test ~ 1, family = binomial())       # Misspecified: intercept only

  list(
    X_test = X_test,
    A_test = A_test,
    Y_test = Y_test,
    pred_test = pred_test,
    ps_model_correct = ps_model_correct,
    ps_model_misspec = ps_model_misspec,
    true_value = true_value  # Fixed true value passed in
  )
}

estimators_sim1 <- list(
  # Naive estimator (biased reference)
  Naive = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "naive",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # CL - Correctly specified outcome model (X, X^2)
  CL_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test, X2 = data$X_test^2),
      treatment_level = 0,
      estimator = "cl",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # CL - Misspecified outcome model (X only)
  CL_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "cl",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Correctly specified propensity model (X)
  IPW_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Misspecified propensity model (intercept only)
  IPW_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both correctly specified
  DR_Both_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test, X2 = data$X_test^2),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Outcome misspecified, PS correct
  DR_OM_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - PS misspecified, Outcome correct
  DR_PS_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test, X2 = data$X_test^2),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both misspecified
  DR_Both_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW with influence function SE (correctly specified)
  IPW_Influence = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  },

  # DR with cross-fitting using grf regression_forest
  DR_CrossFit_GRF = function(data) {
    if (!USE_GRF) {
      return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
    }
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X = data$X_test),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ml_learner("grf", model = "probability_forest"),
      outcome_model = ml_learner("grf", model = "regression_forest"),
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 1
set.seed(CONFIG$seed + 1)  # Different seed for test data (training already used CONFIG$seed)
cat("Running", CONFIG$n_sims, "simulations...\n")
start_time <- Sys.time()

results_sim1 <- run_simulations(
  n_sims = CONFIG$n_sims,
  dgp = dgp_continuous,
  estimators = estimators_sim1,
  use_parallel = CONFIG$use_parallel,
  n_cores = CONFIG$n_cores,
  sim_name = "Sim 1",
  n_train = CONFIG$n_train,  # Not used anymore but kept for compatibility
  n_test = CONFIG$n_test,
  n_truth = CONFIG$n_truth,  # Not used anymore but kept for compatibility
  pred_model = pred_model_fixed,  # Pass fixed prediction model
  true_value = true_mse_fixed     # Pass fixed true value
)

cat(sprintf("\n  Completed in %.1f minutes\n", difftime(Sys.time(), start_time, units = "mins")))

# Compute metrics for each estimator (true value is now fixed across simulations)
table_sim1 <- do.call(rbind, lapply(unique(results_sim1$estimator), function(est) {
  idx <- results_sim1$estimator == est
  m <- compute_metrics(
    estimates = results_sim1$estimate[idx],
    ses = results_sim1$se[idx],
    ci_lowers = results_sim1$ci_lower[idx],
    ci_uppers = results_sim1$ci_upper[idx],
    true_value = results_sim1$true_value[idx],  # Simulation-specific true values
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
#
# NOTE: Training data is generated ONCE and the prediction model is fixed
# across all simulations.

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 2: Counterfactual MSE & AUC - Binary Outcome\n")
cat("(Boyer, Dahabreh & Steingrimsson 2025)\n")
cat(strrep("=", 70), "\n\n")

# Generate FIXED training data and prediction model ONCE before simulations
cat("Generating fixed training data and prediction model...\n")
set.seed(CONFIG$seed + 100)  # Different seed from sim1

# Training data (fixed across all simulations)
X_train_fixed_sim2 <- mvrnorm(CONFIG$n_train, mu = c(0.2, 0, 0.5), Sigma = diag(c(0.2, 0.2, 0.2)))
ps_true_train <- plogis(0.5 - 2 * X_train_fixed_sim2[, 1] + 3 * X_train_fixed_sim2[, 1]^2 + 
                        2 * X_train_fixed_sim2[, 2] - X_train_fixed_sim2[, 3])
A_train_fixed_sim2 <- rbinom(CONFIG$n_train, 1, ps_true_train)
mu_true_train <- plogis(0.2 + 3 * X_train_fixed_sim2[, 1] - 2 * X_train_fixed_sim2[, 1]^2 + 
                        2 * X_train_fixed_sim2[, 2] + X_train_fixed_sim2[, 3] - 2 * A_train_fixed_sim2)
Y_train_fixed_sim2 <- rbinom(CONFIG$n_train, 1, mu_true_train)

# Fit prediction model ONCE on fixed training data (main effects only - slightly misspecified)
train_df_sim2 <- data.frame(X1 = X_train_fixed_sim2[, 1], X2 = X_train_fixed_sim2[, 2],
                            X3 = X_train_fixed_sim2[, 3], Y = Y_train_fixed_sim2)
pred_model_fixed_sim2 <- glm(Y ~ X1 + X2 + X3, data = train_df_sim2, family = binomial())

# Compute true MSE and AUC using very large truth sample (reduces MC error)
cat(sprintf("  Computing true MSE/AUC with n=%d samples...\n", CONFIG$n_truth_fixed))
X_truth_fixed_sim2 <- mvrnorm(CONFIG$n_truth_fixed, mu = c(0.2, 0, 0.5), Sigma = diag(c(0.2, 0.2, 0.2)))
mu_truth <- plogis(0.2 + 3 * X_truth_fixed_sim2[, 1] - 2 * X_truth_fixed_sim2[, 1]^2 + 
                   2 * X_truth_fixed_sim2[, 2] + X_truth_fixed_sim2[, 3])
Y_truth_fixed_sim2 <- rbinom(CONFIG$n_truth_fixed, 1, mu_truth)
truth_df_sim2 <- data.frame(X1 = X_truth_fixed_sim2[, 1], X2 = X_truth_fixed_sim2[, 2],
                            X3 = X_truth_fixed_sim2[, 3])
pred_truth_fixed_sim2 <- predict(pred_model_fixed_sim2, newdata = truth_df_sim2, type = "response")

true_mse_fixed_sim2 <- mean((Y_truth_fixed_sim2 - pred_truth_fixed_sim2)^2)
true_auc_fixed_sim2 <- compute_auc(Y_truth_fixed_sim2, pred_truth_fixed_sim2)

# Estimate MC error
mc_se_mse_sim2 <- sd((Y_truth_fixed_sim2 - pred_truth_fixed_sim2)^2) / sqrt(CONFIG$n_truth_fixed)
cat(sprintf("  True MSE (fixed): %.4f (MC SE: %.5f)\n", true_mse_fixed_sim2, mc_se_mse_sim2))
cat(sprintf("  True AUC (fixed): %.4f\n", true_auc_fixed_sim2))

dgp_binary <- function(n_train, n_test, n_truth, pred_model, true_mse, true_auc) {
  # Test data (varies each simulation - this is what bootstrap estimates SE for)
  X_test <- mvrnorm(n_test, mu = c(0.2, 0, 0.5), Sigma = diag(c(0.2, 0.2, 0.2)))
  ps_true <- plogis(0.5 - 2 * X_test[, 1] + 3 * X_test[, 1]^2 + 2 * X_test[, 2] - X_test[, 3])
  A_test <- rbinom(n_test, 1, ps_true)
  mu_true <- plogis(0.2 + 3 * X_test[, 1] - 2 * X_test[, 1]^2 + 2 * X_test[, 2] + X_test[, 3] - 2 * A_test)
  Y_test <- rbinom(n_test, 1, mu_true)

  # Predictions using FIXED prediction model
  test_df <- data.frame(X1 = X_test[, 1], X2 = X_test[, 2], X3 = X_test[, 3])
  pred_test <- predict(pred_model, newdata = test_df, type = "response")

  # Fit propensity models on test data (these are refit each simulation)
  ps_df <- data.frame(A = A_test, test_df)
  ps_model_correct <- glm(A ~ X1 + I(X1^2) + X2 + X3, data = ps_df, family = binomial())
  ps_model_misspec <- glm(A ~ X1 + X2 + X3, data = ps_df, family = binomial())

  list(
    X_test = X_test,
    A_test = A_test,
    Y_test = Y_test,
    pred_test = pred_test,
    ps_model_correct = ps_model_correct,
    ps_model_misspec = ps_model_misspec,
    true_mse = true_mse,
    true_auc = true_auc,
    true_value = true_mse  # For MSE simulations
  )
}

estimators_sim2_mse <- list(
  # Naive estimator (biased reference)
  Naive = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "naive",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # CL - Correctly specified outcome model (includes X1^2)
  CL_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "cl",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # CL - Misspecified outcome model (no X1^2)
  CL_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "cl",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Correctly specified propensity model (includes X1^2)
  IPW_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Misspecified propensity model (no X1^2)
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
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both correctly specified
  DR_Both_Correct = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Outcome misspecified, PS correct
  DR_OM_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - PS misspecified, Outcome correct
  DR_PS_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both misspecified
  DR_Both_Misspec = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW with influence function SE (correctly specified)
  IPW_Influence = function(data) {
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  },

  # DR with cross-fitting using grf probability_forest (binary outcome)
  DR_CrossFit_GRF = function(data) {
    if (!USE_GRF) {
      return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
    }
    cf_mse(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ml_learner("grf", model = "probability_forest"),
      outcome_model = ml_learner("grf", model = "probability_forest"),
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 2 - MSE
set.seed(CONFIG$seed + 101)  # Different seed for test data
cat("Running", CONFIG$n_sims, "MSE simulations...\n")
start_time <- Sys.time()

results_sim2_mse <- run_simulations(
  n_sims = CONFIG$n_sims,
  dgp = dgp_binary,
  estimators = estimators_sim2_mse,
  use_parallel = CONFIG$use_parallel,
  n_cores = CONFIG$n_cores,
  sim_name = "Sim 2 MSE",
  n_train = CONFIG$n_train,  # Not used anymore but kept for compatibility
  n_test = CONFIG$n_test,
  n_truth = CONFIG$n_truth,  # Not used anymore but kept for compatibility
  pred_model = pred_model_fixed_sim2,
  true_mse = true_mse_fixed_sim2,
  true_auc = true_auc_fixed_sim2
)

cat(sprintf("\n  Completed in %.1f minutes\n", difftime(Sys.time(), start_time, units = "mins")))

# Compute metrics (true value is now fixed across simulations)
table_sim2_mse <- do.call(rbind, lapply(unique(results_sim2_mse$estimator), function(est) {
  idx <- results_sim2_mse$estimator == est
  m <- compute_metrics(
    estimates = results_sim2_mse$estimate[idx],
    ses = results_sim2_mse$se[idx],
    ci_lowers = results_sim2_mse$ci_lower[idx],
    ci_uppers = results_sim2_mse$ci_upper[idx],
    true_value = true_mse_fixed_sim2  # Fixed true value
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

# Modify DGP to return true_auc as true_value (uses same fixed pred_model)
dgp_binary_auc <- function(n_train, n_test, n_truth, pred_model, true_mse, true_auc) {
  data <- dgp_binary(n_train, n_test, n_truth, pred_model, true_mse, true_auc)
  data$true_value <- data$true_auc
  data
}

estimators_sim2_auc <- list(
  # Naive estimator (biased reference)
  Naive = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "naive",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # OM - Correctly specified outcome model (includes X1^2)
  OM_Correct = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "om",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # OM - Misspecified outcome model (no X1^2)
  OM_Misspec = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "om",
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Correctly specified propensity model (includes X1^2)
  IPW_Correct = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW - Misspecified propensity model (no X1^2)
  IPW_Misspec = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both correctly specified
  DR_Both_Correct = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Outcome misspecified, PS correct
  DR_OM_Misspec = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - PS misspecified, Outcome correct
  DR_PS_Misspec = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both misspecified
  DR_Both_Misspec = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = data$ps_model_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IPW with influence function SE (correctly specified)
  IPW_Influence = function(data) {
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X1_sq = data$X_test[, 1]^2,
                              X2 = data$X_test[, 2], X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "ipw",
      propensity_model = data$ps_model_correct,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  },

  # DR with cross-fitting using grf probability_forest
  DR_CrossFit_GRF = function(data) {
    if (!USE_GRF) {
      return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
    }
    cf_auc(
      predictions = data$pred_test,
      outcomes = data$Y_test,
      treatment = data$A_test,
      covariates = data.frame(X1 = data$X_test[, 1], X2 = data$X_test[, 2],
                              X3 = data$X_test[, 3]),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ml_learner("grf", model = "probability_forest"),
      outcome_model = ml_learner("grf", model = "probability_forest"),
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 2B - AUC
set.seed(CONFIG$seed + 102)  # Different seed for test data
cat("Running", CONFIG$n_sims, "AUC simulations...\n")
start_time <- Sys.time()

results_sim2_auc <- run_simulations(
  n_sims = CONFIG$n_sims,
  dgp = dgp_binary_auc,
  estimators = estimators_sim2_auc,
  use_parallel = CONFIG$use_parallel,
  n_cores = CONFIG$n_cores,
  sim_name = "Sim 2B AUC",
  n_train = CONFIG$n_train,  # Not used anymore but kept for compatibility
  n_test = CONFIG$n_test,
  n_truth = CONFIG$n_truth,  # Not used anymore but kept for compatibility
  pred_model = pred_model_fixed_sim2,
  true_mse = true_mse_fixed_sim2,
  true_auc = true_auc_fixed_sim2
)

cat(sprintf("\n  Completed in %.1f minutes\n", difftime(Sys.time(), start_time, units = "mins")))

# Compute metrics (true value is now fixed across simulations)
table_sim2_auc <- do.call(rbind, lapply(unique(results_sim2_auc$estimator), function(est) {
  idx <- results_sim2_auc$estimator == est
  m <- compute_metrics(
    estimates = results_sim2_auc$estimate[idx],
    ses = results_sim2_auc$se[idx],
    ci_lowers = results_sim2_auc$ci_lower[idx],
    ci_uppers = results_sim2_auc$ci_upper[idx],
    true_value = true_auc_fixed_sim2  # Fixed true value
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
#
# NOTE: We use the TRUE probability model for predictions (not a trained model),
# so the true AUC is a population constant computed once before simulations.

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 3: Transportability AUC\n")
cat("(Li, Gatsonis, Dahabreh & Steingrimsson 2022)\n")
cat(strrep("=", 70), "\n\n")

# Compute FIXED true AUC using very large sample (reduces MC error)
cat("Computing fixed true AUC for target population...\n")
set.seed(CONFIG$seed + 200)
X_truth_sim3 <- mvrnorm(CONFIG$n_truth_fixed, mu = rep(0, 3), Sigma = diag(3))
p_source_truth <- plogis(-0.3 + 0.8 * X_truth_sim3[, 1] + 0.8 * X_truth_sim3[, 1]^2 +
                          0.6 * X_truth_sim3[, 2] - 0.4 * X_truth_sim3[, 3])
S_truth <- rbinom(CONFIG$n_truth_fixed, 1, p_source_truth)
mu_truth_sim3 <- plogis(-0.5 + 1.5 * X_truth_sim3[, 1] - 0.8 * X_truth_sim3[, 1]^2 +
                        1.2 * X_truth_sim3[, 2] + 0.6 * X_truth_sim3[, 3])
Y_truth_sim3 <- rbinom(CONFIG$n_truth_fixed, 1, mu_truth_sim3)

# True AUC in target (S=0) - predictions are the true probability function
target_idx_truth <- S_truth == 0
true_auc_fixed_sim3 <- compute_auc(Y_truth_sim3[target_idx_truth], mu_truth_sim3[target_idx_truth])
cat(sprintf("  True AUC (fixed): %.4f (n_target = %d)\n", 
            true_auc_fixed_sim3, sum(target_idx_truth)))

dgp_transport_auc <- function(n_total, true_auc_fixed) {
  # Generate estimation sample (varies each simulation)
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

  # Predictions (use true model - NOT a trained model)
  pred <- mu

  covs <- data.frame(X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])

  # Correctly specified models (fitted on this sample's data)
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
    true_value = true_auc_fixed,  # Fixed true value passed in
    source_auc = compute_auc(Y[S == 1], pred[S == 1])
  )
}

estimators_sim3 <- list(
  # Source only (biased reference)
  Source_Only = function(data) {
    list(estimate = data$source_auc, se = NULL, ci_lower = NULL, ci_upper = NULL)
  },

  # OM - Correctly specified outcome model
  OM_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "om",
      outcome_model = data$om_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # OM - Misspecified outcome model
  OM_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "om",
      outcome_model = data$om_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IOW (IPW) - Both correctly specified
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
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IOW (IPW) - Selection model misspecified
  IOW_Sel_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "ipw",
      selection_model = data$sel_misspec,
      propensity_model = data$ps_model,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both correctly specified
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
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Selection misspecified, Outcome correct
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
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Outcome misspecified, Selection correct
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
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both misspecified
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
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR with cross-fitting using grf probability_forest
  DR_CrossFit_GRF = function(data) {
    if (!USE_GRF) {
      return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
    }
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = data$A,
      source = data$S,
      covariates = data$covs,
      treatment_level = 0,
      analysis = "transport",
      estimator = "dr",
      selection_model = ml_learner("grf", model = "probability_forest"),
      outcome_model = ml_learner("grf", model = "probability_forest"),
      propensity_model = data$ps_model,
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 3
set.seed(CONFIG$seed + 3)
cat("Running", CONFIG$n_sims, "transportability AUC simulations...\n")
start_time <- Sys.time()

results_sim3 <- run_simulations(
  n_sims = CONFIG$n_sims,
  dgp = dgp_transport_auc,
  estimators = estimators_sim3,
  use_parallel = CONFIG$use_parallel,
  n_cores = CONFIG$n_cores,
  sim_name = "Sim 3 Transport",
  n_total = CONFIG$n_test * 2,
  true_auc_fixed = true_auc_fixed_sim3  # Pass fixed true value
)

cat(sprintf("\n  Completed in %.1f minutes\n", difftime(Sys.time(), start_time, units = "mins")))

# Compute metrics (true value is now fixed across simulations)
table_sim3 <- do.call(rbind, lapply(unique(results_sim3$estimator), function(est) {
  idx <- results_sim3$estimator == est
  m <- compute_metrics(
    estimates = results_sim3$estimate[idx],
    ses = results_sim3$se[idx],
    ci_lowers = results_sim3$ci_lower[idx],
    ci_uppers = results_sim3$ci_upper[idx],
    true_value = true_auc_fixed_sim3  # Fixed true value
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 3 Transportability AUC Results:\n")
print(table_sim3, row.names = FALSE)


# =============================================================================
# SIMULATION 4: Factual Transportability AUC (Li et al. 2022)
# =============================================================================
#
# This simulation uses the DGP from:
# Li, Gatsonis, Dahabreh & Steingrimsson (2022). "Estimating the area
# under the ROC curve when transporting a prediction model to a target
# population." Biometrics.
#
# This is FACTUAL transportability (no treatment/counterfactual) - we
# evaluate the AUC of a prediction model in the target population using
# inverse-odds weighting based on the selection model only.
#
# DGP:
#   X = (X1, X2, X3) ~ MVN(mu=(0.2, -1, 2), Sigma=diag(0.1, 0.1, 0.1))
#   S ~ Bernoulli(expit(-4.8 - 3.5*X1 + 5*X1^2 + 2.5*X2 + 2*X3))
#   Y ~ Bernoulli(expit(0.2 + 3*X1 - 5*X1^2 + 2*X2 + X3))
#
# Prediction model (true outcome model):
#   g(X) = expit(0.2 + 3*X1 - 5*X1^2 + 2*X2 + X3)
#
# Correct models: include X1^2 term
# Misspecified models: main effects only (omit X1^2)
#
# Target: AUC in target population (S=0) - no treatment involved

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION 4: Factual Transportability AUC (No Treatment)\n")
cat("(Li, Gatsonis, Dahabreh & Steingrimsson 2022)\n")
cat(strrep("=", 70), "\n\n")

# Compute FIXED true AUC using very large superpopulation
cat("Computing fixed true AUC for target population (factual mode)...\n")
set.seed(CONFIG$seed + 400)

# Generate superpopulation (10,000 observations as in Li et al.)
n_super <- 10000
X_super_sim4 <- mvrnorm(n_super, mu = c(0.2, -1, 2), Sigma = diag(c(0.1, 0.1, 0.1)))

# Selection model: P(S=1|X) - probability of being in source
p_source_sim4 <- plogis(-4.8 - 3.5 * X_super_sim4[, 1] + 5 * X_super_sim4[, 1]^2 + 
                        2.5 * X_super_sim4[, 2] + 2 * X_super_sim4[, 3])
S_super_sim4 <- rbinom(n_super, 1, p_source_sim4)

# Outcome model: P(Y=1|X) - same for source and target (transportability assumption)
mu_super_sim4 <- plogis(0.2 + 3 * X_super_sim4[, 1] - 5 * X_super_sim4[, 1]^2 + 
                        2 * X_super_sim4[, 2] + X_super_sim4[, 3])
Y_super_sim4 <- rbinom(n_super, 1, mu_super_sim4)

# Predictions using TRUE outcome model (as in Li et al.)
pred_super_sim4 <- mu_super_sim4

# True AUC in target (S=0) using very large sample
# Generate target superpopulation separately to reduce MC error
cat(sprintf("  Computing true AUC with n=%d target samples...\n", CONFIG$n_truth_fixed))
X_truth_target_sim4 <- mvrnorm(CONFIG$n_truth_fixed, mu = c(0.2, -1, 2), 
                                Sigma = diag(c(0.1, 0.1, 0.1)))
# For target, we can condition on S=0 by rejection sampling or just generate
# outcomes directly from outcome model
mu_truth_target_sim4 <- plogis(0.2 + 3 * X_truth_target_sim4[, 1] - 
                               5 * X_truth_target_sim4[, 1]^2 + 
                               2 * X_truth_target_sim4[, 2] + X_truth_target_sim4[, 3])
Y_truth_target_sim4 <- rbinom(CONFIG$n_truth_fixed, 1, mu_truth_target_sim4)
pred_truth_target_sim4 <- mu_truth_target_sim4

true_auc_fixed_sim4 <- compute_auc(Y_truth_target_sim4, pred_truth_target_sim4)

# Summary statistics from Li et al.
n_source_super <- sum(S_super_sim4 == 1)
n_target_super <- sum(S_super_sim4 == 0)
prop_y1_source <- mean(Y_super_sim4[S_super_sim4 == 1])
prop_y1_target <- mean(Y_super_sim4[S_super_sim4 == 0])

cat(sprintf("  Superpopulation: n_source=%d, n_target=%d\n", n_source_super, n_target_super))
cat(sprintf("  Proportion Y=1: source=%.1f%%, target=%.1f%%\n", 
            prop_y1_source * 100, prop_y1_target * 100))
cat(sprintf("  True AUC in target (fixed): %.4f\n", true_auc_fixed_sim4))


dgp_factual_transport_auc <- function(n_obs, true_auc_fixed) {
  # Generate superpopulation
  n_super <- 10000
  X_super <- mvrnorm(n_super, mu = c(0.2, -1, 2), Sigma = diag(c(0.1, 0.1, 0.1)))
  
  # Selection model: P(S=1|X)
  p_source <- plogis(-4.8 - 3.5 * X_super[, 1] + 5 * X_super[, 1]^2 + 
                     2.5 * X_super[, 2] + 2 * X_super[, 3])
  S_super <- rbinom(n_super, 1, p_source)
  
  # Outcome model: P(Y=1|X)
  mu_super <- plogis(0.2 + 3 * X_super[, 1] - 5 * X_super[, 1]^2 + 
                     2 * X_super[, 2] + X_super[, 3])
  Y_super <- rbinom(n_super, 1, mu_super)
  
  # Predictions (true model)
  pred_super <- mu_super
  
  # Sample from source (all S=1 observations)
  source_idx <- which(S_super == 1)
  
  # Sample from target with sampling fraction ~0.0325
  target_super_idx <- which(S_super == 0)
  target_sample_size <- round(length(target_super_idx) * 0.0325)
  target_idx <- sample(target_super_idx, min(target_sample_size, length(target_super_idx)))
  
  # Combine samples
  sample_idx <- c(source_idx, target_idx)
  
  X <- X_super[sample_idx, ]
  S <- S_super[sample_idx]
  Y <- Y_super[sample_idx]
  pred <- pred_super[sample_idx]
  
  covs <- data.frame(X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])
  
  # Correctly specified selection model (includes X1^2)
  sel_correct <- glm(I(1 - S) ~ X1 + I(X1^2) + X2 + X3,
                     data = data.frame(S = S, covs), family = binomial())
  
  # Misspecified selection model (omits X1^2)
  sel_misspec <- glm(I(1 - S) ~ X1 + X2 + X3,
                     data = data.frame(S = S, covs), family = binomial())
  
  # Correctly specified outcome model (includes X1^2) - fitted on source data
  om_correct <- glm(Y ~ X1 + I(X1^2) + X2 + X3,
                    data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                    family = binomial())
  
  # Misspecified outcome model (omits X1^2)
  om_misspec <- glm(Y ~ X1 + X2 + X3,
                    data = data.frame(Y = Y[S == 1], covs[S == 1, ]),
                    family = binomial())
  
  list(
    X = X,
    S = S,
    Y = Y,
    pred = pred,
    covs = covs,
    sel_correct = sel_correct,
    sel_misspec = sel_misspec,
    om_correct = om_correct,
    om_misspec = om_misspec,
    true_value = true_auc_fixed,
    source_auc = compute_auc(Y[S == 1], pred[S == 1]),
    n_source = sum(S == 1),
    n_target = sum(S == 0)
  )
}

estimators_sim4 <- list(
  # Source only (biased reference - ignores covariate shift)
  Source_Only = function(data) {
    list(estimate = data$source_auc, se = NULL, ci_lower = NULL, ci_upper = NULL)
  },

  # OM - Correctly specified outcome model (factual mode: no treatment)
  OM_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,  # Factual mode - no treatment
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "om",
      outcome_model = data$om_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # OM - Misspecified outcome model
  OM_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "om",
      outcome_model = data$om_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IOW (IPW) - Correctly specified selection model
  IOW_Sel_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "ipw",
      selection_model = data$sel_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # IOW (IPW) - Misspecified selection model
  IOW_Sel_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "ipw",
      selection_model = data$sel_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both correctly specified
  DR_Both_Correct = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_correct,
      outcome_model = data$om_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Selection misspecified, Outcome correct
  DR_Sel_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_misspec,
      outcome_model = data$om_correct,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Outcome misspecified, Selection correct
  DR_OM_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_correct,
      outcome_model = data$om_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR - Both misspecified
  DR_Both_Misspec = function(data) {
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "dr",
      selection_model = data$sel_misspec,
      outcome_model = data$om_misspec,
      se_method = "bootstrap",
      n_boot = CONFIG$n_boot,
      boot_ci_type = CONFIG$boot_ci_type,
      conf_level = CONFIG$conf_level
    )
  },

  # DR with cross-fitting using grf probability_forest
  DR_CrossFit_GRF = function(data) {
    if (!USE_GRF) {
      return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
    }
    tr_auc(
      predictions = data$pred,
      outcomes = data$Y,
      treatment = NULL,
      source = data$S,
      covariates = data$covs,
      treatment_level = NULL,
      analysis = "transport",
      estimator = "dr",
      selection_model = ml_learner("grf", model = "probability_forest"),
      outcome_model = ml_learner("grf", model = "probability_forest"),
      cross_fit = TRUE,
      n_folds = 5,
      se_method = "influence",
      conf_level = CONFIG$conf_level
    )
  }
)

# Run simulation 4
set.seed(CONFIG$seed + 4)
cat("Running", CONFIG$n_sims, "factual transportability AUC simulations...\n")
start_time <- Sys.time()

results_sim4 <- run_simulations(
  n_sims = CONFIG$n_sims,
  dgp = dgp_factual_transport_auc,
  estimators = estimators_sim4,
  use_parallel = CONFIG$use_parallel,
  n_cores = CONFIG$n_cores,
  sim_name = "Sim 4 Factual",
  n_obs = CONFIG$n_test,  # Not directly used; DGP generates ~1000 obs
  true_auc_fixed = true_auc_fixed_sim4
)

cat(sprintf("\n  Completed in %.1f minutes\n", difftime(Sys.time(), start_time, units = "mins")))

# Compute metrics (true value is now fixed across simulations)
table_sim4 <- do.call(rbind, lapply(unique(results_sim4$estimator), function(est) {
  idx <- results_sim4$estimator == est
  m <- compute_metrics(
    estimates = results_sim4$estimate[idx],
    ses = results_sim4$se[idx],
    ci_lowers = results_sim4$ci_lower[idx],
    ci_uppers = results_sim4$ci_upper[idx],
    true_value = true_auc_fixed_sim4  # Fixed true value
  )
  format_metrics_row(m, est)
}))

cat("\nSimulation 4 Factual Transportability AUC Results:\n")
print(table_sim4, row.names = FALSE)


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
write.csv(table_sim4, file.path(CONFIG$output_dir, "sim4_factual_transport_auc.csv"),
          row.names = FALSE)

# Save raw results for further analysis
saveRDS(list(
  config = CONFIG,
  sim1 = results_sim1,
  sim2_mse = results_sim2_mse,
  sim2_auc = results_sim2_auc,
  sim3 = results_sim3,
  sim4 = results_sim4
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

  "## Simulation 4: Factual Transportability AUC (No Treatment)\n\n",
  "DGP from Li, Gatsonis, Dahabreh & Steingrimsson (2022)\n",
  "Evaluates prediction model AUC in target population without treatment/counterfactual.\n\n",
  "```\n",
  paste(capture.output(print(table_sim4, row.names = FALSE)), collapse = "\n"),
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
cat("- sim4_factual_transport_auc.csv\n")
cat("- simulation_results_raw.rds\n")
cat("- simulation_summary.md\n")

cat("\n", strrep("=", 70), "\n")
cat("Simulation Study Complete\n")
cat(strrep("=", 70), "\n")
