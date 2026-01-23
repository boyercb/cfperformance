#' Variance Estimation for Counterfactual Performance Measures
#'
#' Functions for computing standard errors and confidence intervals using
#' influence function-based and bootstrap methods.
#'
#' @name variance
#' @keywords internal
NULL


#' Influence Function-Based Standard Errors for MSE
#'
#' Computes standard errors using the efficient influence function for
#' the doubly robust MSE estimator.
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of observed outcomes.
#' @param treatment Numeric vector of treatment indicators.
#' @param covariates Matrix or data frame of covariates.
#' @param treatment_level Counterfactual treatment level.
#' @param estimator Character string specifying estimator type.
#' @param propensity_model Fitted propensity model.
#' @param outcome_model Fitted outcome model.
#'
#' @return Standard error estimate.
#'
#' @details
#' The efficient influence function for the MSE under treatment level \eqn{a} is:
#' \deqn{\chi(O) = \frac{I(A=a)}{e_a(X)}[L(Y, \mu(X^*)) - h_a(X)] + h_a(X) - \psi(a)}
#'
#' where \eqn{e_a(X) = P(A=a|X)} is the propensity score, \eqn{h_a(X) = E[L|X,A=a]}
#' is the conditional loss, and \eqn{\psi(a)} is the target estimand.
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @keywords internal
.influence_se_mse <- function(predictions, outcomes, treatment, covariates,
                               treatment_level, estimator, propensity_model,
                               outcome_model) {

  n <- length(outcomes)
  loss <- (outcomes - predictions)^2

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Treatment indicator
  I_a <- as.numeric(treatment == treatment_level)

  if (estimator == "naive") {
    # For naive estimator, influence function is just centered loss
    psi_hat <- mean(loss)
    phi <- loss - psi_hat
    se <- sqrt(var(phi) / n)
    return(se)
  }

  # Get propensity scores
  ps <- predict(propensity_model, newdata = covariates, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps  # P(A = 0 | X)
  }
  ps <- pmax(pmin(ps, 0.99), 0.01)  # Truncate for stability

  # Get conditional loss predictions
  pY <- predict(outcome_model, newdata = covariates, type = "response")
  h <- pY - 2 * predictions * pY + predictions^2

  if (estimator == "cl") {
    # Conditional loss estimator
    psi_hat <- mean(h)
    phi <- h - psi_hat
    se <- sqrt(var(phi) / n)
    return(se)

  } else if (estimator == "ipw") {
    # IPW estimator - use weighted variance
    psi_hat <- sum(I_a / ps * loss) / sum(I_a)
    phi <- I_a / ps * loss - psi_hat
    se <- sqrt(var(phi) / n)
    return(se)

  } else if (estimator == "dr") {
    # Doubly robust estimator
    # Efficient influence function:
    # phi = (I_a / e_a) * (L - h) + h - psi
    augmentation <- I_a / ps * (loss - h)
    psi_hat <- mean(h + augmentation)

    # Influence function for each observation
    phi <- h + augmentation - psi_hat

    # Note: This is the "plug-in" variance that ignores uncertainty
    # in nuisance function estimation. With sample splitting/cross-fitting,
    # this provides valid inference even without nuisance function correction.
    se <- sqrt(var(phi) / n)
    return(se)
  }

  return(NA_real_)
}


#' Bootstrap Standard Errors for MSE
#'
#' Computes bootstrap standard errors and confidence intervals for
#' counterfactual MSE estimators.
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of observed outcomes.
#' @param treatment Numeric vector of treatment indicators.
#' @param covariates Matrix or data frame of covariates.
#' @param treatment_level Counterfactual treatment level.
#' @param estimator Character string specifying estimator type.
#' @param n_boot Number of bootstrap replications.
#' @param conf_level Confidence level.
#' @param parallel Logical indicating whether to use parallel processing.
#' @param ncores Number of cores for parallel processing.
#' @param ... Additional arguments.
#'
#' @return List containing:
#'   \item{se}{Bootstrap standard error}
#'   \item{ci_lower}{Lower confidence interval bound}
#'   \item{ci_upper}{Upper confidence interval bound}
#'   \item{boot_estimates}{Vector of bootstrap estimates}
#'
#' @keywords internal
.bootstrap_mse <- function(predictions, outcomes, treatment, covariates,
                           treatment_level, estimator, n_boot = 500,
                           conf_level = 0.95, parallel = FALSE,
                           ncores = NULL, ...) {

  n <- length(outcomes)
  alpha <- 1 - conf_level

  # Define single bootstrap iteration
  boot_one <- function(idx) {
    # Refit nuisance models on bootstrap sample
    nuisance_b <- .fit_nuisance_models(
      treatment = treatment[idx],
      outcomes = outcomes[idx],
      covariates = covariates[idx, , drop = FALSE],
      treatment_level = treatment_level,
      propensity_model = NULL,
      outcome_model = NULL
    )

    .compute_mse(
      predictions = predictions[idx],
      outcomes = outcomes[idx],
      treatment = treatment[idx],
      covariates = covariates[idx, , drop = FALSE],
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = nuisance_b$propensity,
      outcome_model = nuisance_b$outcome
    )
  }

  # Generate bootstrap samples
  boot_indices <- lapply(seq_len(n_boot), function(b) {
    sample(n, n, replace = TRUE)
  })

  # Run bootstrap
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    ncores <- ncores %||% (parallel::detectCores() - 1)
    boot_estimates <- parallel::mclapply(boot_indices, boot_one,
                                          mc.cores = ncores)
    boot_estimates <- unlist(boot_estimates)
  } else {
    boot_estimates <- vapply(boot_indices, boot_one, numeric(1))
  }

  # Remove any NA/NaN values
  boot_estimates <- boot_estimates[is.finite(boot_estimates)]

  list(
    se = sd(boot_estimates),
    ci_lower = unname(quantile(boot_estimates, alpha / 2)),
    ci_upper = unname(quantile(boot_estimates, 1 - alpha / 2)),
    boot_estimates = boot_estimates
  )
}


#' Bootstrap Standard Errors for AUC
#'
#' Computes bootstrap standard errors and confidence intervals for
#' counterfactual AUC estimators.
#'
#' @inheritParams .bootstrap_mse
#'
#' @return List containing bootstrap results.
#'
#' @keywords internal
.bootstrap_auc <- function(predictions, outcomes, treatment, covariates,
                           treatment_level, estimator, n_boot = 500,
                           conf_level = 0.95, parallel = FALSE,
                           ncores = NULL, ...) {

  n <- length(outcomes)
  alpha <- 1 - conf_level

  boot_one <- function(idx) {
    nuisance_b <- .fit_nuisance_models(
      treatment = treatment[idx],
      outcomes = outcomes[idx],
      covariates = covariates[idx, , drop = FALSE],
      treatment_level = treatment_level,
      propensity_model = NULL,
      outcome_model = NULL
    )

    .compute_auc(
      predictions = predictions[idx],
      outcomes = outcomes[idx],
      treatment = treatment[idx],
      covariates = covariates[idx, , drop = FALSE],
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = nuisance_b$propensity,
      outcome_model = nuisance_b$outcome
    )
  }

  boot_indices <- lapply(seq_len(n_boot), function(b) {
    sample(n, n, replace = TRUE)
  })

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    ncores <- ncores %||% (parallel::detectCores() - 1)
    boot_estimates <- parallel::mclapply(boot_indices, boot_one,
                                          mc.cores = ncores)
    boot_estimates <- unlist(boot_estimates)
  } else {
    boot_estimates <- vapply(boot_indices, boot_one, numeric(1))
  }

  boot_estimates <- boot_estimates[is.finite(boot_estimates)]

  list(
    se = sd(boot_estimates),
    ci_lower = unname(quantile(boot_estimates, alpha / 2)),
    ci_upper = unname(quantile(boot_estimates, 1 - alpha / 2)),
    boot_estimates = boot_estimates
  )
}


#' Influence Function-Based Standard Errors for AUC
#'
#' Computes standard errors for counterfactual AUC using influence functions.
#'
#' @inheritParams .influence_se_mse
#'
#' @return Standard error estimate.
#'
#' @details
#' The influence function for the AUC is based on the U-statistic representation
#' of the concordance probability. See Li et al. (2024) for details on the
#' influence function-based variance estimator for counterfactual AUC.
#'
#' @references
#' Li, B., Steingrimsson, J. A., & Dahabreh, I. J. (2024).
#' "Efficient inference for counterfactual area under the ROC curve."
#'
#' @keywords internal
.influence_se_auc <- function(predictions, outcomes, treatment, covariates,
                               treatment_level, estimator, propensity_model,
                               outcome_model) {

  n <- length(outcomes)

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Treatment indicator
  I_a <- as.numeric(treatment == treatment_level)

  if (estimator == "naive") {
    # Naive AUC - use DeLong's method
    se <- .delong_se(predictions, outcomes)
    return(se)
  }

  # Get propensity scores
  ps <- predict(propensity_model, newdata = covariates, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }
  ps <- pmax(pmin(ps, 0.99), 0.01)

  # Get outcome predictions
  pY <- predict(outcome_model, newdata = covariates, type = "response")

  # For AUC, we need to compute the influence function based on the

  # concordance probability estimator. This is more complex than MSE
  # due to the pairwise nature of the estimand.

  # Simplified approach: use a jackknife-like variance estimator
  # based on the efficient influence function structure

  if (estimator == "dr") {
    # Compute AUC estimate first
    auc_hat <- .compute_auc(predictions, outcomes, treatment, covariates,
                            treatment_level, estimator, propensity_model,
                            outcome_model)

    # Compute influence function contributions
    # This is an approximation based on the structure of the DR estimator

    # Cases and controls under counterfactual treatment
    cases <- which(outcomes == 1 & I_a == 1)
    controls <- which(outcomes == 0 & I_a == 1)

    n1 <- length(cases)
    n0 <- length(controls)

    if (n1 < 2 || n0 < 2) {
      warning("Too few cases or controls for influence function SE")
      return(NA_real_)
    }

    # Weight contributions
    w_cases <- 1 / ps[cases]
    w_controls <- 1 / ps[controls]

    # Compute pairwise concordance contributions
    phi <- numeric(n)

    for (i in cases) {
      # Contribution from case i
      for (j in controls) {
        conc <- as.numeric(predictions[i] > predictions[j]) +
          0.5 * as.numeric(predictions[i] == predictions[j])
        phi[i] <- phi[i] + (1/ps[i]) * (conc - auc_hat) / (n1 * n0)
        phi[j] <- phi[j] + (1/ps[j]) * (conc - auc_hat) / (n1 * n0)
      }
    }

    # Add outcome model contributions for non-treated
    non_treated <- which(I_a == 0)
    for (i in non_treated) {
      # Expected concordance contribution
      for (j in controls) {
        if (I_a[j] == 1) {
          conc_exp <- pY[i] * (1 - pY[j]) + 0.5 * pY[i] * pY[j]
          phi[i] <- phi[i] + (conc_exp - auc_hat * pY[i]) / (n1 * n0)
        }
      }
    }

    se <- sqrt(var(phi) / n)
    return(se)

  } else if (estimator == "ipw") {
    # IPW variance
    auc_hat <- .compute_auc(predictions, outcomes, treatment, covariates,
                            treatment_level, "ipw", propensity_model,
                            outcome_model)

    cases <- which(outcomes == 1 & I_a == 1)
    controls <- which(outcomes == 0 & I_a == 1)
    n1 <- length(cases)
    n0 <- length(controls)

    if (n1 < 2 || n0 < 2) {
      return(NA_real_)
    }

    phi <- numeric(n)
    for (i in cases) {
      for (j in controls) {
        conc <- as.numeric(predictions[i] > predictions[j]) +
          0.5 * as.numeric(predictions[i] == predictions[j])
        phi[i] <- phi[i] + (1/ps[i]) * (conc - auc_hat) / (n1 * n0)
        phi[j] <- phi[j] + (1/ps[j]) * (conc - auc_hat) / (n1 * n0)
      }
    }

    se <- sqrt(var(phi) / n)
    return(se)

  } else if (estimator == "cl") {
    # Outcome model estimator
    auc_hat <- .compute_auc(predictions, outcomes, treatment, covariates,
                            treatment_level, "cl", propensity_model,
                            outcome_model)

    # Influence function based on expected concordance
    phi <- numeric(n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) {
          conc_exp <- pY[i] * (1 - pY[j]) + 0.5 * pY[i] * pY[j]
          phi[i] <- phi[i] + (conc_exp - auc_hat) / (n * (n - 1))
        }
      }
    }

    se <- sqrt(var(phi) / n)
    return(se)
  }

  return(NA_real_)
}


#' DeLong's Method for AUC Standard Error
#'
#' Computes the standard error of the empirical AUC using DeLong's method.
#'
#' @param predictions Numeric vector of predictions.
#' @param outcomes Binary outcome vector.
#'
#' @return Standard error of AUC.
#'
#' @references
#' DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988).
#' "Comparing the areas under two or more correlated receiver operating
#' characteristic curves: a nonparametric approach."
#' *Biometrics*, 44(3), 837-845.
#'
#' @keywords internal
.delong_se <- function(predictions, outcomes) {
  cases <- predictions[outcomes == 1]
  controls <- predictions[outcomes == 0]

  n1 <- length(cases)
  n0 <- length(controls)

  if (n1 < 2 || n0 < 2) {
    return(NA_real_)
  }

  # Compute placement values (Mann-Whitney statistics)
  V10 <- numeric(n1)  # For each case, proportion of controls with lower score
  V01 <- numeric(n0)  # For each control, proportion of cases with higher score

  for (i in seq_len(n1)) {
    V10[i] <- mean(as.numeric(cases[i] > controls) +
                     0.5 * as.numeric(cases[i] == controls))
  }

  for (j in seq_len(n0)) {
    V01[j] <- mean(as.numeric(controls[j] < cases) +
                     0.5 * as.numeric(controls[j] == cases))
  }

  # DeLong variance
  S10 <- var(V10)
  S01 <- var(V01)

  var_auc <- S10 / n1 + S01 / n0
  se <- sqrt(var_auc)

  return(se)
}


#' Cross-Fitting for Nuisance Model Estimation
#'
#' Implements K-fold cross-fitting for nuisance model estimation to enable
#' valid inference with flexible machine learning estimators.
#'
#' @param treatment Numeric vector of treatment indicators.
#' @param outcomes Numeric vector of observed outcomes.
#' @param covariates Matrix or data frame of covariates.
#' @param treatment_level Counterfactual treatment level.
#' @param predictions Numeric vector of model predictions.
#' @param K Number of folds for cross-fitting (default: 5).
#' @param propensity_formula Formula for propensity model (optional).
#' @param outcome_formula Formula for outcome model (optional).
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return List containing:
#'   \item{ps}{Cross-fitted propensity scores}
#'   \item{h}{Cross-fitted conditional loss predictions}
#'   \item{folds}{Fold assignments}
#'
#' @details
#' Cross-fitting (sample splitting) allows the use of flexible machine learning
#' methods for nuisance function estimation while maintaining valid inference.
#' Each observation's nuisance function predictions are made using models
#' trained on data excluding that observation's fold.
#'
#' When `propensity_learner` or `outcome_learner` is provided as an `ml_learner`
#' object, the corresponding ML method is used for cross-fitted estimation.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., et al. (2018).
#' "Double/debiased machine learning for treatment and structural parameters."
#' *The Econometrics Journal*, 21(1), C1-C68.
#'
#' @keywords internal
.cross_fit_nuisance <- function(treatment, outcomes, covariates,
                                 treatment_level, predictions,
                                 K = 5, propensity_formula = NULL,
                                 outcome_formula = NULL,
                                 propensity_learner = NULL,
                                 outcome_learner = NULL,
                                 parallel = FALSE,
                                 ...) {

  n <- length(outcomes)

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Create fold assignments
  folds <- sample(rep(1:K, length.out = n))

  # Initialize output vectors
  ps_cf <- numeric(n)  # Cross-fitted propensity scores
  h_cf <- numeric(n)   # Cross-fitted conditional loss predictions

  for (k in 1:K) {
    # Training and validation indices
    train_idx <- which(folds != k)
    val_idx <- which(folds == k)

    # Fit propensity model on training fold
    ps_data <- cbind(A = treatment[train_idx], covariates[train_idx, , drop = FALSE])

    if (!is.null(propensity_learner) && is_ml_learner(propensity_learner)) {
      # Use ML learner
      ps_formula <- if (!is.null(propensity_formula)) propensity_formula else A ~ .
      ps_model <- .fit_ml_learner(propensity_learner, ps_formula,
                                   data = ps_data, family = "binomial")
      ps_pred <- .predict_ml_learner(ps_model, covariates[val_idx, , drop = FALSE])
    } else if (is.null(propensity_formula)) {
      ps_model <- glm(A ~ ., data = ps_data, family = binomial())
      ps_pred <- predict(ps_model, newdata = covariates[val_idx, , drop = FALSE],
                         type = "response")
    } else {
      ps_model <- glm(propensity_formula, data = ps_data, family = binomial())
      ps_pred <- predict(ps_model, newdata = covariates[val_idx, , drop = FALSE],
                         type = "response")
    }

    if (treatment_level == 0) {
      ps_pred <- 1 - ps_pred
    }
    ps_cf[val_idx] <- pmax(pmin(ps_pred, 0.99), 0.01)

    # Fit outcome model on training fold (among treated)
    subset_train <- train_idx[treatment[train_idx] == treatment_level]
    outcome_data <- cbind(Y = outcomes[subset_train],
                          covariates[subset_train, , drop = FALSE])

    if (!is.null(outcome_learner) && is_ml_learner(outcome_learner)) {
      # Use ML learner
      out_formula <- if (!is.null(outcome_formula)) outcome_formula else Y ~ .
      outcome_model <- .fit_ml_learner(outcome_learner, out_formula,
                                        data = outcome_data, family = "binomial")
      pY <- .predict_ml_learner(outcome_model, covariates[val_idx, , drop = FALSE])
    } else if (is.null(outcome_formula)) {
      outcome_model <- glm(Y ~ ., data = outcome_data, family = binomial())
      pY <- predict(outcome_model, newdata = covariates[val_idx, , drop = FALSE],
                    type = "response")
    } else {
      outcome_model <- glm(outcome_formula, data = outcome_data, family = binomial())
      pY <- predict(outcome_model, newdata = covariates[val_idx, , drop = FALSE],
                    type = "response")
    }

    # Compute conditional loss: E[(Y - pred)^2 | X, A=a] = p(1-p) + (p - pred)^2
    # For binary Y: E[Y^2] = p, so E[(Y - pred)^2] = p - 2*p*pred + pred^2
    h_cf[val_idx] <- pY - 2 * predictions[val_idx] * pY + predictions[val_idx]^2
  }

  list(
    ps = ps_cf,
    h = h_cf,
    folds = folds
  )
}


#' Compute DR MSE with Cross-Fitting
#'
#' Computes the doubly robust MSE estimator using cross-fitted nuisance functions.
#'
#' @inheritParams cf_mse
#' @param K Number of folds for cross-fitting.
#'
#' @return Doubly robust MSE estimate with cross-fitting.
#'
#' @keywords internal
.compute_mse_crossfit <- function(predictions, outcomes, treatment, covariates,
                                   treatment_level, K = 5,
                                   propensity_learner = NULL,
                                   outcome_learner = NULL,
                                   parallel = FALSE,
                                   ...) {

  n <- length(outcomes)
  loss <- (outcomes - predictions)^2

  # Get cross-fitted nuisance functions
  cf_nuisance <- .cross_fit_nuisance(
    treatment = treatment,
    outcomes = outcomes,
    covariates = covariates,
    treatment_level = treatment_level,
    predictions = predictions,
    K = K,
    propensity_learner = propensity_learner,
    outcome_learner = outcome_learner,
    parallel = parallel,
    ...
  )

  ps <- cf_nuisance$ps
  h <- cf_nuisance$h

  # Treatment indicator
  I_a <- as.numeric(treatment == treatment_level)

  # DR estimator with cross-fitted nuisance
  augmentation <- I_a / ps * (loss - h)
  estimate <- mean(h + augmentation)

  # Influence function for SE
  phi <- h + augmentation - estimate
  se <- sqrt(var(phi) / n)

  list(
    estimate = estimate,
    se = se,
    ps = ps,
    h = h,
    folds = cf_nuisance$folds
  )
}
