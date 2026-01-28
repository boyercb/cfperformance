#' Estimate Transportable Calibration in the Target Population
#'
#' Estimates the calibration of a prediction model in a target population
#' using data transported from a source population. Supports both
#' **counterfactual** (under hypothetical intervention) and **factual**
#' (observational) prediction model transportability.
#'
#' @param predictions Numeric vector of model predictions (typically probabilities).
#' @param outcomes Numeric vector of observed outcomes (must be binary 0/1).
#' @param treatment Numeric vector of treatment indicators (0/1), or `NULL` for
#'   factual prediction model transportability (no treatment/intervention).
#'   When `NULL`, only the selection model is used for weighting.
#' @param source Numeric vector of population indicators (1=source/RCT, 0=target).
#' @param covariates A matrix or data frame of baseline covariates.
#' @param treatment_level The treatment level of interest (default: `NULL`).
#'   Required when `treatment` is provided; should be `NULL` when `treatment`
#'   is `NULL` (factual mode).
#' @param analysis Character string specifying the type of analysis:
#'   - `"transport"`: Use source outcomes for target estimation (default)
#'   - `"joint"`: Pool source and target data
#' @param estimator Character string specifying the estimator:
#'   - `"dr"`: Doubly robust estimator (default)
#'   - `"ipw"`: Inverse probability weighting estimator
#'   - `"om"`: Outcome model estimator
#' @param selection_model Optional fitted selection model for P(S=0|X). If NULL,
#'   a logistic regression model is fit using the covariates.
#' @param propensity_model Optional fitted propensity score model for P(A=1|X,S=1).
#'   If NULL, a logistic regression model is fit using source data.
#' @param outcome_model Optional fitted outcome model for E\[Y|X,A,S\].
#'   If NULL, a regression model is fit using the relevant data.
#' @param smoother Smoothing method for the calibration curve:
#'   - `"loess"`: Local polynomial regression (default)
#'   - `"binned"`: Binned calibration
#' @param n_bins Number of bins for binned calibration (default: 10).
#' @param span Span parameter for LOESS smoothing (default: 0.75).
#' @param se_method Method for standard error estimation:
#'   - `"none"`: No standard error estimation (default, fastest)
#'   - `"bootstrap"`: Bootstrap standard errors
#' @param n_boot Number of bootstrap replications (default: 200).
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param stratified_boot Logical indicating whether to use stratified bootstrap
#'   that preserves the source/target ratio (default: TRUE).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class `c("tr_calibration", "tr_performance")` containing:
#'   \item{predicted}{Vector of predicted probabilities}
#'   \item{observed}{Vector of smoothed observed probabilities}
#'   \item{smoother}{Smoothing method used}
#'   \item{ici}{Integrated calibration index}
#'   \item{e50}{Median absolute calibration error}
#'   \item{e90}{90th percentile absolute calibration error}
#'   \item{emax}{Maximum absolute calibration error}
#'   \item{se}{Named list of standard errors for calibration metrics (if computed)}
#'   \item{ci_lower}{Named list of lower confidence bounds}
#'   \item{ci_upper}{Named list of upper confidence bounds}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{n_target}{Number of target observations}
#'   \item{n_source}{Number of source observations}
#'   \item{treatment_level}{Treatment level (NULL for factual mode)}
#'
#' @details
#' This function implements estimators for transporting prediction model
#' calibration from a source population to a target population. It supports
#' two modes:
#'
#' ## Counterfactual Mode (treatment provided)
#' When `treatment` is specified, estimates calibration for counterfactual
#' outcomes under a hypothetical intervention. Requires selection, propensity,
#' and outcome models.
#'
#' ## Factual Mode (treatment = NULL)
#' When `treatment` is `NULL`, estimates calibration for observed outcomes in
#' the target population using only the selection model for inverse-odds
#' weighting. This is appropriate for factual prediction model
#' transportability.
#'
#' ## Analysis Types
#' **Transportability Analysis** (`analysis = "transport"`): Uses outcome data
#' from the source population to estimate calibration in the target population.
#'
#' **Joint Analysis** (`analysis = "joint"`): Pools source and target data
#' to estimate calibration in the target population.
#'
#' ## Calibration Metrics
#' The function computes several calibration metrics:
#' - **ICI** (Integrated Calibration Index): Mean absolute difference between
#'   predicted and observed probabilities
#' - **E50**: Median absolute calibration error
#' - **E90**: 90th percentile absolute calibration error
#' - **Emax**: Maximum absolute calibration error
#'
#' For observational analysis (single population), use [cf_calibration()] instead.
#'
#' @references
#' Steingrimsson, J. A., et al. (2023). "Transporting a Prediction Model for
#' Use in a New Target Population." *American Journal of Epidemiology*,
#' 192(2), 296-304. \doi{10.1093/aje/kwac128}
#'
#' Voter, S. R., et al. (2025). "Transportability of machine learning-based
#' counterfactual prediction models with application to CASS."
#' *Diagnostic and Prognostic Research*, 9(4).
#' \doi{10.1186/s41512-025-00201-y}
#'
#' Austin, P. C., & Steyerberg, E. W. (2019). "The Integrated Calibration
#' Index (ICI) and related metrics for quantifying the calibration of
#' logistic regression models."
#' *Statistics in Medicine*, 38(21), 4051-4065.
#'
#' @seealso [cf_calibration()], [tr_mse()], [tr_auc()], [plot.tr_calibration()]
#'
#' @export
#'
#' @examples
#' # Generate example data with source (RCT) and target populations
#' set.seed(123)
#' n <- 1000
#' # Covariates
#' x <- rnorm(n)
#' # Source indicator (S=1 for RCT, S=0 for target)
#' s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
#' # Treatment (randomized in source, confounded in target)
#' a <- ifelse(s == 1,
#'             rbinom(n, 1, 0.5),  # Randomized in RCT
#'             rbinom(n, 1, plogis(-0.5 + 0.5 * x)))  # Confounded in target
#' # Outcome (binary)
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' # Prediction model
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate transportable calibration with IPW
#' result <- tr_calibration(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   analysis = "transport",
#'   estimator = "ipw",
#'   se_method = "none"  # Use "bootstrap" for SEs
#' )
#' print(result)
tr_calibration <- function(predictions,
                           outcomes,
                           treatment = NULL,
                           source,
                           covariates,
                           treatment_level = NULL,
                           analysis = c("transport", "joint"),
                           estimator = c("dr", "ipw", "om"),
                           selection_model = NULL,
                           propensity_model = NULL,
                           outcome_model = NULL,
                           smoother = c("loess", "binned"),
                           n_bins = 10,
                           span = 0.75,
                           se_method = c("none", "bootstrap"),
                           n_boot = 200,
                           conf_level = 0.95,
                           stratified_boot = TRUE,
                           ...) {

  # Match arguments
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  smoother <- match.arg(smoother)
  se_method <- match.arg(se_method)

  # Validate inputs
  .validate_tr_calibration_inputs(predictions, outcomes, treatment, source, covariates, treatment_level)
  
  # Determine mode: factual (no treatment) or counterfactual
  factual_mode <- is.null(treatment)
  
  # Default treatment_level to 1 for backward compatibility if treatment is provided
  if (!factual_mode && is.null(treatment_level)) {
    treatment_level <- 1
  }

  # Compute calibration
  result <- .compute_tr_calibration(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    source = source,
    covariates = covariates,
    treatment_level = treatment_level,
    analysis = analysis,
    estimator = estimator,
    selection_model = selection_model,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    smoother = smoother,
    n_bins = n_bins,
    span = span,
    ...
  )

  # Add standard errors if requested
  if (se_method == "bootstrap") {
    boot_results <- .bootstrap_tr_calibration(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = estimator,
      smoother = smoother,
      n_bins = n_bins,
      span = span,
      n_boot = n_boot,
      stratified = stratified_boot,
      conf_level = conf_level,
      eval_points = result$predicted,  # Evaluate at same points as main estimate
      ...
    )

    # Standard errors for metrics
    result$se <- list(
      ici = boot_results$se_ici,
      e50 = boot_results$se_e50,
      e90 = boot_results$se_e90,
      emax = boot_results$se_emax
    )

    # Confidence intervals (percentile method)
    result$ci_lower <- list(
      ici = boot_results$ci_ici[1],
      e50 = boot_results$ci_e50[1],
      e90 = boot_results$ci_e90[1],
      emax = boot_results$ci_emax[1]
    )
    result$ci_upper <- list(
      ici = boot_results$ci_ici[2],
      e50 = boot_results$ci_e50[2],
      e90 = boot_results$ci_e90[2],
      emax = boot_results$ci_emax[2]
    )

    result$boot_curves <- boot_results$curves
    result$boot_samples <- boot_results$samples
  } else {
    result$se <- NULL
    result$ci_lower <- NULL
    result$ci_upper <- NULL
    result$boot_curves <- NULL
    result$boot_samples <- NULL
  }

  # Add metadata
  result$metric <- "calibration"
  result$analysis <- analysis
  result$estimator <- estimator
  result$treatment_level <- treatment_level
  result$conf_level <- conf_level
  result$se_method <- se_method
  result$call <- match.call()

  class(result) <- c("tr_calibration", "tr_performance")
  return(result)
}


#' Validate Inputs for Transportability Calibration
#' @keywords internal
#' @noRd
.validate_tr_calibration_inputs <- function(predictions, outcomes, treatment,
                                            source, covariates, treatment_level = NULL) {
  n <- length(outcomes)

  # Check lengths match
  if (length(predictions) != n) {
    stop("'predictions' must have the same length as 'outcomes'")
  }
  if (length(source) != n) {
    stop("'source' must have the same length as 'outcomes'")
  }
  if (NROW(covariates) != n) {
    stop("'covariates' must have the same number of rows as length of 'outcomes'")
  }

  # Check binary outcomes
  if (!all(outcomes[!is.na(outcomes)] %in% c(0, 1))) {
    stop("Calibration requires binary outcomes (0/1)")
  }

  # Check source indicator
  if (!all(source %in% c(0, 1))) {
    stop("'source' must be a binary indicator (0=target, 1=source)")
  }
  
  # Factual mode validation: both treatment and treatment_level should be NULL
  # Counterfactual mode: treatment must be provided, treatment_level can default to 1
  if (is.null(treatment) && !is.null(treatment_level)) {
    stop("If treatment is NULL (factual mode), treatment_level must also be NULL")
  }
  
  # If treatment provided, validate it
  if (!is.null(treatment)) {
    if (length(treatment) != n) {
      stop("'treatment' must have the same length as 'outcomes'")
    }
    if (!all(treatment %in% c(0, 1))) {
      stop("'treatment' must be a binary indicator (0/1)")
    }
  }

  # Check predictions are probabilities
  if (any(predictions < 0 | predictions > 1, na.rm = TRUE)) {
    stop("'predictions' must be probabilities between 0 and 1")
  }

  # Check for sufficient data
  n_source <- sum(source == 1)
  n_target <- sum(source == 0)

  if (n_source < 10) {
    stop("Insufficient source data (need at least 10 observations)")
  }
  if (n_target < 10) {
    stop("Insufficient target data (need at least 10 observations)")
  }

  invisible(TRUE)
}


#' Compute Transportable Calibration
#' @keywords internal
#' @noRd
.compute_tr_calibration <- function(predictions, outcomes, treatment, source,
                                    covariates, treatment_level, analysis,
                                    estimator, selection_model, propensity_model,
                                    outcome_model, smoother, n_bins, span, ...) {

  n <- length(outcomes)
  n_source <- sum(source == 1)
  n_target <- sum(source == 0)

  # Determine if factual mode (no treatment)
  factual_mode <- is.null(treatment)
  
  # Indicator for treatment level (or all observations in factual mode)
  if (factual_mode) {
    I_a <- rep(TRUE, n)
  } else {
    I_a <- treatment == treatment_level
  }

  # Fit nuisance models
  nuisance <- .fit_calibration_nuisance(
    outcomes = outcomes,
    treatment = treatment,
    source = source,
    covariates = covariates,
    treatment_level = treatment_level,
    analysis = analysis,
    estimator = estimator,
    selection_model = selection_model,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    predictions = predictions
  )

  # Dispatch to appropriate analysis type
  if (analysis == "transport") {
    calib_result <- .tr_calibration_transport(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      nuisance = nuisance,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  } else {
    calib_result <- .tr_calibration_joint(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      nuisance = nuisance,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  }

  # Compute naive estimate (source only, no weighting)
  naive_result <- .naive_calibration_source(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    source = source,
    treatment_level = treatment_level,
    smoother = smoother,
    n_bins = n_bins,
    span = span
  )

  list(
    predicted = calib_result$predicted,
    observed = calib_result$observed,
    predictions_raw = predictions,  # Raw predictions for histogram
    weights = calib_result$weights,
    smoother = smoother,
    ici = calib_result$ici,
    e50 = calib_result$e50,
    e90 = calib_result$e90,
    emax = calib_result$emax,
    naive_ici = naive_result$ici,
    n_source = n_source,
    n_target = n_target,
    selection_model = nuisance$selection_model,
    propensity_model = nuisance$propensity_model,
    outcome_model = nuisance$outcome_model
  )
}


#' Fit Nuisance Models for Calibration
#' @keywords internal
#' @noRd
.fit_calibration_nuisance <- function(outcomes, treatment, source, covariates,
                                      treatment_level, analysis, estimator,
                                      selection_model, propensity_model,
                                      outcome_model, predictions) {

  n <- length(outcomes)
  cov_df <- as.data.frame(covariates)
  
  # Determine if factual mode (no treatment)
  factual_mode <- is.null(treatment)

  # Selection model: P(S=0|X)
  if (is.null(selection_model)) {
    sel_data <- cbind(S0 = 1 - source, cov_df)
    selection_model <- glm(S0 ~ ., data = sel_data, family = binomial())
  }

  # Get selection probabilities
  p_s0 <- predict(selection_model, newdata = cov_df, type = "response")
  p_s1 <- 1 - p_s0

  # Propensity model - only needed for counterfactual mode
  if (factual_mode) {
    propensity_model <- NULL
    ps <- rep(1, n)  # No propensity needed
  } else {
    if (is.null(propensity_model)) {
      if (analysis == "transport") {
        # P(A=1|X, S=1) using source data only
        source_idx <- source == 1
        ps_data <- cbind(A = treatment[source_idx], cov_df[source_idx, , drop = FALSE])
        propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
      } else {
        # P(A=1|X) using all data
        ps_data <- cbind(A = treatment, cov_df)
        propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
      }
    }
    
    # Get propensity scores
    ps <- predict(propensity_model, newdata = cov_df, type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
  }

  # Outcome model for OM and DR estimators
  if (estimator %in% c("om", "dr") && is.null(outcome_model)) {
    if (factual_mode) {
      # Factual mode: E[Y|X, S=1] using all source data
      if (analysis == "transport") {
        source_idx <- source == 1
      } else {
        source_idx <- rep(TRUE, n)
      }
      if (sum(source_idx) > 5) {
        om_data <- cbind(Y = outcomes[source_idx],
                         pred = predictions[source_idx],
                         cov_df[source_idx, , drop = FALSE])
        outcome_model <- glm(Y ~ pred, data = om_data, family = binomial())
      }
    } else {
      # Counterfactual mode: E[Y|X, A=a, S=1]
      if (analysis == "transport") {
        # E[Y|X, A=a, S=1] using source data
        source_a_idx <- source == 1 & treatment == treatment_level
        if (sum(source_a_idx) > 5) {
          om_data <- cbind(Y = outcomes[source_a_idx],
                           pred = predictions[source_a_idx],
                           cov_df[source_a_idx, , drop = FALSE])
          outcome_model <- glm(Y ~ pred, data = om_data, family = binomial())
        }
      } else {
        # E[Y|X, A=a] using all data with A=a
        a_idx <- treatment == treatment_level
        if (sum(a_idx) > 5) {
          om_data <- cbind(Y = outcomes[a_idx],
                           pred = predictions[a_idx],
                           cov_df[a_idx, , drop = FALSE])
          outcome_model <- glm(Y ~ pred, data = om_data, family = binomial())
        }
      }
    }
  }

  list(
    selection_model = selection_model,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    p_s0 = p_s0,
    p_s1 = p_s1,
    ps = ps
  )
}


#' Transportability Calibration - Transport Analysis
#' @keywords internal
#' @noRd
.tr_calibration_transport <- function(predictions, outcomes, treatment, source,
                                      covariates, treatment_level, estimator,
                                      nuisance, smoother, n_bins, span) {

  n <- length(outcomes)
  n_target <- sum(source == 0)

  # Determine if factual mode (no treatment)
  factual_mode <- is.null(treatment)
  
  # Indicator for treatment level in source (or all source in factual mode)
  if (factual_mode) {
    I_a <- rep(TRUE, n)
  } else {
    I_a <- treatment == treatment_level
  }
  I_source <- source == 1
  I_target <- source == 0

  if (estimator == "ipw") {
    # IPW weights: P(S=0|X) / [P(S=1|X) * P(A=a|X,S=1)]
    # For factual mode: P(S=0|X) / P(S=1|X)
    # Applied to source observations (with A=a for counterfactual)
    weights <- rep(0, n)
    idx <- I_source & I_a

    # Avoid division by very small numbers
    if (factual_mode) {
      denom <- nuisance$p_s1[idx]
    } else {
      denom <- nuisance$p_s1[idx] * nuisance$ps[idx]
    }
    denom <- pmax(denom, 1e-10)

    weights[idx] <- nuisance$p_s0[idx] / denom

    # Subset to source with treatment level
    pred_sub <- predictions[idx]
    out_sub <- outcomes[idx]
    w_sub <- weights[idx]

    # Normalize weights to sum to n_target
    w_sub <- w_sub / sum(w_sub) * n_target

    # Compute weighted calibration curve
    calib <- .compute_weighted_calibration(
      predictions = pred_sub,
      outcomes = out_sub,
      weights = w_sub,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )

  } else if (estimator == "om") {
    # OM estimator: use outcome model predictions
    if (is.null(nuisance$outcome_model)) {
      stop("Outcome model fitting failed - insufficient data")
    }

    # Predict E[Y|X, pred] for target population
    target_df <- data.frame(pred = predictions[I_target])
    h_a <- predict(nuisance$outcome_model, newdata = target_df, type = "response")

    # Calibration based on model predictions
    pred_sub <- predictions[I_target]
    out_sub <- h_a
    w_sub <- rep(1, n_target)

    calib <- .compute_weighted_calibration(
      predictions = pred_sub,
      outcomes = out_sub,
      weights = w_sub,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  } else {
    # DR estimator: augmented IPW
    if (is.null(nuisance$outcome_model)) {
      stop("Outcome model fitting failed - insufficient data for DR estimator")
    }

    # Get outcome model predictions for all observations
    all_df <- data.frame(pred = predictions)
    mu_hat <- predict(nuisance$outcome_model, newdata = all_df, type = "response")

    # IPW weights for source observations (with A=a for counterfactual)
    idx <- I_source & I_a
    if (factual_mode) {
      denom <- nuisance$p_s1[idx]
    } else {
      denom <- nuisance$p_s1[idx] * nuisance$ps[idx]
    }
    denom <- pmax(denom, 1e-10)
    ipw_weight <- nuisance$p_s0[idx] / denom

    # Augmented pseudo-outcomes for source observations
    augmentation <- ipw_weight * (outcomes[idx] - mu_hat[idx])

    # Pseudo-outcomes: mu_hat for target + augmentation for source
    # We compute over ALL target observations using mu_hat, plus augmentation from source
    pred_all <- c(predictions[I_target], predictions[idx])
    pseudo_out <- c(mu_hat[I_target], mu_hat[idx] + augmentation)
    w_all <- c(rep(1, n_target), rep(1, sum(idx)))

    # Normalize weights
    w_all <- w_all / sum(w_all) * n_target

    calib <- .compute_weighted_calibration(
      predictions = pred_all,
      outcomes = pseudo_out,
      weights = w_all,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  }

  calib
}


#' Transportability Calibration - Joint Analysis
#' @keywords internal
#' @noRd
.tr_calibration_joint <- function(predictions, outcomes, treatment, source,
                                  covariates, treatment_level, estimator,
                                  nuisance, smoother, n_bins, span) {

  n <- length(outcomes)
  n_target <- sum(source == 0)

  # Determine if factual mode (no treatment)
  factual_mode <- is.null(treatment)

  # Indicator for treatment level (or all observations in factual mode)
  if (factual_mode) {
    I_a <- rep(TRUE, n)
  } else {
    I_a <- treatment == treatment_level
  }
  I_target <- source == 0

  if (estimator == "ipw") {
    # IPW weights: P(S=0|X) / P(A=a|X)
    # For factual mode: Just use P(S=0|X) as weight
    # Applied to all observations (with A=a for counterfactual)
    weights <- rep(0, n)
    idx <- I_a

    if (factual_mode) {
      # In factual mode for joint, we use all data weighted by P(S=0|X)
      weights[idx] <- nuisance$p_s0[idx]
    } else {
      # Avoid division by very small numbers
      denom <- nuisance$ps[idx]
      denom <- pmax(denom, 1e-10)
      weights[idx] <- nuisance$p_s0[idx] / denom
    }

    # Subset to treatment level
    pred_sub <- predictions[idx]
    out_sub <- outcomes[idx]
    w_sub <- weights[idx]

    # Normalize weights to sum to n_target
    w_sub <- w_sub / sum(w_sub) * n_target

    # Compute weighted calibration curve
    calib <- .compute_weighted_calibration(
      predictions = pred_sub,
      outcomes = out_sub,
      weights = w_sub,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )

  } else if (estimator == "om") {
    # OM estimator: use outcome model predictions for target
    if (is.null(nuisance$outcome_model)) {
      stop("Outcome model fitting failed - insufficient data")
    }

    # Predict for target population
    target_df <- data.frame(pred = predictions[I_target])
    h_a <- predict(nuisance$outcome_model, newdata = target_df, type = "response")

    pred_sub <- predictions[I_target]
    out_sub <- h_a
    w_sub <- rep(1, n_target)

    calib <- .compute_weighted_calibration(
      predictions = pred_sub,
      outcomes = out_sub,
      weights = w_sub,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  } else {
    # DR estimator: augmented IPW for joint analysis
    if (is.null(nuisance$outcome_model)) {
      stop("Outcome model fitting failed - insufficient data for DR estimator")
    }

    # Get outcome model predictions for all observations
    all_df <- data.frame(pred = predictions)
    mu_hat <- predict(nuisance$outcome_model, newdata = all_df, type = "response")

    # IPW weights: P(S=0|X) / P(A=a|X) for all observations with A=a
    # For factual mode: just P(S=0|X) as weight
    idx <- I_a
    if (factual_mode) {
      ipw_weight <- nuisance$p_s0[idx]
    } else {
      denom <- nuisance$ps[idx]
      denom <- pmax(denom, 1e-10)
      ipw_weight <- nuisance$p_s0[idx] / denom
    }

    # Augmented pseudo-outcomes
    augmentation <- ipw_weight * (outcomes[idx] - mu_hat[idx])

    # Combine: mu_hat for target + augmented for those with A=a
    pred_all <- c(predictions[I_target], predictions[idx])
    pseudo_out <- c(mu_hat[I_target], mu_hat[idx] + augmentation)
    w_all <- c(rep(1, n_target), rep(1, sum(idx)))

    # Normalize weights
    w_all <- w_all / sum(w_all) * n_target

    calib <- .compute_weighted_calibration(
      predictions = pred_all,
      outcomes = pseudo_out,
      weights = w_all,
      smoother = smoother,
      n_bins = n_bins,
      span = span
    )
  }

  calib
}


#' Compute Weighted Calibration Curve
#' @keywords internal
#' @noRd
.compute_weighted_calibration <- function(predictions, outcomes, weights,
                                          smoother, n_bins, span) {

  if (smoother == "loess") {
    # LOESS smoothing
    fit <- stats::loess(outcomes ~ predictions, weights = weights, span = span)
    # Get unique sorted predictions for curve
    pred_unique <- sort(unique(predictions))
    observed <- predict(fit, newdata = pred_unique)

    # Handle NAs at edges
    valid <- !is.na(observed)
    predicted <- pred_unique[valid]
    observed <- observed[valid]

  } else {
    # Binned calibration
    breaks <- seq(0, 1, length.out = n_bins + 1)
    bins <- cut(predictions, breaks = breaks, include.lowest = TRUE)

    # Weighted means within bins
    predicted <- tapply(predictions * weights, bins, sum) /
                 tapply(weights, bins, sum)
    observed <- tapply(outcomes * weights, bins, sum) /
                tapply(weights, bins, sum)

    # Remove empty bins
    valid <- !is.na(predicted) & !is.na(observed)
    predicted <- predicted[valid]
    observed <- observed[valid]
  }

  # Compute calibration metrics
  abs_diff <- abs(observed - predicted)
  ici <- mean(abs_diff, na.rm = TRUE)
  e50 <- stats::median(abs_diff, na.rm = TRUE)
  e90 <- stats::quantile(abs_diff, 0.9, na.rm = TRUE, names = FALSE)
  emax <- max(abs_diff, na.rm = TRUE)

  list(
    predicted = as.numeric(predicted),
    observed = as.numeric(observed),
    weights = weights,
    ici = ici,
    e50 = e50,
    e90 = e90,
    emax = emax
  )
}


#' Naive Calibration from Source Only
#' @keywords internal
#' @noRd
.naive_calibration_source <- function(predictions, outcomes, treatment, source,
                                      treatment_level, smoother, n_bins, span) {

  # Determine if factual mode (no treatment)
  factual_mode <- is.null(treatment)
  
  # Source observations with treatment level (no weighting)
  # For factual mode, use all source data
  if (factual_mode) {
    idx <- source == 1
  } else {
    idx <- source == 1 & treatment == treatment_level
  }

  pred_sub <- predictions[idx]
  out_sub <- outcomes[idx]
  w_sub <- rep(1, sum(idx))

  .compute_weighted_calibration(
    predictions = pred_sub,
    outcomes = out_sub,
    weights = w_sub,
    smoother = smoother,
    n_bins = n_bins,
    span = span
  )
}


#' Bootstrap Standard Errors for Transportable Calibration
#' @keywords internal
#' @noRd
.bootstrap_tr_calibration <- function(predictions, outcomes, treatment, source,
                                      covariates, treatment_level, analysis,
                                      estimator, smoother, n_bins, span,
                                      n_boot, stratified, conf_level = 0.95,
                                      eval_points = NULL, ...) {

  n <- length(outcomes)

  # Function to generate bootstrap indices
  get_boot_idx <- function() {
    if (stratified) {
      # Stratified bootstrap: sample separately within source and target
      source_idx <- which(source == 1)
      target_idx <- which(source == 0)
      c(sample(source_idx, replace = TRUE),
        sample(target_idx, replace = TRUE))
    } else {
      sample(n, replace = TRUE)
    }
  }

  # Bootstrap samples for metrics
  boot_ici <- numeric(n_boot)
  boot_e50 <- numeric(n_boot)
  boot_e90 <- numeric(n_boot)
  boot_emax <- numeric(n_boot)

  # Storage for bootstrap curves (evaluated at common points)
  n_eval <- if (!is.null(eval_points)) length(eval_points) else 0
  boot_observed <- if (n_eval > 0) matrix(NA, nrow = n_boot, ncol = n_eval) else NULL

  for (b in seq_len(n_boot)) {
    boot_idx <- get_boot_idx()

    tryCatch({
      # Handle NULL treatment for factual mode
      boot_treatment <- if (is.null(treatment)) NULL else treatment[boot_idx]
      
      boot_result <- .compute_tr_calibration(
        predictions = predictions[boot_idx],
        outcomes = outcomes[boot_idx],
        treatment = boot_treatment,
        source = source[boot_idx],
        covariates = covariates[boot_idx, , drop = FALSE],
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = estimator,
        selection_model = NULL,
        propensity_model = NULL,
        outcome_model = NULL,
        smoother = smoother,
        n_bins = n_bins,
        span = span,
        ...
      )

      boot_ici[b] <- boot_result$ici
      boot_e50[b] <- boot_result$e50
      boot_e90[b] <- boot_result$e90
      boot_emax[b] <- boot_result$emax

      # Interpolate bootstrap curve to common evaluation points
      if (n_eval > 0 && length(boot_result$predicted) > 1) {
        boot_observed[b, ] <- approx(
          x = boot_result$predicted,
          y = boot_result$observed,
          xout = eval_points,
          rule = 2  # Use nearest value for extrapolation
        )$y
      }
    }, error = function(e) {
      boot_ici[b] <<- NA
      boot_e50[b] <<- NA
      boot_e90[b] <<- NA
      boot_emax[b] <<- NA
    })
  }

  # Remove failed bootstraps
  valid <- !is.na(boot_ici)

  # Compute SEs
  se_ici <- sd(boot_ici[valid], na.rm = TRUE)
  se_e50 <- sd(boot_e50[valid], na.rm = TRUE)
  se_e90 <- sd(boot_e90[valid], na.rm = TRUE)
  se_emax <- sd(boot_emax[valid], na.rm = TRUE)

  # Confidence intervals (percentile method)
  alpha <- 1 - conf_level
  ci_ici <- quantile(boot_ici[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e50 <- quantile(boot_e50[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e90 <- quantile(boot_e90[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_emax <- quantile(boot_emax[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  # Compute pointwise CIs for the calibration curve
  curves <- NULL
  if (n_eval > 0 && !is.null(boot_observed)) {
    ci_lower_curve <- apply(boot_observed[valid, , drop = FALSE], 2,
                            quantile, probs = alpha/2, na.rm = TRUE)
    ci_upper_curve <- apply(boot_observed[valid, , drop = FALSE], 2,
                            quantile, probs = 1 - alpha/2, na.rm = TRUE)
    curves <- data.frame(
      predicted = eval_points,
      ci_lower = ci_lower_curve,
      ci_upper = ci_upper_curve
    )
  }

  list(
    se_ici = se_ici,
    se_e50 = se_e50,
    se_e90 = se_e90,
    se_emax = se_emax,
    ci_ici = as.numeric(ci_ici),
    ci_e50 = as.numeric(ci_e50),
    ci_e90 = as.numeric(ci_e90),
    ci_emax = as.numeric(ci_emax),
    curves = curves,
    samples = data.frame(
      ici = boot_ici[valid],
      e50 = boot_e50[valid],
      e90 = boot_e90[valid],
      emax = boot_emax[valid]
    )
  )
}


#' Plot Method for tr_calibration Objects
#'
#' Creates a calibration plot showing predicted vs observed probabilities
#' in the target population, with optional confidence bands and histogram.
#'
#' @param x A `tr_calibration` object.
#' @param add_reference Logical; add 45-degree reference line (default: TRUE).
#' @param show_metrics Logical; show calibration metrics on plot (default: TRUE).
#' @param add_histogram Logical; add histogram of predictions below the
#'   calibration curve (default: TRUE).
#' @param add_rug Logical; add rug plot to show individual predictions
#'   (default: FALSE).
#' @param add_ci Logical; add bootstrap confidence bands if available
#'   (default: TRUE if bootstrap was run).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object (if ggplot2 available) or base R plot.
#'
#' @export
plot.tr_calibration <- function(x, add_reference = TRUE, show_metrics = TRUE, 
                                 add_histogram = TRUE, add_rug = FALSE,
                                 add_ci = !is.null(x$boot_curves), ...) {

  estimator_label <- toupper(x$estimator)
  has_ci <- !is.null(x$boot_curves) && add_ci
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R plot - simple version without histogram
    if (add_histogram && !is.null(x$predictions_raw)) {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
      par(mar = c(0, 4, 3, 2))
    }
    
    plot(x$predicted, x$observed,
         type = "l", lwd = 2, col = "#2E86AB",
         xlab = if (add_histogram) "" else "Predicted probability",
         ylab = sprintf("Probability in Target (%s)", estimator_label),
         main = "Calibration Curve in the Target Population",
         xlim = c(0, 1), ylim = c(0, 1),
         xaxt = if (add_histogram) "n" else "s")
    if (add_reference) {
      abline(0, 1, lty = 2, col = "gray")
    }
    
    # Add CI band if available
    if (has_ci) {
      polygon(c(x$boot_curves$predicted, rev(x$boot_curves$predicted)),
              c(x$boot_curves$ci_lower, rev(x$boot_curves$ci_upper)),
              col = rgb(46/255, 134/255, 171/255, 0.2), border = NA)
    }
    
    if (add_histogram && !is.null(x$predictions_raw)) {
      par(mar = c(4, 4, 0, 2))
      hist(x$predictions_raw, breaks = 30, col = "#2E86AB", border = "white",
           main = "", xlab = "Predicted probability", xlim = c(0, 1))
    }
    
    return(invisible(NULL))
  }

  # ggplot2 version
  df <- data.frame(
    predicted = x$predicted,
    observed = x$observed
  )

  subtitle_text <- NULL
  if (show_metrics) {
    subtitle_text <- sprintf("ICI = %.3f, E50 = %.3f, Emax = %.3f%s",
                             x$ici, x$e50, x$emax,
                             if (has_ci) sprintf(" (%d%% CI)", round(x$conf_level * 100)) else "")
  }

  p_cal <- ggplot2::ggplot(df, ggplot2::aes(x = .data$predicted, y = .data$observed))
  
  # Add confidence band first (so it's behind the line)
  if (has_ci) {
    df_ci <- x$boot_curves
    p_cal <- p_cal +
      ggplot2::geom_ribbon(
        data = df_ci,
        ggplot2::aes(x = .data$predicted, ymin = .data$ci_lower, ymax = .data$ci_upper),
        fill = "#2E86AB", alpha = 0.2, inherit.aes = FALSE
      )
  }
  
  p_cal <- p_cal +
    ggplot2::geom_line(linewidth = 1.2, color = "#2E86AB") +
    ggplot2::labs(
      y = sprintf("Probability in Target (%s)", estimator_label),
      title = "Calibration Curve in the Target Population",
      subtitle = subtitle_text
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw()

  if (add_reference) {
    p_cal <- p_cal + ggplot2::geom_abline(slope = 1, intercept = 0,
                                           linetype = "dashed", color = "gray50")
  }

  # Add rug if requested
  if (add_rug && !is.null(x$predictions_raw)) {
    df_raw <- data.frame(predictions_raw = x$predictions_raw)
    p_cal <- p_cal + 
      ggplot2::geom_rug(data = df_raw, 
                        ggplot2::aes(x = .data$predictions_raw, y = NULL),
                        alpha = 0.3, color = "#2E86AB")
  }

  # Add histogram subplot if requested
  if (add_histogram && !is.null(x$predictions_raw)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      # Fall back to just the calibration plot with rug
      message("Install 'patchwork' package for histogram subplot. Showing rug plot instead.")
      df_raw <- data.frame(predictions_raw = x$predictions_raw)
      p_cal <- p_cal + 
        ggplot2::geom_rug(data = df_raw, 
                          ggplot2::aes(x = .data$predictions_raw, y = NULL),
                          alpha = 0.3, color = "#2E86AB") +
        ggplot2::labs(x = "Predicted probability")
      return(p_cal)
    }
    
    # Remove x-axis label from calibration plot
    p_cal <- p_cal + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    # Create histogram
    df_raw <- data.frame(predictions_raw = x$predictions_raw)
    p_hist <- ggplot2::ggplot(df_raw, ggplot2::aes(x = .data$predictions_raw)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2E86AB", color = "white", alpha = 0.8) +
      ggplot2::labs(x = "Predicted probability", y = "Count") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 5.5, l = 5.5))
    
    # Combine plots
    p_combined <- patchwork::wrap_plots(p_cal, p_hist, ncol = 1, heights = c(3, 1))
    return(p_combined)
  }
  
  # No histogram - add x-axis label
  p_cal <- p_cal + ggplot2::labs(x = "Predicted probability")
  return(p_cal)
}
