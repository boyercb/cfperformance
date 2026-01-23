#' Estimate (Counterfactual) Calibration in the Target Population
#'
#' Estimates the calibration of a prediction model in a target population
#' using data transported from a source population (typically an RCT).
#'
#' @param predictions Numeric vector of model predictions (typically probabilities).
#' @param outcomes Numeric vector of observed outcomes (must be binary 0/1).
#' @param treatment Numeric vector of treatment indicators (0/1).
#' @param source Numeric vector of population indicators (1=source/RCT, 0=target).
#' @param covariates A matrix or data frame of baseline covariates.
#' @param treatment_level The treatment level of interest (default: 0).
#' @param analysis Character string specifying the type of analysis:
#'   - `"transport"`: Use source outcomes for target estimation (default)
#'   - `"joint"`: Pool source and target data
#' @param estimator Character string specifying the estimator:
#'   - `"ipw"`: Inverse probability weighting estimator (default)
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
#'   - `"bootstrap"`: Bootstrap standard errors (default)
#'   - `"none"`: No standard error estimation
#' @param n_boot Number of bootstrap replications (default: 500).
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
#'
#' @details
#' This function implements estimators for transporting prediction model
#' calibration from a source population (typically an RCT) to a target
#' population.
#'
#' **Transportability Analysis**: Uses outcome data from the source/RCT
#' population to estimate calibration in the target population. The IPW
#' estimator weights source observations to represent the target population.
#'
#' **Joint Analysis**: Pools source and target data to estimate calibration
#' in the target population.
#'
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
                           treatment,
                           source,
                           covariates,
                           treatment_level = 0,
                           analysis = c("transport", "joint"),
                           estimator = c("ipw", "om"),
                           selection_model = NULL,
                           propensity_model = NULL,
                           outcome_model = NULL,
                           smoother = c("loess", "binned"),
                           n_bins = 10,
                           span = 0.75,
                           se_method = c("bootstrap", "none"),
                           n_boot = 500,
                           conf_level = 0.95,
                           stratified_boot = TRUE,
                           ...) {

  # Match arguments
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  smoother <- match.arg(smoother)
  se_method <- match.arg(se_method)

  # Validate inputs
  .validate_tr_calibration_inputs(predictions, outcomes, treatment, source, covariates)

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
    alpha <- 1 - conf_level
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

    result$boot_samples <- boot_results$samples
  } else {
    result$se <- NULL
    result$ci_lower <- NULL
    result$ci_upper <- NULL
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
                                            source, covariates) {
  n <- length(outcomes)

  # Check lengths match
  if (length(predictions) != n) {
    stop("'predictions' must have the same length as 'outcomes'")
  }
  if (length(treatment) != n) {
    stop("'treatment' must have the same length as 'outcomes'")
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

  # Check treatment indicator
  if (!all(treatment %in% c(0, 1))) {
    stop("'treatment' must be a binary indicator (0/1)")
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

  # Indicator for treatment level
  I_a <- treatment == treatment_level

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

  # Selection model: P(S=0|X)
  if (is.null(selection_model)) {
    sel_data <- cbind(S0 = 1 - source, cov_df)
    selection_model <- glm(S0 ~ ., data = sel_data, family = binomial())
  }

  # Get selection probabilities
  p_s0 <- predict(selection_model, newdata = cov_df, type = "response")
  p_s1 <- 1 - p_s0

  # Propensity model
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

  # Outcome model for OM estimator
  if (estimator == "om" && is.null(outcome_model)) {
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

  # Indicator for treatment level in source
  I_a <- treatment == treatment_level
  I_source <- source == 1
  I_target <- source == 0

  if (estimator == "ipw") {
    # IPW weights: P(S=0|X) / [P(S=1|X) * P(A=a|X,S=1)]
    # Applied to source observations with A=a
    weights <- rep(0, n)
    idx <- I_source & I_a

    # Avoid division by very small numbers
    denom <- nuisance$p_s1[idx] * nuisance$ps[idx]
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

  } else {
    # OM estimator: use outcome model predictions
    if (is.null(nuisance$outcome_model)) {
      stop("Outcome model fitting failed - insufficient data")
    }

    # Predict E[Y|X, pred, A=a, S=1] for target population
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

  # Indicator for treatment level
  I_a <- treatment == treatment_level
  I_target <- source == 0

  if (estimator == "ipw") {
    # IPW weights: P(S=0|X) / P(A=a|X)
    # Applied to all observations with A=a
    weights <- rep(0, n)
    idx <- I_a

    # Avoid division by very small numbers
    denom <- nuisance$ps[idx]
    denom <- pmax(denom, 1e-10)

    weights[idx] <- nuisance$p_s0[idx] / denom

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

  } else {
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

  # Source observations with treatment level (no weighting)
  idx <- source == 1 & treatment == treatment_level

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
                                      n_boot, stratified, ...) {

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

  # Bootstrap samples
  boot_ici <- numeric(n_boot)
  boot_e50 <- numeric(n_boot)
  boot_e90 <- numeric(n_boot)
  boot_emax <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    boot_idx <- get_boot_idx()

    tryCatch({
      boot_result <- .compute_tr_calibration(
        predictions = predictions[boot_idx],
        outcomes = outcomes[boot_idx],
        treatment = treatment[boot_idx],
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
  alpha <- 0.05  # Will be overridden by conf_level in main function
  ci_ici <- quantile(boot_ici[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e50 <- quantile(boot_e50[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e90 <- quantile(boot_e90[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_emax <- quantile(boot_emax[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  list(
    se_ici = se_ici,
    se_e50 = se_e50,
    se_e90 = se_e90,
    se_emax = se_emax,
    ci_ici = as.numeric(ci_ici),
    ci_e50 = as.numeric(ci_e50),
    ci_e90 = as.numeric(ci_e90),
    ci_emax = as.numeric(ci_emax),
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
#' in the target population.
#'
#' @param x A `tr_calibration` object.
#' @param add_reference Logical; add 45-degree reference line (default: TRUE).
#' @param show_metrics Logical; show calibration metrics on plot (default: TRUE).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object (if ggplot2 available) or base R plot.
#'
#' @export
plot.tr_calibration <- function(x, add_reference = TRUE, show_metrics = TRUE, ...) {

  estimator_label <- toupper(x$estimator)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R plot
    plot(x$predicted, x$observed,
         type = "l", lwd = 2,
         xlab = "Predicted probability",
         ylab = sprintf("Probability in Target Population (%s)", estimator_label),
         main = "Calibration Curve in the Target Population",
         xlim = c(0, 1), ylim = c(0, 1))
    if (add_reference) {
      abline(0, 1, lty = 2, col = "gray")
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
    subtitle_text <- sprintf("ICI = %.3f, E50 = %.3f, Emax = %.3f",
                             x$ici, x$e50, x$emax)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_line(linewidth = 1.2, color = "#2E86AB") +
    ggplot2::labs(
      x = "Predicted probability",
      y = sprintf("Probability in Target Population (%s)", estimator_label),
      title = "Calibration Curve in the Target Population",
      subtitle = subtitle_text
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw()

  if (add_reference) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0,
                                   linetype = "dashed", color = "gray50")
  }

  return(p)
}
