#' Estimate Counterfactual Calibration Curve
#'
#' Estimates the calibration curve of a prediction model under a hypothetical
#' intervention where treatment is set to a specific level.
#'
#' @inheritParams cf_mse
#' @param smoother Smoothing method for the calibration curve:
#'   - `"loess"`: Local polynomial regression (default)
#'   - `"binned"`: Binned calibration
#' @param n_bins Number of bins for binned calibration (default: 10).
#' @param span Span parameter for LOESS smoothing (default: 0.75).
#' @param se_method Method for standard error estimation:
#'   - `"none"`: No standard errors (default, fastest)
#'   - `"bootstrap"`: Bootstrap standard errors
#' @param n_boot Number of bootstrap replications (default: 200).
#' @param conf_level Confidence level for intervals (default: 0.95).
#'
#' @return An object of class `c("cf_calibration", "cf_performance")` containing:
#'   \item{predicted}{Vector of predicted probabilities}
#'   \item{observed}{Vector of smoothed observed probabilities}
#'   \item{smoother}{Smoothing method used}
#'   \item{ici}{Integrated calibration index}
#'   \item{e50}{Median absolute calibration error}
#'   \item{e90}{90th percentile absolute calibration error}
#'   \item{emax}{Maximum absolute calibration error}
#'   \item{se}{List of standard errors (if se_method = "bootstrap")}
#'   \item{ci_lower}{List of lower CI bounds (if se_method = "bootstrap")}
#'   \item{ci_upper}{List of upper CI bounds (if se_method = "bootstrap")}
#'   \item{boot_curves}{Bootstrap calibration curves for CI bands (if se_method = "bootstrap")}
#'
#' @details
#' The counterfactual calibration curve estimates the relationship between
#' predicted risk and observed risk under the counterfactual intervention.
#'
#' The function implements three estimators:
#'
#' **IPW Estimator**: Weights observations by the inverse probability of
#' receiving the counterfactual treatment. Requires a correctly specified
#' propensity score model.
#'
#' **Conditional Loss (CL) Estimator**: Uses the fitted outcome model
#' \eqn{E[Y | X, A=a]} to estimate calibration over all observations. Requires a
#' correctly specified outcome model.
#'
#' **Doubly Robust (DR) Estimator**: Combines CL and IPW approaches. Consistent
#' if either the propensity or outcome model is correctly specified.
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' Steingrimsson, J. A., Gatsonis, C., Li, B., & Dahabreh, I. J. (2023).
#' "Transporting a prediction model for use in a new target population."
#' *American Journal of Epidemiology*, 192(2), 296-304.
#'
#' @seealso [cf_mse()], [cf_auc()], [plot.cf_calibration()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate counterfactual calibration curve with different estimators
#' result_ipw <- cf_calibration(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "ipw"
#' )
#'
#' result_dr <- cf_calibration(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr"
#' )
#' print(result_dr)
#' # plot(result_dr)  # If ggplot2 is available
#'
#' # With bootstrap confidence bands
#' # result_boot <- cf_calibration(
#' #   predictions = pred, outcomes = y, treatment = a,
#' #   covariates = data.frame(x = x), treatment_level = 0,
#' #   se_method = "bootstrap", n_boot = 200
#' # )
#' # plot(result_boot)  # Shows confidence bands
cf_calibration <- function(predictions,
                           outcomes,
                           treatment,
                           covariates,
                           treatment_level = 0,
                           estimator = c("dr", "ipw", "cl"),
                           propensity_model = NULL,
                           outcome_model = NULL,
                           smoother = c("loess", "binned"),
                           n_bins = 10,
                           span = 0.75,
                           se_method = c("none", "bootstrap"),
                           n_boot = 200,
                           conf_level = 0.95,
                           ps_trim = NULL,
                           ...) {

  estimator <- match.arg(estimator)
  smoother <- match.arg(smoother)
  se_method <- match.arg(se_method)

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  # Validate inputs
  .validate_inputs(predictions, outcomes, treatment, covariates)

  if (!all(outcomes %in% c(0, 1))) {
    stop("Calibration requires binary outcomes (0/1)")
  }

  n <- length(outcomes)

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Fit propensity model if needed for IPW or DR
  if (estimator %in% c("ipw", "dr") && is.null(propensity_model)) {
    ps_data <- cbind(A = treatment, covariates)
    propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
  }

  # Fit outcome model if needed for CL or DR
  if (estimator %in% c("cl", "dr") && is.null(outcome_model)) {
    subset_idx <- treatment == treatment_level
    outcome_data <- cbind(Y = outcomes, covariates)[subset_idx, ]
    outcome_model <- glm(Y ~ ., data = outcome_data, family = binomial())
  }

  # Get propensity scores
  if (estimator %in% c("ipw", "dr")) {
    ps <- .predict_nuisance(propensity_model, covariates, type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
    # Trim extreme propensities for stability
    ps <- .trim_propensity(ps, ps_trim_spec$method, ps_trim_spec$bounds)
  }

  # Get outcome model predictions E[Y|X, A=a]
  if (estimator %in% c("cl", "dr")) {
    mu_hat <- .predict_nuisance(outcome_model, covariates, type = "response")
  }

  # Indicator for treatment level
  I_a <- as.numeric(treatment == treatment_level)

  # Compute pseudo-outcomes based on estimator
  if (estimator == "ipw") {
    # IPW: weight observations in treatment group
    # Use only observations with A = a
    pred_use <- predictions[I_a == 1]
    pseudo_outcomes <- outcomes[I_a == 1]
    weights <- 1 / ps[I_a == 1]
    # Normalize weights
    weights <- weights / mean(weights)

  } else if (estimator == "cl") {
    # CL: use outcome model predictions for all observations
    pred_use <- predictions
    pseudo_outcomes <- mu_hat
    weights <- rep(1, n)

  } else if (estimator == "dr") {
    # DR: augmented IPW over all observations
    # Pseudo-outcome: mu_hat + I(A=a)/ps * (Y - mu_hat)
    pred_use <- predictions
    augmentation <- I_a / ps * (outcomes - mu_hat)
    pseudo_outcomes <- mu_hat + augmentation
    weights <- rep(1, n)
  }

  # Compute calibration curve
  if (smoother == "loess") {
    fit <- loess(pseudo_outcomes ~ pred_use, weights = weights, span = span)
    predicted <- sort(unique(pred_use))
    observed <- predict(fit, newdata = predicted)
  } else if (smoother == "binned") {
    bins <- cut(pred_use, breaks = n_bins, include.lowest = TRUE)
    predicted <- tapply(pred_use, bins, mean)
    observed <- tapply(pseudo_outcomes * weights, bins, sum) / tapply(weights, bins, sum)
  }

  # Compute calibration metrics
  abs_diff <- abs(observed - predicted)
  ici <- mean(abs_diff, na.rm = TRUE)
  e50 <- median(abs_diff, na.rm = TRUE)
  e90 <- quantile(abs_diff, 0.9, na.rm = TRUE)
  emax <- max(abs_diff, na.rm = TRUE)

  result <- list(
    predicted = predicted,
    observed = observed,
    predictions_raw = pred_use,  # Raw predictions for histogram
    weights = weights,
    smoother = smoother,
    estimator = estimator,
    metric = "calibration",
    treatment_level = treatment_level,
    n_obs = if (estimator == "ipw") sum(I_a) else n,
    ici = ici,
    e50 = e50,
    e90 = e90,
    emax = emax,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    se_method = se_method,
    conf_level = conf_level,
    call = match.call()
  )

  # Add bootstrap standard errors and confidence bands if requested
  if (se_method == "bootstrap") {
    boot_results <- .bootstrap_cf_calibration(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      smoother = smoother,
      n_bins = n_bins,
      span = span,
      n_boot = n_boot,
      conf_level = conf_level,
      eval_points = predicted,  # Evaluate at same points as main estimate
      ...
    )

    # Standard errors for metrics
    result$se <- list(
      ici = boot_results$se_ici,
      e50 = boot_results$se_e50,
      e90 = boot_results$se_e90,
      emax = boot_results$se_emax
    )

    # Confidence intervals for metrics (percentile method)
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

    # Store bootstrap curves for plotting confidence bands
    result$boot_curves <- boot_results$curves
    result$boot_samples <- boot_results$samples
  } else {
    result$se <- NULL
    result$ci_lower <- NULL
    result$ci_upper <- NULL
    result$boot_curves <- NULL
    result$boot_samples <- NULL
  }

  class(result) <- c("cf_calibration", "cf_performance")
  return(result)
}


#' Bootstrap Standard Errors for Counterfactual Calibration
#' @keywords internal
#' @noRd
.bootstrap_cf_calibration <- function(predictions, outcomes, treatment,
                                       covariates, treatment_level, estimator,
                                       smoother, n_bins, span, n_boot,
                                       conf_level, eval_points, ...) {

  n <- length(outcomes)

  # Storage for bootstrap results
  boot_ici <- numeric(n_boot)
  boot_e50 <- numeric(n_boot)
  boot_e90 <- numeric(n_boot)
  boot_emax <- numeric(n_boot)

  # Storage for bootstrap curves (evaluated at common points)
  n_eval <- length(eval_points)
  boot_observed <- matrix(NA, nrow = n_boot, ncol = n_eval)

  for (b in seq_len(n_boot)) {
    boot_idx <- sample(n, replace = TRUE)

    tryCatch({
      # Refit on bootstrap sample
      boot_result <- cf_calibration(
        predictions = predictions[boot_idx],
        outcomes = outcomes[boot_idx],
        treatment = treatment[boot_idx],
        covariates = covariates[boot_idx, , drop = FALSE],
        treatment_level = treatment_level,
        estimator = estimator,
        propensity_model = NULL,  # Refit on bootstrap
        outcome_model = NULL,
        smoother = smoother,
        n_bins = n_bins,
        span = span,
        se_method = "none",  # Don't recurse
        ...
      )

      boot_ici[b] <- boot_result$ici
      boot_e50[b] <- boot_result$e50
      boot_e90[b] <- boot_result$e90
      boot_emax[b] <- boot_result$emax

      # Interpolate bootstrap curve to common evaluation points
      if (length(boot_result$predicted) > 1) {
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

  # Compute SEs for metrics
  se_ici <- sd(boot_ici[valid], na.rm = TRUE)
  se_e50 <- sd(boot_e50[valid], na.rm = TRUE)
  se_e90 <- sd(boot_e90[valid], na.rm = TRUE)
  se_emax <- sd(boot_emax[valid], na.rm = TRUE)

  # Confidence intervals for metrics (percentile method)
  alpha <- 1 - conf_level
  ci_ici <- quantile(boot_ici[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e50 <- quantile(boot_e50[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_e90 <- quantile(boot_e90[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  ci_emax <- quantile(boot_emax[valid], c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  # Compute pointwise CIs for the calibration curve
  ci_lower_curve <- apply(boot_observed[valid, , drop = FALSE], 2,
                          quantile, probs = alpha/2, na.rm = TRUE)
  ci_upper_curve <- apply(boot_observed[valid, , drop = FALSE], 2,
                          quantile, probs = 1 - alpha/2, na.rm = TRUE)

  list(
    se_ici = se_ici,
    se_e50 = se_e50,
    se_e90 = se_e90,
    se_emax = se_emax,
    ci_ici = as.numeric(ci_ici),
    ci_e50 = as.numeric(ci_e50),
    ci_e90 = as.numeric(ci_e90),
    ci_emax = as.numeric(ci_emax),
    curves = data.frame(
      predicted = eval_points,
      ci_lower = ci_lower_curve,
      ci_upper = ci_upper_curve
    ),
    samples = data.frame(
      ici = boot_ici[valid],
      e50 = boot_e50[valid],
      e90 = boot_e90[valid],
      emax = boot_emax[valid]
    )
  )
}
