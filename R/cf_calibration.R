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
#'
#' @return An object of class `c("cf_calibration", "cf_performance")` containing:
#'   \item{predicted}{Vector of predicted probabilities}
#'   \item{observed}{Vector of smoothed observed probabilities}
#'   \item{smoother}{Smoothing method used}
#'   \item{ici}{Integrated calibration index}
#'   \item{e50}{Median absolute calibration error}
#'   \item{e90}{90th percentile absolute calibration error}
#'   \item{emax}{Maximum absolute calibration error}
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
#' E[Y|X, A=a] to estimate calibration over all observations. Requires a
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
                           ...) {

  estimator <- match.arg(estimator)
  smoother <- match.arg(smoother)

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
    # Truncate extreme propensities for stability
    ps <- pmax(pmin(ps, 0.975), 0.025)
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
    call = match.call()
  )

  class(result) <- c("cf_calibration", "cf_performance")
  return(result)
}
