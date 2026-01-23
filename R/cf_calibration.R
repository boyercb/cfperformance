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
#' This is done by applying inverse probability weights to the calibration
#' curve estimation.
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
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
#' # Estimate counterfactual calibration curve
#' result <- cf_calibration(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0
#' )
#' print(result)
#' # plot(result)  # If ggplot2 is available
cf_calibration <- function(predictions,
                           outcomes,
                           treatment,
                           covariates,
                           treatment_level = 0,
                           estimator = c("ipw", "cl"),
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

  # Fit propensity model if needed for IPW
  if (estimator == "ipw" && is.null(propensity_model)) {
    ps_data <- cbind(A = treatment, as.data.frame(covariates))
    propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
  }

  # Get propensity scores
  if (estimator == "ipw") {
    ps <- .predict_nuisance(propensity_model, as.data.frame(covariates), type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
  }

  # Indicator for treatment level
  I_a <- treatment == treatment_level

  # Compute weights for IPW calibration
  if (estimator == "ipw") {
    weights <- rep(0, n)
    weights[I_a] <- 1 / ps[I_a]
    # Normalize weights
    weights <- weights / sum(weights) * sum(I_a)
  } else {
    weights <- rep(1, n)
  }

  # Subset to counterfactual treatment group
  pred_sub <- predictions[I_a]
  out_sub <- outcomes[I_a]
  w_sub <- weights[I_a]

  # Compute calibration curve
  if (smoother == "loess") {
    fit <- loess(out_sub ~ pred_sub, weights = w_sub, span = span)
    predicted <- sort(unique(pred_sub))
    observed <- predict(fit, newdata = predicted)
  } else if (smoother == "binned") {
    bins <- cut(pred_sub, breaks = n_bins, include.lowest = TRUE)
    predicted <- tapply(pred_sub, bins, mean)
    observed <- tapply(out_sub * w_sub, bins, sum) / tapply(w_sub, bins, sum)
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
    weights = w_sub,
    smoother = smoother,
    estimator = estimator,
    metric = "calibration",
    treatment_level = treatment_level,
    n_obs = sum(I_a),
    ici = ici,
    e50 = e50,
    e90 = e90,
    emax = emax,
    propensity_model = propensity_model,
    call = match.call()
  )

  class(result) <- c("cf_calibration", "cf_performance")
  return(result)
}
