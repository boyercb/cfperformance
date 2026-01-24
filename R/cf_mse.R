#' Estimate Counterfactual Mean Squared Error
#'
#' Estimates the mean squared error (MSE) of a prediction model under a
#' hypothetical intervention where treatment is set to a specific level.
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of observed outcomes.
#' @param treatment Numeric vector of treatment indicators (0/1).
#' @param covariates A matrix or data frame of baseline covariates (confounders).
#' @param treatment_level The counterfactual treatment level (default: 0).
#' @param estimator Character string specifying the estimator:
#'   - `"naive"`: Naive estimator (biased)
#'   - `"cl"`: Conditional loss estimator
#'   - `"ipw"`: Inverse probability weighting estimator
#'   - `"dr"`: Doubly robust estimator (default)
#' @param propensity_model Optional fitted propensity score model. If NULL,
#'   a logistic regression model is fit using the covariates.
#' @param outcome_model Optional fitted outcome model. If NULL, a regression
#'   model is fit using the covariates among treated/untreated. For binary
#'   outcomes, this should be a model for E\[Y|X,A\] (binomial family). For
#'   continuous outcomes, this should be a model for E\[L|X,A\] (gaussian family).
#' @param outcome_type Character string specifying the outcome type:
#'   - `"auto"`: Auto-detect from data (default)
#'   - `"binary"`: Binary outcome (0/1) - uses efficient transformation
#'   - `"continuous"`: Continuous outcome - models loss directly
#' @param se_method Method for standard error estimation:
#'   - `"bootstrap"`: Bootstrap standard errors (default)
#'   - `"influence"`: Influence function-based standard errors
#'   - `"none"`: No standard error estimation
#' @param n_boot Number of bootstrap replications (default: 500).
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param cross_fit Logical indicating whether to use cross-fitting for
#'   nuisance model estimation (default: FALSE). Cross-fitting enables
#'   valid inference when using flexible machine learning estimators.
#' @param n_folds Number of folds for cross-fitting (default: 5).
#' @param parallel Logical indicating whether to use parallel processing
#'   for bootstrap (default: FALSE).
#' @param ncores Number of cores for parallel processing (default: NULL,
#'   which uses all available cores minus one).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class `c("cf_mse", "cf_performance")` containing:
#'   \item{estimate}{Point estimate of counterfactual MSE}
#'   \item{se}{Standard error (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound}
#'   \item{ci_upper}{Upper confidence interval bound}
#'   \item{estimator}{Estimator used}
#'   \item{naive_estimate}{Naive MSE for comparison}
#'   \item{n_obs}{Number of observations}
#'   \item{treatment_level}{Counterfactual treatment level}
#'
#' @details
#' The function implements four estimators for the counterfactual MSE:
#'
#' **Naive Estimator**: Simply computes the empirical MSE using observed
#' outcomes. This is biased for the counterfactual estimand when treatment
#' affects outcomes.
#'
#' **Conditional Loss (CL) Estimator**: Models the expected loss conditional
#' on covariates and treatment, then marginalizes. Requires a correctly

#' specified outcome model.
#'
#' **IPW Estimator**: Weights observations by the inverse probability of
#' receiving the counterfactual treatment. Requires a correctly specified
#' propensity score model.
#'
#' **Doubly Robust (DR) Estimator**: Combines CL and IPW approaches. Consistent
#' if either the propensity or outcome model is correctly specified.
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @seealso [cf_auc()], [cf_calibration()], [fit_nuisance()]
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
#' pred <- plogis(-1 + 0.8 * x)  # Predictions from some model
#'
#' # Estimate counterfactual MSE under no treatment
#' result <- cf_mse(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"  # Skip SE for speed in example
#' )
#' print(result)
cf_mse <- function(predictions,
                   outcomes,
                   treatment,
                   covariates,
                   treatment_level = 0,
                   estimator = c("dr", "cl", "ipw", "naive"),
                   propensity_model = NULL,
                   outcome_model = NULL,
                   outcome_type = c("auto", "binary", "continuous"),
                   se_method = c("bootstrap", "influence", "none"),
                   n_boot = 500,
                   conf_level = 0.95,
                   cross_fit = FALSE,
                   n_folds = 5,
                   parallel = FALSE,
                   ncores = NULL,
                   ...) {


  # Input validation
  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)
  outcome_type <- match.arg(outcome_type)

  # Validate inputs

  .validate_inputs(predictions, outcomes, treatment, covariates)

  n <- length(outcomes)

  # Auto-detect outcome type if not specified
  if (outcome_type == "auto") {
    outcome_type <- if (all(outcomes %in% c(0, 1))) "binary" else "continuous"
  }

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Detect if ml_learners are provided
  use_ml_propensity <- is_ml_learner(propensity_model)
  use_ml_outcome <- is_ml_learner(outcome_model)

  # Fit nuisance models if not provided
  if (estimator != "naive") {
    if (cross_fit && estimator == "dr") {
      # Use cross-fitting for DR estimator
      cf_result <- .compute_mse_crossfit(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        treatment_level = treatment_level,
        K = n_folds,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        outcome_type = outcome_type,
        parallel = parallel,
        ...
      )
      estimate <- cf_result$estimate
      nuisance <- list(propensity = NULL, outcome = NULL,
                       cross_fitted = TRUE,
                       ps = cf_result$ps,
                       h = cf_result$h,
                       folds = cf_result$folds)

      # SE from cross-fitting
      if (se_method == "influence") {
        se <- cf_result$se
        z <- qnorm(1 - (1 - conf_level) / 2)
        ci_lower <- estimate - z * se
        ci_upper <- estimate + z * se
      }
    } else {
      nuisance <- .fit_nuisance_models(
        treatment = treatment,
        outcomes = outcomes,
        covariates = covariates,
        treatment_level = treatment_level,
        propensity_model = propensity_model,
        outcome_model = outcome_model,
        estimator = estimator,
        outcome_type = outcome_type,
        predictions = predictions
      )
      estimate <- NULL
    }
  } else {
    nuisance <- list(propensity = NULL, outcome = NULL)
    estimate <- NULL
  }

  # Compute point estimate (if not already computed via cross-fitting)
  if (is.null(estimate)) {
    estimate <- .compute_mse(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = nuisance$propensity,
      outcome_model = nuisance$outcome,
      outcome_type = outcome_type
    )
  }

  # Compute naive estimate for comparison
  naive_estimate <- .compute_mse(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    covariates = covariates,
    treatment_level = treatment_level,
    estimator = "naive",
    propensity_model = NULL,
    outcome_model = NULL
  )

  # Compute standard errors (if not already computed via cross-fitting)
  if (se_method == "bootstrap") {
    boot_result <- .bootstrap_mse(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      n_boot = n_boot,
      conf_level = conf_level,
      parallel = parallel,
      ncores = ncores,
      ...
    )
    se <- boot_result$se
    ci_lower <- boot_result$ci_lower
    ci_upper <- boot_result$ci_upper
  } else if (se_method == "influence") {
    # Check if SE was already computed via cross-fitting
    if (!(cross_fit && estimator == "dr")) {
      se <- .influence_se_mse(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        treatment_level = treatment_level,
        estimator = estimator,
        propensity_model = nuisance$propensity,
        outcome_model = nuisance$outcome
      )
      z <- qnorm(1 - (1 - conf_level) / 2)
      ci_lower <- estimate - z * se
      ci_upper <- estimate + z * se
    }
  }

  # Construct result object
  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    estimator = estimator,
    metric = "mse",
    treatment_level = treatment_level,
    n_obs = n,
    naive_estimate = naive_estimate,
    propensity_model = nuisance$propensity,
    outcome_model = nuisance$outcome,
    call = match.call()
  )

  class(result) <- c("cf_mse", "cf_performance")
  return(result)
}


# Internal function to compute MSE
.compute_mse <- function(predictions, outcomes, treatment, covariates,
                         treatment_level, estimator, propensity_model,
                         outcome_model, outcome_type = "binary") {

  n <- length(outcomes)
  loss <- (outcomes - predictions)^2

  if (estimator == "naive") {
    return(mean(loss))
  }

  # Convert covariates to data frame if needed for prediction
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Get propensity scores for ALL observations
  if (!is.null(propensity_model)) {
    ps <- .predict_nuisance(propensity_model, covariates)
    if (treatment_level == 0) {
      ps <- 1 - ps  # P(A = 0 | X)
    }
    # Truncate extreme propensities for stability
    ps <- pmax(pmin(ps, 0.99), 0.01)
  }

  # Get outcome predictions (conditional loss) for ALL observations
  if (!is.null(outcome_model)) {
    # For binary outcomes, the outcome model predicts E[Y|X,A] = pY
    # and we transform to E[L|X,A] = pY - 2*pred*pY + pred^2
    # For continuous outcomes, the model directly predicts E[L|X,A]
    pY <- .predict_nuisance(outcome_model, covariates)
    if (outcome_type == "binary") {
      # E[(Y - pred)^2 | X] = E[Y | X] - 2*pred*E[Y | X] + pred^2
      # since Y^2 = Y for binary Y
      h <- pY - 2 * predictions * pY + predictions^2
    } else {
      h <- pY
    }
  }

  # Indicator for counterfactual treatment
  I_a <- as.numeric(treatment == treatment_level)

  if (estimator == "cl") {
    # Conditional loss estimator
    return(mean(h))

  } else if (estimator == "ipw") {
    # IPW estimator (Horvitz-Thompson style)
    weights <- I_a / ps
    return(mean(weights * loss))

  } else if (estimator == "dr") {
    # Doubly robust estimator
    augmentation <- I_a / ps * (loss - h)
    return(mean(h + augmentation))
  }
}


# Internal function to fit nuisance models
.fit_nuisance_models <- function(treatment, outcomes, covariates,
                                  treatment_level, propensity_model,
                                  outcome_model, estimator = "dr",
                                  outcome_type = "binary",
                                  predictions = NULL) {

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Fit propensity model if not provided (needed for ipw and dr)
  if (estimator %in% c("ipw", "dr") && is.null(propensity_model)) {
    ps_data <- cbind(A = treatment, covariates)
    propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
  } else if (is_ml_learner(propensity_model)) {
    # ml_learner spec should use cross-fitting - warn user
    warning("ml_learner passed to propensity_model without cross_fit=TRUE. ",
            "Using cross_fit=TRUE is recommended for ML learners.",
            call. = FALSE)
    ps_data <- cbind(A = treatment, covariates)
    propensity_model <- .fit_ml_learner(propensity_model, A ~ .,
                                         data = ps_data, family = "binomial")
  }

  # Fit outcome model if not provided (needed for cl and dr only)
  # For binary outcomes: model E[Y | X, A=a] and transform to loss later
  # For continuous outcomes: model E[L | X, A=a] directly
  if (estimator %in% c("cl", "dr") && is.null(outcome_model)) {
    subset_idx <- treatment == treatment_level

    if (outcome_type == "binary") {
      # Model E[Y | X, A=a] - the transformation to loss happens in .compute_mse
      outcome_data <- cbind(Y = outcomes, covariates)[subset_idx, ]
      outcome_model <- glm(Y ~ ., data = outcome_data, family = binomial())
    } else {
      # Model E[L | X, A=a] directly for continuous outcomes
      loss <- (outcomes - predictions)^2
      outcome_data <- cbind(L = loss, covariates)[subset_idx, ]
      outcome_model <- glm(L ~ ., data = outcome_data, family = gaussian())
    }
    # Store the full data for prediction
    attr(outcome_model, "full_data") <- cbind(Y = outcomes, covariates)
  } else if (is_ml_learner(outcome_model)) {
    # ml_learner spec should use cross-fitting - warn user
    warning("ml_learner passed to outcome_model without cross_fit=TRUE. ",
            "Using cross_fit=TRUE is recommended for ML learners.",
            call. = FALSE)
    subset_idx <- treatment == treatment_level

    if (outcome_type == "binary") {
      outcome_data <- cbind(Y = outcomes, covariates)[subset_idx, ]
      outcome_model <- .fit_ml_learner(outcome_model, Y ~ .,
                                        data = outcome_data, family = "binomial")
    } else {
      loss <- (outcomes - predictions)^2
      outcome_data <- cbind(L = loss, covariates)[subset_idx, ]
      outcome_model <- .fit_ml_learner(outcome_model, L ~ .,
                                        data = outcome_data, family = "gaussian")
    }
  }

  list(propensity = propensity_model, outcome = outcome_model)
}


# Input validation
.validate_inputs <- function(predictions, outcomes, treatment, covariates) {
  n <- length(outcomes)

  if (length(predictions) != n) {
    stop("predictions and outcomes must have the same length")
  }
  if (length(treatment) != n) {
    stop("treatment and outcomes must have the same length")
  }
  if (nrow(covariates) != n) {
    stop("covariates must have the same number of rows as outcomes")
  }
  if (!all(treatment %in% c(0, 1))) {
    stop("treatment must be binary (0/1)")
  }
  if (any(is.na(predictions)) || any(is.na(outcomes)) || any(is.na(treatment))) {
    stop("Missing values not allowed in predictions, outcomes, or treatment")
  }
}
