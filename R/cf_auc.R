#' Estimate Counterfactual Area Under the ROC Curve
#'
#' Estimates the area under the receiver operating characteristic curve (AUC)
#' of a prediction model under a hypothetical intervention where treatment is
#' set to a specific level.
#'
#' @inheritParams cf_mse
#'
#' @return An object of class `c("cf_auc", "cf_performance")` containing:
#'   \item{estimate}{Point estimate of counterfactual AUC}
#'   \item{se}{Standard error (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound}
#'   \item{ci_upper}{Upper confidence interval bound}
#'   \item{estimator}{Estimator used}
#'   \item{naive_estimate}{Naive AUC for comparison}
#'   \item{n_obs}{Number of observations}
#'   \item{treatment_level}{Counterfactual treatment level}
#'
#' @details
#' The counterfactual AUC is defined as the probability that a randomly
#' selected individual with the outcome under the counterfactual intervention
#' has a higher predicted risk than a randomly selected individual without
#' the outcome.
#'
#' The function implements three estimators:
#'
#' **Outcome Model (OM/CL) Estimator**: Weights concordant pairs by the
#' predicted probability of case/non-case status under the counterfactual.
#'
#' **IPW Estimator**: Weights concordant pairs by the inverse probability
#' of treatment.
#'
#' **Doubly Robust (DR) Estimator**: Combines OM and IPW for double robustness.
#' When `cross_fit = TRUE`, uses cross-fitting for valid inference with flexible
#' ML methods (see [ml_learner()]).
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
#' "Estimating the area under the ROC curve when transporting a prediction
#' model to a target population." *Biometrics*, 79(3), 2343-2356.
#' \doi{10.1111/biom.13796}
#'
#' @seealso [cf_mse()], [cf_calibration()], [ml_learner()]
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
#' # Estimate counterfactual AUC under no treatment
#' result <- cf_auc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"
#' )
#' print(result)
cf_auc <- function(predictions,
                   outcomes,
                   treatment,
                   covariates,
                   treatment_level = 0,
                   estimator = c("dr", "cl", "ipw", "naive"),
                   propensity_model = NULL,
                   outcome_model = NULL,
                   se_method = c("bootstrap", "influence", "none"),
                   n_boot = 500,
                   conf_level = 0.95,
                   cross_fit = FALSE,
                   n_folds = 5,
                   parallel = FALSE,
                   ncores = NULL,
                   ps_trim = NULL,
                   ...) {

  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  # Validate inputs
  .validate_inputs(predictions, outcomes, treatment, covariates)

  if (!all(outcomes %in% c(0, 1))) {
    stop("AUC requires binary outcomes (0/1)")
  }

  n <- length(outcomes)

  # Detect if ml_learners are provided
  use_ml_propensity <- is_ml_learner(propensity_model)
  use_ml_outcome <- is_ml_learner(outcome_model)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Fit nuisance models if needed
  if (estimator != "naive") {
    if (cross_fit && estimator == "dr") {
      # Use cross-fitting for DR estimator
      cf_result <- .compute_auc_crossfit(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        treatment_level = treatment_level,
        K = n_folds,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        parallel = parallel,
        ps_trim_spec = ps_trim_spec,
        ...
      )
      estimate <- cf_result$estimate
      nuisance <- list(propensity = NULL, outcome = NULL,
                       cross_fitted = TRUE,
                       ps = cf_result$ps,
                       q = cf_result$q,
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
        outcome_model = outcome_model
      )
      estimate <- NULL
    }
  } else {
    nuisance <- list(propensity = NULL, outcome = NULL)
    estimate <- NULL
  }

  # Compute point estimate (if not already computed via cross-fitting)
  if (is.null(estimate)) {
    estimate <- .compute_auc(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = nuisance$propensity,
      outcome_model = nuisance$outcome,
      ps_trim_spec = ps_trim_spec
    )
  }

  # Naive estimate
  naive_estimate <- .compute_auc_naive(predictions, outcomes)

  # Standard errors (if not already computed via cross-fitting)
  if (se_method == "bootstrap") {
    boot_result <- .bootstrap_auc(
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
      ps_trim = ps_trim
    )
    se <- boot_result$se
    ci_lower <- boot_result$ci_lower
    ci_upper <- boot_result$ci_upper
  } else if (se_method == "influence" && !(cross_fit && estimator == "dr")) {
    se <- .influence_se_auc(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = nuisance$propensity,
      outcome_model = nuisance$outcome,
      ps_trim_spec = ps_trim_spec
    )
    z <- qnorm(1 - (1 - conf_level) / 2)
    ci_lower <- estimate - z * se
    ci_upper <- estimate + z * se
  }

  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    estimator = estimator,
    metric = "auc",
    treatment_level = treatment_level,
    n_obs = n,
    naive_estimate = naive_estimate,
    propensity_model = nuisance$propensity,
    outcome_model = nuisance$outcome,
    call = match.call()
  )

  class(result) <- c("cf_auc", "cf_performance")
  return(result)
}


# Compute naive AUC using rank method
.compute_auc_naive <- function(predictions, outcomes) {
  n1 <- sum(outcomes == 1)
  n0 <- sum(outcomes == 0)

  if (n1 == 0 || n0 == 0) {
    warning("AUC undefined: no cases or no non-cases")
    return(NA_real_)
  }

  r <- rank(c(predictions[outcomes == 1], predictions[outcomes == 0]))
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


# Compute counterfactual AUC
.compute_auc <- function(predictions, outcomes, treatment, covariates,
                         treatment_level, estimator, propensity_model,
                         outcome_model, ps_trim_spec = NULL) {

  if (estimator == "naive") {
    return(.compute_auc_naive(predictions, outcomes))
  }

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)

  # Propensity scores
  if (!is.null(propensity_model)) {
    ps <- .predict_nuisance(propensity_model, as.data.frame(covariates), type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
    # Trim propensity scores for stability
    ps <- .trim_propensity(ps, ps_trim_spec$method, ps_trim_spec$bounds)
  }

  # Outcome probabilities
  if (!is.null(outcome_model)) {
    full_data <- attr(outcome_model, "full_data")
    if (is.null(full_data)) {
      full_data <- cbind(Y = outcomes, as.data.frame(covariates))
    }
    q_hat <- .predict_nuisance(outcome_model, full_data, type = "response")
  }

  # Indicator for treatment level
  I_a <- as.numeric(treatment == treatment_level)

  # Concordance indicator matrix
  ind_f <- outer(predictions, predictions, ">")

  if (estimator == "cl") {
    # Outcome model estimator
    mat_om0 <- outer(q_hat, 1 - q_hat, "*")
    mat_om1 <- mat_om0 * ind_f
    return(sum(mat_om1) / sum(mat_om0))

  } else if (estimator == "ipw") {
    # IPW estimator
    pi_ratio <- I_a / ps
    mat_ipw0 <- outer(I_a * (outcomes == 1), I_a * (outcomes == 0), "*") *
                outer(pi_ratio, pi_ratio, "*")
    mat_ipw1 <- mat_ipw0 * ind_f
    return(sum(mat_ipw1) / sum(mat_ipw0))

  } else if (estimator == "dr") {
    # Doubly robust estimator
    pi_ratio <- I_a / ps

    # IPW matrices
    mat_ipw0 <- outer(I_a * (outcomes == 1), I_a * (outcomes == 0), "*") *
                outer(pi_ratio, pi_ratio, "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    # OM matrices
    mat_om0 <- outer(q_hat, 1 - q_hat, "*")
    mat_om1 <- mat_om0 * ind_f

    # DR matrices
    mat_dr0 <- outer(I_a * pi_ratio * q_hat, I_a * pi_ratio * (1 - q_hat), "*")
    diag(mat_dr0) <- 0
    mat_dr1 <- mat_dr0 * ind_f

    return((sum(mat_ipw1) + sum(mat_om1) - sum(mat_dr1)) /
           (sum(mat_ipw0) + sum(mat_om0) - sum(mat_dr0)))
  }
}
