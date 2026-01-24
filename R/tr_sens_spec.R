#' Estimate (Counterfactual) Sensitivity in the Target Population
#'
#' Estimates the sensitivity (true positive rate) of a binary classifier at
#' one or more thresholds in a target population using data transported from
#' a source population (typically an RCT).
#'
#' @inheritParams tr_mse
#' @param threshold Numeric vector of classification thresholds. Predictions
#'   above this value are classified as positive. Can be a single value or
#'   a vector for computing sensitivity at multiple thresholds simultaneously.
#'   Default is 0.5.
#'
#' @return An object of class `c("tr_sensitivity", "tr_performance")` containing:
#'   \item{estimate}{Point estimate(s) of transportable sensitivity}
#'   \item{se}{Standard error(s) (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound(s)}
#'   \item{ci_upper}{Upper confidence interval bound(s)}
#'   \item{threshold}{Threshold value(s) used}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{naive_estimate}{Naive sensitivity for comparison}
#'   \item{n_target}{Number of target observations}
#'   \item{n_source}{Number of source observations}
#'   \item{treatment_level}{Treatment level}
#'
#' @details
#' Sensitivity (also known as true positive rate or recall) is defined as:
#' \deqn{Sensitivity(c) = P(\hat{Y} > c | Y = 1)}
#'
#' In the transportability setting, we estimate sensitivity in the target
#' population using outcome data from the source population. The function
#' implements three estimators following Steingrimsson et al. (2024):
#'
#' **Outcome Model (OM) Estimator**:
#' \deqn{\hat{\psi}_{sens,om} = \frac{\sum_i I(S_i=0) I(\hat{h}(X_i) > c) \hat{m}(X_i)}{\sum_i I(S_i=0) \hat{m}(X_i)}}
#' where \eqn{\hat{m}(X) \approx P(Y=1|X, R=1)}.
#'
#' **IPW Estimator**:
#' \deqn{\hat{\psi}_{sens,ipw} = \frac{\sum_i I(\hat{h}(X_i) > c, Y_i=1, R_i=1) \hat{w}(X_i)}{\sum_i I(Y_i=1, R_i=1) \hat{w}(X_i)}}
#' where \eqn{\hat{w}(X) \approx P(R=0|X)/P(R=1|X)}.
#'
#' **Doubly Robust (DR) Estimator**: Combines OM and IPW for protection against
#' misspecification of either model.
#'
#' @references
#' Steingrimsson, J. A., Wen, L., Voter, S., & Dahabreh, I. J. (2024).
#' "Interpretable meta-analysis of model or marker performance."
#' *arXiv preprint arXiv:2409.13458*.
#'
#' Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
#' "Counterfactual risk assessments, evaluation, and fairness."
#' *Proceedings of the 2020 Conference on Fairness, Accountability, and
#' Transparency*, 582-593.
#'
#' @seealso [tr_specificity()], [tr_fpr()], [tr_auc()]
#'
#' @export
#'
#' @examples
#' # Generate example data with source (RCT) and target populations
#' set.seed(123)
#' n <- 1000
#' x <- rnorm(n)
#' s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
#' a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate transportable sensitivity at default threshold (0.5)
#' result <- tr_sensitivity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"
#' )
#' print(result)
#'
#' # Estimate at multiple thresholds
#' result_multi <- tr_sensitivity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   threshold = c(0.3, 0.5, 0.7),
#'   se_method = "none"
#' )
#' print(result_multi)
tr_sensitivity <- function(predictions,
                           outcomes,
                           treatment,
                           source,
                           covariates,
                           threshold = 0.5,
                           treatment_level = 0,
                           analysis = c("transport", "joint"),
                           estimator = c("dr", "om", "ipw", "naive"),
                           selection_model = NULL,
                           propensity_model = NULL,
                           outcome_model = NULL,
                           se_method = c("none", "bootstrap", "influence"),
                           n_boot = 200,
                           conf_level = 0.95,
                           stratified_boot = TRUE,
                           cross_fit = FALSE,
                           n_folds = 5,
                           ps_trim = NULL,
                           parallel = FALSE,
                           ncores = NULL,
                           ...) {

  # Input validation
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  # Influence function SE requires cross-fitting
  if (se_method == "influence" && !cross_fit) {
    stop("Influence function standard errors require cross_fit = TRUE")
  }

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  .validate_transport_inputs(predictions, outcomes, treatment, source, covariates)

  # Check binary outcomes
  if (!all(outcomes %in% c(0, 1))) {
    stop("Sensitivity requires binary outcomes (0/1)")
  }

  # Validate threshold
  if (!is.numeric(threshold) || any(threshold < 0) || any(threshold > 1)) {
    stop("threshold must be numeric value(s) between 0 and 1")
  }

  n <- length(outcomes)
  n_target <- sum(source == 0)
  n_source <- sum(source == 1)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Detect if ml_learners are provided
  use_ml_selection <- is_ml_learner(selection_model)
  use_ml_propensity <- is_ml_learner(propensity_model)
  use_ml_outcome <- is_ml_learner(outcome_model)

  # Fit nuisance models or use cross-fitting
  if (estimator != "naive") {
    if (cross_fit && estimator == "dr") {
      # Use cross-fitting for DR estimator
      cf_nuisance <- .cross_fit_transport_nuisance_sens_spec(
        treatment = treatment,
        outcomes = outcomes,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        K = n_folds,
        selection_learner = if (use_ml_selection) selection_model else NULL,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        ps_trim_spec = ps_trim_spec,
        parallel = parallel,
        ...
      )
      nuisance <- list(selection = NULL, propensity = NULL, outcome = NULL,
                       cross_fitted = TRUE,
                       pi_s0 = cf_nuisance$pi_s0,
                       ps = cf_nuisance$ps,
                       m_hat = cf_nuisance$m_hat,
                       folds = cf_nuisance$folds)

      # Compute point estimates and SE using cross-fitted nuisance
      if (se_method == "influence") {
        # Compute estimate and SE together via influence function
        results <- lapply(threshold, function(c) {
          .compute_tr_sensitivity_crossfit(
            predictions = predictions,
            outcomes = outcomes,
            treatment = treatment,
            source = source,
            threshold = c,
            treatment_level = treatment_level,
            analysis = analysis,
            pi_s0 = nuisance$pi_s0,
            ps = nuisance$ps,
            m_hat = nuisance$m_hat,
            return_se = TRUE
          )
        })
        estimate <- sapply(results, function(x) x$estimate)
        se <- sapply(results, function(x) x$se)
        z <- qnorm(1 - (1 - conf_level) / 2)
        ci_lower <- estimate - z * se
        ci_upper <- estimate + z * se
      } else {
        # Just compute point estimates
        estimate <- sapply(threshold, function(c) {
          .compute_tr_sensitivity_crossfit(
            predictions = predictions,
            outcomes = outcomes,
            treatment = treatment,
            source = source,
            threshold = c,
            treatment_level = treatment_level,
            analysis = analysis,
            pi_s0 = nuisance$pi_s0,
            ps = nuisance$ps,
            m_hat = nuisance$m_hat,
            return_se = FALSE
          )
        })
      }
    } else {
      nuisance <- .fit_transport_nuisance_sens_spec(
        treatment = treatment,
        outcomes = outcomes,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        selection_model = selection_model,
        propensity_model = propensity_model,
        outcome_model = outcome_model
      )
      # Compute point estimates for non-crossfit case
      estimate <- sapply(threshold, function(c) {
        .compute_tr_sensitivity(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          source = source,
          covariates = covariates,
          threshold = c,
          treatment_level = treatment_level,
          analysis = analysis,
          estimator = estimator,
          selection_model = nuisance$selection,
          propensity_model = nuisance$propensity,
          outcome_model = nuisance$outcome
        )
      })
    }
  } else {
    nuisance <- list(selection = NULL, propensity = NULL, outcome = NULL)
    # Compute naive estimate
    estimate <- sapply(threshold, function(c) {
      .compute_tr_sensitivity(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        source = source,
        covariates = covariates,
        threshold = c,
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = "naive",
        selection_model = NULL,
        propensity_model = NULL,
        outcome_model = NULL
      )
    })
  }

  # Compute naive estimate for comparison
  naive_estimate <- sapply(threshold, function(c) {
    .compute_tr_sensitivity(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = c,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = "naive",
      selection_model = NULL,
      propensity_model = NULL,
      outcome_model = NULL
    )
  })

  # Compute standard errors via bootstrap
  if (se_method == "bootstrap") {
    boot_result <- .bootstrap_tr_sens_spec(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = estimator,
      metric = "sensitivity",
      n_boot = n_boot,
      conf_level = conf_level,
      stratified = stratified_boot,
      ps_trim_spec = ps_trim_spec,
      parallel = parallel,
      ncores = ncores,
      ...
    )
    se <- boot_result$se
    ci_lower <- boot_result$ci_lower
    ci_upper <- boot_result$ci_upper
  }

  # Construct result object
  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    threshold = threshold,
    estimator = estimator,
    analysis = analysis,
    metric = "sensitivity",
    treatment_level = treatment_level,
    n_target = n_target,
    n_source = n_source,
    n_obs = n,
    naive_estimate = naive_estimate,
    selection_model = nuisance$selection,
    propensity_model = nuisance$propensity,
    outcome_model = nuisance$outcome,
    call = match.call()
  )

  class(result) <- c("tr_sensitivity", "tr_performance")
  return(result)
}


#' Estimate (Counterfactual) Specificity in the Target Population
#'
#' Estimates the specificity (true negative rate) of a binary classifier at
#' one or more thresholds in a target population using data transported from
#' a source population (typically an RCT).
#'
#' @inheritParams tr_sensitivity
#'
#' @return An object of class `c("tr_specificity", "tr_performance")` containing:
#'   \item{estimate}{Point estimate(s) of transportable specificity}
#'   \item{se}{Standard error(s) (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound(s)}
#'   \item{ci_upper}{Upper confidence interval bound(s)}
#'   \item{threshold}{Threshold value(s) used}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{naive_estimate}{Naive specificity for comparison}
#'   \item{n_target}{Number of target observations}
#'   \item{n_source}{Number of source observations}
#'   \item{treatment_level}{Treatment level}
#'
#' @details
#' Specificity (also known as true negative rate) is defined as:
#' \deqn{Specificity(c) = P(\hat{Y} \le c | Y = 0)}
#'
#' In the transportability setting, we estimate specificity in the target
#' population using outcome data from the source population. The estimators
#' mirror those for sensitivity (see [tr_sensitivity()]).
#'
#' @references
#' Steingrimsson, J. A., Wen, L., Voter, S., & Dahabreh, I. J. (2024).
#' "Interpretable meta-analysis of model or marker performance."
#' *arXiv preprint arXiv:2409.13458*.
#'
#' @seealso [tr_sensitivity()], [tr_fpr()], [tr_auc()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 1000
#' x <- rnorm(n)
#' s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
#' a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate transportable specificity
#' result <- tr_specificity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"
#' )
#' print(result)
tr_specificity <- function(predictions,
                           outcomes,
                           treatment,
                           source,
                           covariates,
                           threshold = 0.5,
                           treatment_level = 0,
                           analysis = c("transport", "joint"),
                           estimator = c("dr", "om", "ipw", "naive"),
                           selection_model = NULL,
                           propensity_model = NULL,
                           outcome_model = NULL,
                           se_method = c("none", "bootstrap", "influence"),
                           n_boot = 200,
                           conf_level = 0.95,
                           stratified_boot = TRUE,
                           cross_fit = FALSE,
                           n_folds = 5,
                           ps_trim = NULL,
                           parallel = FALSE,
                           ncores = NULL,
                           ...) {

  # Input validation
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  # Influence function SE requires cross-fitting
  if (se_method == "influence" && !cross_fit) {
    stop("Influence function standard errors require cross_fit = TRUE")
  }

  # Parse ps_trim specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  .validate_transport_inputs(predictions, outcomes, treatment, source, covariates)

  # Check binary outcomes
  if (!all(outcomes %in% c(0, 1))) {
    stop("Specificity requires binary outcomes (0/1)")
  }

  # Validate threshold
  if (!is.numeric(threshold) || any(threshold < 0) || any(threshold > 1)) {
    stop("threshold must be numeric value(s) between 0 and 1")
  }

  n <- length(outcomes)
  n_target <- sum(source == 0)
  n_source <- sum(source == 1)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Detect if ml_learners are provided
  use_ml_selection <- is_ml_learner(selection_model)
  use_ml_propensity <- is_ml_learner(propensity_model)
  use_ml_outcome <- is_ml_learner(outcome_model)

  # Fit nuisance models or use cross-fitting
  if (estimator != "naive") {
    if (cross_fit && estimator == "dr") {
      # Use cross-fitting for DR estimator
      cf_nuisance <- .cross_fit_transport_nuisance_sens_spec(
        treatment = treatment,
        outcomes = outcomes,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        K = n_folds,
        selection_learner = if (use_ml_selection) selection_model else NULL,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        ps_trim_spec = ps_trim_spec,
        parallel = parallel,
        ...
      )
      nuisance <- list(selection = NULL, propensity = NULL, outcome = NULL,
                       cross_fitted = TRUE,
                       pi_s0 = cf_nuisance$pi_s0,
                       ps = cf_nuisance$ps,
                       m_hat = cf_nuisance$m_hat,
                       folds = cf_nuisance$folds)

      # Compute point estimates and SE using cross-fitted nuisance
      if (se_method == "influence") {
        # Compute estimate and SE together via influence function
        results <- lapply(threshold, function(c) {
          .compute_tr_specificity_crossfit(
            predictions = predictions,
            outcomes = outcomes,
            treatment = treatment,
            source = source,
            threshold = c,
            treatment_level = treatment_level,
            analysis = analysis,
            pi_s0 = nuisance$pi_s0,
            ps = nuisance$ps,
            m_hat = nuisance$m_hat,
            return_se = TRUE
          )
        })
        estimate <- sapply(results, function(x) x$estimate)
        se <- sapply(results, function(x) x$se)
        z <- qnorm(1 - (1 - conf_level) / 2)
        ci_lower <- estimate - z * se
        ci_upper <- estimate + z * se
      } else {
        # Just compute point estimates
        estimate <- sapply(threshold, function(c) {
          .compute_tr_specificity_crossfit(
            predictions = predictions,
            outcomes = outcomes,
            treatment = treatment,
            source = source,
            threshold = c,
            treatment_level = treatment_level,
            analysis = analysis,
            pi_s0 = nuisance$pi_s0,
            ps = nuisance$ps,
            m_hat = nuisance$m_hat,
            return_se = FALSE
          )
        })
      }
    } else {
      nuisance <- .fit_transport_nuisance_sens_spec(
        treatment = treatment,
        outcomes = outcomes,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        selection_model = selection_model,
        propensity_model = propensity_model,
        outcome_model = outcome_model
      )
      # Compute point estimates for non-crossfit case
      estimate <- sapply(threshold, function(c) {
        .compute_tr_specificity(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          source = source,
          covariates = covariates,
          threshold = c,
          treatment_level = treatment_level,
          analysis = analysis,
          estimator = estimator,
          selection_model = nuisance$selection,
          propensity_model = nuisance$propensity,
          outcome_model = nuisance$outcome
        )
      })
    }
  } else {
    nuisance <- list(selection = NULL, propensity = NULL, outcome = NULL)
    # Compute naive estimate
    estimate <- sapply(threshold, function(c) {
      .compute_tr_specificity(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        source = source,
        covariates = covariates,
        threshold = c,
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = "naive",
        selection_model = NULL,
        propensity_model = NULL,
        outcome_model = NULL
      )
    })
  }

  # Compute naive estimate for comparison
  naive_estimate <- sapply(threshold, function(c) {
    .compute_tr_specificity(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = c,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = "naive",
      selection_model = NULL,
      propensity_model = NULL,
      outcome_model = NULL
    )
  })

  # Compute standard errors via bootstrap
  if (se_method == "bootstrap") {
    boot_result <- .bootstrap_tr_sens_spec(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = estimator,
      metric = "specificity",
      n_boot = n_boot,
      conf_level = conf_level,
      stratified = stratified_boot,
      ps_trim_spec = ps_trim_spec,
      parallel = parallel,
      ncores = ncores,
      ...
    )
    se <- boot_result$se
    ci_lower <- boot_result$ci_lower
    ci_upper <- boot_result$ci_upper
  }

  # Construct result object
  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    threshold = threshold,
    estimator = estimator,
    analysis = analysis,
    metric = "specificity",
    treatment_level = treatment_level,
    n_target = n_target,
    n_source = n_source,
    n_obs = n,
    naive_estimate = naive_estimate,
    selection_model = nuisance$selection,
    propensity_model = nuisance$propensity,
    outcome_model = nuisance$outcome,
    call = match.call()
  )

  class(result) <- c("tr_specificity", "tr_performance")
  return(result)
}


# ==============================================================================
# Aliases for ML terminology
# ==============================================================================

#' @rdname tr_sensitivity
#' @export
tr_tpr <- tr_sensitivity

#' @rdname tr_specificity
#' @export
tr_tnr <- tr_specificity

#' Estimate (Counterfactual) False Positive Rate in the Target Population
#'
#' Estimates the false positive rate (1 - specificity) of a binary classifier.
#' This is a convenience wrapper around [tr_specificity()].
#'
#' @inheritParams tr_specificity
#'
#' @return An object of class `c("tr_fpr", "tr_performance")` with the same
#'   structure as [tr_specificity()], but with `estimate` containing 1 - specificity.
#'
#' @seealso [tr_specificity()], [tr_sensitivity()]
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' s <- rbinom(n, 1, 0.6)
#' a <- rbinom(n, 1, 0.5)
#' y <- rbinom(n, 1, plogis(-1 + x))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' result <- tr_fpr(
#'   predictions = pred, outcomes = y, treatment = a,
#'   source = s, covariates = data.frame(x = x),
#'   se_method = "none"
#' )
#' print(result)
tr_fpr <- function(predictions,
                   outcomes,
                   treatment,
                   source,
                   covariates,
                   threshold = 0.5,
                   treatment_level = 0,
                   analysis = c("transport", "joint"),
                   estimator = c("dr", "om", "ipw", "naive"),
                   selection_model = NULL,
                   propensity_model = NULL,
                   outcome_model = NULL,
                   se_method = c("none", "bootstrap", "influence"),
                   n_boot = 200,
                   conf_level = 0.95,
                   stratified_boot = TRUE,
                   cross_fit = FALSE,
                   n_folds = 5,
                   ps_trim = NULL,
                   parallel = FALSE,
                   ncores = NULL,
                   ...) {

  # Call specificity and transform

  spec_result <- tr_specificity(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    source = source,
    covariates = covariates,
    threshold = threshold,
    treatment_level = treatment_level,
    analysis = analysis,
    estimator = estimator,
    selection_model = selection_model,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    se_method = se_method,
    n_boot = n_boot,
    conf_level = conf_level,
    stratified_boot = stratified_boot,
    ps_trim = ps_trim,
    parallel = parallel,
    ncores = ncores,
    ...
  )

  # Transform to FPR = 1 - specificity
  spec_result$estimate <- 1 - spec_result$estimate
  spec_result$naive_estimate <- 1 - spec_result$naive_estimate

  # Swap CI bounds since we're subtracting from 1
  if (!is.null(spec_result$ci_lower)) {
    old_lower <- spec_result$ci_lower
    spec_result$ci_lower <- 1 - spec_result$ci_upper
    spec_result$ci_upper <- 1 - old_lower
  }

  spec_result$metric <- "fpr"
  class(spec_result) <- c("tr_fpr", "tr_performance")

  return(spec_result)
}


# ==============================================================================
# Internal Functions: Cross-Fitting for Nuisance Models
# ==============================================================================

#' Cross-fit nuisance models for transportability sensitivity/specificity
#'
#' Implements K-fold cross-fitting for nuisance model estimation in
#' transportability settings for sensitivity and specificity estimation.
#'
#' @param treatment Numeric vector of treatment indicators.
#' @param outcomes Numeric vector of observed outcomes.
#' @param source Numeric vector indicating source (1) or target (0) population.
#' @param covariates Matrix or data frame of covariates.
#' @param treatment_level Counterfactual treatment level.
#' @param analysis Either "transport" or "joint".
#' @param K Number of folds for cross-fitting (default: 5).
#' @param selection_learner Optional ml_learner for selection model.
#' @param propensity_learner Optional ml_learner for propensity model.
#' @param outcome_learner Optional ml_learner for outcome model.
#' @param parallel Logical for parallel processing.
#' @param ps_trim_spec Parsed propensity score trimming specification.
#' @param ... Additional arguments.
#'
#' @return List containing cross-fitted nuisance function predictions.
#'
#' @keywords internal
.cross_fit_transport_nuisance_sens_spec <- function(treatment, outcomes, source,
                                                     covariates, treatment_level,
                                                     analysis, K = 5,
                                                     selection_learner = NULL,
                                                     propensity_learner = NULL,
                                                     outcome_learner = NULL,
                                                     ps_trim_spec = NULL,
                                                     parallel = FALSE,
                                                     ...) {

  n <- length(outcomes)

  # Parse ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Create fold assignments (stratified by source)
  folds <- integer(n)
  idx_s0 <- which(source == 0)
  idx_s1 <- which(source == 1)
  folds[idx_s0] <- sample(rep(1:K, length.out = length(idx_s0)))
  folds[idx_s1] <- sample(rep(1:K, length.out = length(idx_s1)))

  # Initialize output vectors
  pi_s0_cf <- numeric(n)  # Cross-fitted P(S=0|X)
  ps_cf <- numeric(n)     # Cross-fitted propensity scores
  m_hat_cf <- numeric(n)  # Cross-fitted outcome probabilities P(Y=1|X,A=a)

  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)
  I_s1 <- as.numeric(source == 1)

  for (k in 1:K) {
    train_idx <- which(folds != k)
    val_idx <- which(folds == k)

    # --- Selection model: P(S=0|X) ---
    sel_data <- data.frame(S0 = I_s0[train_idx], covariates[train_idx, , drop = FALSE])

    if (!is.null(selection_learner) && is_ml_learner(selection_learner)) {
      sel_model <- .fit_ml_learner(selection_learner, S0 ~ .,
                                    data = sel_data, family = "binomial")
      pi_s0_pred <- .predict_ml_learner(sel_model, covariates[val_idx, , drop = FALSE])
    } else {
      sel_model <- glm(S0 ~ ., data = sel_data, family = binomial())
      pi_s0_pred <- predict(sel_model, newdata = covariates[val_idx, , drop = FALSE],
                            type = "response")
    }
    pi_s0_cf[val_idx] <- .trim_propensity(pi_s0_pred, ps_trim_spec$method, ps_trim_spec$bounds)

    # --- Propensity model ---
    if (analysis == "transport") {
      # P(A=1|X, S=1) - fit on source population only
      train_s1 <- train_idx[source[train_idx] == 1]
      ps_data <- data.frame(A = treatment[train_s1], covariates[train_s1, , drop = FALSE])
    } else {
      # P(A=1|X) - fit on all training data
      ps_data <- data.frame(A = treatment[train_idx], covariates[train_idx, , drop = FALSE])
    }

    if (!is.null(propensity_learner) && is_ml_learner(propensity_learner)) {
      ps_model <- .fit_ml_learner(propensity_learner, A ~ .,
                                   data = ps_data, family = "binomial")
      ps_pred <- .predict_ml_learner(ps_model, covariates[val_idx, , drop = FALSE])
    } else {
      ps_model <- glm(A ~ ., data = ps_data, family = binomial())
      ps_pred <- predict(ps_model, newdata = covariates[val_idx, , drop = FALSE],
                         type = "response")
    }

    if (treatment_level == 0) {
      ps_pred <- 1 - ps_pred
    }
    ps_cf[val_idx] <- .trim_propensity(ps_pred, ps_trim_spec$method, ps_trim_spec$bounds)

    # --- Outcome model: P(Y=1|X, A=a) ---
    if (analysis == "transport") {
      # Model on source data with treatment_level
      train_s1_a <- train_idx[source[train_idx] == 1 & treatment[train_idx] == treatment_level]
      om_data <- data.frame(Y = outcomes[train_s1_a], covariates[train_s1_a, , drop = FALSE])
    } else {
      # Model on all data with treatment_level
      train_a <- train_idx[treatment[train_idx] == treatment_level]
      om_data <- data.frame(Y = outcomes[train_a], covariates[train_a, , drop = FALSE])
    }

    if (!is.null(outcome_learner) && is_ml_learner(outcome_learner)) {
      om_model <- .fit_ml_learner(outcome_learner, Y ~ .,
                                   data = om_data, family = "binomial")
      m_pred <- .predict_ml_learner(om_model, covariates[val_idx, , drop = FALSE])
    } else {
      om_model <- glm(Y ~ ., data = om_data, family = binomial())
      m_pred <- predict(om_model, newdata = covariates[val_idx, , drop = FALSE],
                        type = "response")
    }
    # Clip to [0, 1] bounds (not propensity trimming - this is an outcome model)
    m_hat_cf[val_idx] <- pmin(pmax(m_pred, 0), 1)
  }

  list(
    pi_s0 = pi_s0_cf,
    ps = ps_cf,
    m_hat = m_hat_cf,
    folds = folds
  )
}


#' Compute transportable sensitivity with cross-fitted nuisance functions
#'
#' Computes the DR transportable sensitivity estimator using cross-fitted
#' nuisance functions and returns both point estimate and influence
#' function-based standard error.
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of observed outcomes.
#' @param treatment Numeric vector of treatment indicators.
#' @param source Numeric vector indicating source (1) or target (0) population.
#' @param threshold Classification threshold.
#' @param treatment_level Counterfactual treatment level.
#' @param analysis Either "transport" or "joint".
#' @param pi_s0 Cross-fitted P(S=0|X).
#' @param ps Cross-fitted propensity scores.
#' @param m_hat Cross-fitted outcome probabilities P(Y=1|X,A=a).
#' @param return_se Logical; if TRUE, returns list with estimate and SE.
#'
#' @return If return_se=FALSE, numeric estimate. If TRUE, list with estimate and se.
#' @noRd
.compute_tr_sensitivity_crossfit <- function(predictions, outcomes, treatment, source,
                                              threshold, treatment_level, analysis,
                                              pi_s0, ps, m_hat, return_se = FALSE) {

  n <- length(outcomes)
  n0 <- sum(source == 0)

  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)
  I_s1 <- as.numeric(source == 1)
  I_pos <- as.numeric(predictions > threshold)  # Positive predictions
  I_y1 <- as.numeric(outcomes == 1)             # Cases

  pi_s1 <- 1 - pi_s0

  # DR estimator for sensitivity: P(pred > c | Y=1, target)
  # = P(pred > c, Y=1 | target) / P(Y=1 | target)

  if (analysis == "transport") {
    # Transport DR estimator
    # w(X) = I(S=1) * I(A=a) * P(S=0|X) / (P(S=1|X) * P(A=a|X,S=1))
    w <- I_s1 * I_a * pi_s0 / (pi_s1 * ps)

    # Numerator: E[I(pred>c) * m(X) | S=0] + IPW augmentation
    term1_num <- sum(I_s0 * I_pos * m_hat)
    term2_num <- sum(w * (I_pos * I_y1 - I_pos * m_hat))
    mu_1 <- (term1_num + term2_num) / n0

    # Denominator: E[m(X) | S=0] + IPW augmentation
    term1_den <- sum(I_s0 * m_hat)
    term2_den <- sum(w * (I_y1 - m_hat))
    mu_0 <- (term1_den + term2_den) / n0

  } else {
    # Joint DR estimator
    w <- I_a * pi_s0 / ps

    term1_num <- sum(I_s0 * I_pos * m_hat)
    term2_num <- sum(w * (I_pos * I_y1 - I_pos * m_hat))
    mu_1 <- (term1_num + term2_num) / n0

    term1_den <- sum(I_s0 * m_hat)
    term2_den <- sum(w * (I_y1 - m_hat))
    mu_0 <- (term1_den + term2_den) / n0
  }

  # Handle edge cases
  if (mu_0 <= 0) {
    if (return_se) return(list(estimate = NA_real_, se = NA_real_))
    return(NA_real_)
  }

  estimate <- mu_1 / mu_0

  if (!return_se) return(estimate)

  # Influence function for ratio estimator: phi = (phi_1 - psi * phi_0) / mu_0
  # phi_1: influence function for numerator (transportability DR)
  # phi_0: influence function for denominator (transportability DR)

  # For transportability: phi_i = (1/n0) * [I(S=0)*g(X) + w(X)*(obs - g(X)) - psi*I(S=0)]
  # where g(X) is outcome model, obs is observed quantity

  # phi_1_i: numerator influence function
  phi_1 <- (I_s0 * I_pos * m_hat + w * (I_pos * I_y1 - I_pos * m_hat) - mu_1 * I_s0) / mean(I_s0)

  # phi_0_i: denominator influence function
  phi_0 <- (I_s0 * m_hat + w * (I_y1 - m_hat) - mu_0 * I_s0) / mean(I_s0)

  # Influence function for ratio (delta method)
  phi <- (phi_1 - estimate * phi_0) / mu_0

  se <- sqrt(var(phi) / n)

  list(estimate = estimate, se = se)
}


#' Compute transportable specificity with cross-fitted nuisance functions
#'
#' Computes the DR transportable specificity estimator using cross-fitted
#' nuisance functions and returns both point estimate and influence
#' function-based standard error.
#'
#' @inheritParams .compute_tr_sensitivity_crossfit
#' @param return_se Logical; if TRUE, returns list with estimate and SE.
#'
#' @return If return_se=FALSE, numeric estimate. If TRUE, list with estimate and se.
#' @noRd
.compute_tr_specificity_crossfit <- function(predictions, outcomes, treatment, source,
                                              threshold, treatment_level, analysis,
                                              pi_s0, ps, m_hat, return_se = FALSE) {

  n <- length(outcomes)
  n0 <- sum(source == 0)

  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)
  I_s1 <- as.numeric(source == 1)
  I_neg <- as.numeric(predictions <= threshold)  # Negative predictions
  I_y0 <- as.numeric(outcomes == 0)              # Non-cases

  pi_s1 <- 1 - pi_s0
  one_minus_m <- 1 - m_hat  # P(Y=0|X)

  # DR estimator for specificity: P(pred <= c | Y=0, target)

  if (analysis == "transport") {
    # Transport DR estimator
    w <- I_s1 * I_a * pi_s0 / (pi_s1 * ps)

    # Numerator: E[I(pred <= c) * (1-m(X)) | S=0] + IPW augmentation
    term1_num <- sum(I_s0 * I_neg * one_minus_m)
    term2_num <- sum(w * (I_neg * I_y0 - I_neg * one_minus_m))
    mu_1 <- (term1_num + term2_num) / n0

    # Denominator: E[1-m(X) | S=0] + IPW augmentation
    term1_den <- sum(I_s0 * one_minus_m)
    term2_den <- sum(w * (I_y0 - one_minus_m))
    mu_0 <- (term1_den + term2_den) / n0

  } else {
    # Joint DR estimator
    w <- I_a * pi_s0 / ps

    term1_num <- sum(I_s0 * I_neg * one_minus_m)
    term2_num <- sum(w * (I_neg * I_y0 - I_neg * one_minus_m))
    mu_1 <- (term1_num + term2_num) / n0

    term1_den <- sum(I_s0 * one_minus_m)
    term2_den <- sum(w * (I_y0 - one_minus_m))
    mu_0 <- (term1_den + term2_den) / n0
  }

  # Handle edge cases
  if (mu_0 <= 0) {
    if (return_se) return(list(estimate = NA_real_, se = NA_real_))
    return(NA_real_)
  }

  estimate <- mu_1 / mu_0

  if (!return_se) return(estimate)

  # Influence function for ratio estimator
  # phi_1_i: numerator influence function
  phi_1 <- (I_s0 * I_neg * one_minus_m + w * (I_neg * I_y0 - I_neg * one_minus_m) - mu_1 * I_s0) / mean(I_s0)

  # phi_0_i: denominator influence function
  phi_0 <- (I_s0 * one_minus_m + w * (I_y0 - one_minus_m) - mu_0 * I_s0) / mean(I_s0)

  # Influence function for ratio (delta method)
  phi <- (phi_1 - estimate * phi_0) / mu_0

  se <- sqrt(var(phi) / n)

  list(estimate = estimate, se = se)
}


# ==============================================================================
# Internal Functions: Nuisance Model Fitting
# ==============================================================================

#' Fit nuisance models for sensitivity/specificity estimation
#' @noRd
.fit_transport_nuisance_sens_spec <- function(treatment, outcomes, source,
                                               covariates, treatment_level, analysis,
                                               selection_model, propensity_model,
                                               outcome_model) {

  df <- cbind(
    Y = outcomes,
    A = treatment,
    S = source,
    as.data.frame(covariates)
  )
  covariate_names <- names(as.data.frame(covariates))


  # Selection model: P(S=0|X) - target probability
  if (is.null(selection_model)) {
    selection_formula <- as.formula(
      paste("I(S == 0) ~", paste(covariate_names, collapse = " + "))
    )
    selection_model <- glm(selection_formula, data = df, family = binomial())
  }

  # Propensity model: P(A=a|X, S=1)
  if (is.null(propensity_model)) {
    if (analysis == "transport") {
      propensity_formula <- as.formula(
        paste("A ~", paste(covariate_names, collapse = " + "))
      )
      propensity_model <- glm(propensity_formula,
                              data = df[df$S == 1, ],
                              family = binomial())
    } else {
      propensity_formula <- as.formula(
        paste("A ~", paste(covariate_names, collapse = " + "))
      )
      propensity_model <- glm(propensity_formula, data = df, family = binomial())
    }
  }

  # Outcome model: E[Y|X, A=a, S=1] = P(Y=1|X, A=a, S=1)
  # This is m_hat(X) in the Steingrimsson paper
  if (is.null(outcome_model)) {
    outcome_formula <- as.formula(
      paste("Y ~", paste(covariate_names, collapse = " + "))
    )

    if (analysis == "transport") {
      subset_data <- df[df$S == 1 & df$A == treatment_level, ]
      if (nrow(subset_data) < 5) {
        warning("Very few observations for outcome model fitting")
      }
      outcome_model <- glm(outcome_formula, data = subset_data, family = binomial())
    } else {
      subset_data <- df[df$A == treatment_level, ]
      outcome_model <- glm(outcome_formula, data = subset_data, family = binomial())
    }
  }

  list(
    selection = selection_model,
    propensity = propensity_model,
    outcome = outcome_model
  )
}


# ==============================================================================
# Internal Functions: Sensitivity/Specificity Computation
# ==============================================================================

#' Compute transportable sensitivity
#' @noRd
.compute_tr_sensitivity <- function(predictions, outcomes, treatment, source,
                                     covariates, threshold, treatment_level,
                                     analysis, estimator, selection_model,
                                     propensity_model, outcome_model) {

  if (estimator == "naive") {
    return(.compute_tr_sens_spec_naive(
      predictions, outcomes, treatment, source,
      threshold, treatment_level, analysis, metric = "sensitivity"
    ))
  }

  # Refit nuisance models if they are NULL (e.g., in bootstrap)
  if (is.null(selection_model) || is.null(propensity_model) || is.null(outcome_model)) {
    models <- .fit_transport_nuisance_sens_spec(
      treatment = treatment,
      outcomes = outcomes,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      analysis = analysis,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
    selection_model <- models$selection
    propensity_model <- models$propensity
    outcome_model <- models$outcome
  }

  if (analysis == "transport") {
    return(.tr_sens_transport(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      metric = "sensitivity"
    ))
  } else {
    return(.tr_sens_joint(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      metric = "sensitivity"
    ))
  }
}


#' Compute transportable specificity
#' @noRd
.compute_tr_specificity <- function(predictions, outcomes, treatment, source,
                                     covariates, threshold, treatment_level,
                                     analysis, estimator, selection_model,
                                     propensity_model, outcome_model) {

  if (estimator == "naive") {
    return(.compute_tr_sens_spec_naive(
      predictions, outcomes, treatment, source,
      threshold, treatment_level, analysis, metric = "specificity"
    ))
  }

  # Refit nuisance models if they are NULL (e.g., in bootstrap)
  if (is.null(selection_model) || is.null(propensity_model) || is.null(outcome_model)) {
    models <- .fit_transport_nuisance_sens_spec(
      treatment = treatment,
      outcomes = outcomes,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      analysis = analysis,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
    selection_model <- models$selection
    propensity_model <- models$propensity
    outcome_model <- models$outcome
  }

  if (analysis == "transport") {
    return(.tr_sens_transport(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      metric = "specificity"
    ))
  } else {
    return(.tr_sens_joint(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      threshold = threshold,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      metric = "specificity"
    ))
  }
}


#' Compute naive sensitivity/specificity
#' @noRd
.compute_tr_sens_spec_naive <- function(predictions, outcomes, treatment, source,
                                         threshold, treatment_level, analysis,
                                         metric) {
  # For transport: use source data with treatment level
  # For joint: use all data with treatment level
  if (analysis == "transport") {
    idx <- source == 1 & treatment == treatment_level
  } else {
    idx <- treatment == treatment_level
  }

  preds_subset <- predictions[idx]
  outcomes_subset <- outcomes[idx]

  if (metric == "sensitivity") {
    # P(pred > c | Y = 1)
    cases <- outcomes_subset == 1
    if (sum(cases) == 0) {
      warning("No cases in relevant subset")
      return(NA_real_)
    }
    return(mean(preds_subset[cases] > threshold))
  } else {
    # specificity: P(pred <= c | Y = 0)
    controls <- outcomes_subset == 0
    if (sum(controls) == 0) {
      warning("No controls in relevant subset")
      return(NA_real_)
    }
    return(mean(preds_subset[controls] <= threshold))
  }
}


#' Transport sensitivity/specificity estimators
#' @noRd
.tr_sens_transport <- function(predictions, outcomes, treatment, source,
                                covariates, threshold, treatment_level, estimator,
                                selection_model, propensity_model, outcome_model,
                                metric) {

  n <- length(outcomes)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X)
  pi_s0 <- predict(selection_model, newdata = df, type = "response")
  pi_s1 <- 1 - pi_s0

  # Get propensity scores P(A=a|X, S=1)
  ps_s1 <- predict(propensity_model, newdata = df, type = "response")
  if (treatment_level == 0) {
    ps_s1 <- 1 - ps_s1
  }

  # Get outcome probabilities m_hat = P(Y=1|X, A=a, S=1)
  m_hat <- predict(outcome_model, newdata = df, type = "response")

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_s1 <- as.numeric(source == 1)
  I_s0 <- as.numeric(source == 0)

  # Classifier predictions above threshold
  I_pos <- as.numeric(predictions > threshold)
  I_neg <- as.numeric(predictions <= threshold)

  # IPW weights: w(X) = P(S=0|X) / P(S=1|X)
  w_hat <- pi_s0 / pi_s1

  # For sensitivity: condition on Y=1

  # For specificity: condition on Y=0
  if (metric == "sensitivity") {
    I_outcome <- as.numeric(outcomes == 1)
    m_outcome <- m_hat  # P(Y=1|X)
    I_class <- I_pos    # pred > c
  } else {
    I_outcome <- as.numeric(outcomes == 0)
    m_outcome <- 1 - m_hat  # P(Y=0|X)
    I_class <- I_neg        # pred <= c
  }

  if (estimator == "om") {
    # Outcome model estimator (equation 5 in Steingrimsson et al.)
    # sum_i I(S_i=0) * I(pred_i > c) * m_hat(X_i) / sum_i I(S_i=0) * m_hat(X_i)
    num <- sum(I_s0 * I_class * m_outcome)
    denom <- sum(I_s0 * m_outcome)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "ipw") {
    # IPW estimator (equation 6 in Steingrimsson et al.)
    # Weight source observations by w(X) and condition on outcome
    ipw_weight <- I_s1 * I_a * w_hat / ps_s1
    # Clip extreme weights
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    num <- sum(I_class * I_outcome * ipw_weight)
    denom <- sum(I_outcome * ipw_weight)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "dr") {
    # Doubly robust estimator (from influence function derivation)
    # Combines OM and IPW with augmentation term

    ipw_weight <- I_s1 * I_a * w_hat / ps_s1
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    # Numerator: DR for P(pred > c, Y=1) in target
    # = OM + IPW - augmentation
    num_om <- sum(I_s0 * I_class * m_outcome)
    num_ipw <- sum(ipw_weight * I_class * I_outcome)
    num_aug <- sum(ipw_weight * I_class * m_outcome)

    # Denominator: DR for P(Y=1) in target
    denom_om <- sum(I_s0 * m_outcome)
    denom_ipw <- sum(ipw_weight * I_outcome)
    denom_aug <- sum(ipw_weight * m_outcome)

    num <- num_om + num_ipw - num_aug
    denom <- denom_om + denom_ipw - denom_aug

    if (denom == 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


#' Joint sensitivity/specificity estimators
#' @noRd
.tr_sens_joint <- function(predictions, outcomes, treatment, source,
                            covariates, threshold, treatment_level, estimator,
                            selection_model, propensity_model, outcome_model,
                            metric) {

  n <- length(outcomes)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X)
  pi_s0 <- predict(selection_model, newdata = df, type = "response")

  # Get propensity scores P(A=a|X)
  ps <- predict(propensity_model, newdata = df, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }

  # Get outcome probabilities m_hat = P(Y=1|X, A=a)
  m_hat <- predict(outcome_model, newdata = df, type = "response")

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)

  # Classifier predictions
  I_pos <- as.numeric(predictions > threshold)
  I_neg <- as.numeric(predictions <= threshold)

  # IPW weights using target selection
  w_hat <- pi_s0 / (1 - pi_s0)
  w_hat[source == 0] <- 1  # Target observations get weight 1

  if (metric == "sensitivity") {
    I_outcome <- as.numeric(outcomes == 1)
    m_outcome <- m_hat
    I_class <- I_pos
  } else {
    I_outcome <- as.numeric(outcomes == 0)
    m_outcome <- 1 - m_hat
    I_class <- I_neg
  }

  if (estimator == "om") {
    # Outcome model: use all target observations weighted by outcome probability
    num <- sum(pi_s0 * I_class * m_outcome)
    denom <- sum(pi_s0 * m_outcome)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "ipw") {
    # IPW: weight by selection and propensity
    ipw_weight <- I_a * pi_s0 / ((1 - pi_s0) * ps)
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    num <- sum(I_class * I_outcome * ipw_weight)
    denom <- sum(I_outcome * ipw_weight)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "dr") {
    # Doubly robust
    ipw_weight <- I_a * pi_s0 / ((1 - pi_s0) * ps)
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    num_om <- sum(pi_s0 * I_class * m_outcome)
    num_ipw <- sum(ipw_weight * I_class * I_outcome)
    num_aug <- sum(ipw_weight * I_class * m_outcome)

    denom_om <- sum(pi_s0 * m_outcome)
    denom_ipw <- sum(ipw_weight * I_outcome)
    denom_aug <- sum(ipw_weight * m_outcome)

    num <- num_om + num_ipw - num_aug
    denom <- denom_om + denom_ipw - denom_aug

    if (denom == 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


# ==============================================================================
# Internal Functions: Bootstrap
# ==============================================================================

#' Bootstrap for sensitivity/specificity
#' @noRd
.bootstrap_tr_sens_spec <- function(predictions, outcomes, treatment, source,
                                     covariates, threshold, treatment_level,
                                     analysis, estimator, metric, n_boot,
                                     conf_level, stratified, parallel, ncores, ...) {

  n <- length(outcomes)
  n_thresholds <- length(threshold)
  alpha <- 1 - conf_level

  # Bootstrap function for single replicate
  boot_fn <- function(idx) {
    sapply(threshold, function(c) {
      if (metric == "sensitivity") {
        .compute_tr_sensitivity(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          source = source[idx],
          covariates = covariates[idx, , drop = FALSE],
          threshold = c,
          treatment_level = treatment_level,
          analysis = analysis,
          estimator = estimator,
          selection_model = NULL,  # Refit in each bootstrap
          propensity_model = NULL,
          outcome_model = NULL
        )
      } else {
        .compute_tr_specificity(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          source = source[idx],
          covariates = covariates[idx, , drop = FALSE],
          threshold = c,
          treatment_level = treatment_level,
          analysis = analysis,
          estimator = estimator,
          selection_model = NULL,
          propensity_model = NULL,
          outcome_model = NULL
        )
      }
    })
  }

  # Generate bootstrap samples
  if (stratified) {
    # Stratify by source
    idx_source <- which(source == 1)
    idx_target <- which(source == 0)

    boot_indices <- lapply(seq_len(n_boot), function(b) {
      c(
        sample(idx_source, length(idx_source), replace = TRUE),
        sample(idx_target, length(idx_target), replace = TRUE)
      )
    })
  } else {
    boot_indices <- lapply(seq_len(n_boot), function(b) {
      sample(n, n, replace = TRUE)
    })
  }

  # Run bootstrap
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    boot_estimates <- parallel::mclapply(boot_indices, boot_fn, mc.cores = ncores)
  } else {
    boot_estimates <- lapply(boot_indices, boot_fn)
  }

  # Convert to matrix (n_boot x n_thresholds)
  boot_matrix <- do.call(rbind, boot_estimates)

  # Compute SE and CI for each threshold
  se <- apply(boot_matrix, 2, sd, na.rm = TRUE)
  ci_lower <- apply(boot_matrix, 2, quantile, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- apply(boot_matrix, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  list(
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    boot_estimates = boot_matrix
  )
}


# ==============================================================================
# Print Methods
# ==============================================================================

#' @export
print.tr_sensitivity <- function(x, digits = 4, ...) {
  cat("\nTransportable Sensitivity Estimate\n")
  cat("===================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Analysis:", x$analysis, "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N (source):", x$n_source, "\n")
  cat("N (target):", x$n_target, "\n\n")

  if (length(x$threshold) == 1) {
    cat("Threshold:", x$threshold, "\n")
    cat("Estimate:", round(x$estimate, digits), "\n")
    if (!is.null(x$se)) {
      cat("SE:", round(x$se, digits), "\n")
      cat(sprintf("%.0f%% CI: [%.4f, %.4f]\n",
                  x$conf_level * 100, x$ci_lower, x$ci_upper))
    }
    cat("Naive estimate:", round(x$naive_estimate, digits), "\n")
  } else {
    cat("Results by threshold:\n")
    df <- data.frame(
      Threshold = x$threshold,
      Estimate = round(x$estimate, digits),
      Naive = round(x$naive_estimate, digits)
    )
    if (!is.null(x$se)) {
      df$SE <- round(x$se, digits)
      df$CI_Lower <- round(x$ci_lower, digits)
      df$CI_Upper <- round(x$ci_upper, digits)
    }
    print(df, row.names = FALSE)
  }

  invisible(x)
}


#' @export
print.tr_specificity <- function(x, digits = 4, ...) {
  cat("\nTransportable Specificity Estimate\n")
  cat("===================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Analysis:", x$analysis, "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N (source):", x$n_source, "\n")
  cat("N (target):", x$n_target, "\n\n")

  if (length(x$threshold) == 1) {
    cat("Threshold:", x$threshold, "\n")
    cat("Estimate:", round(x$estimate, digits), "\n")
    if (!is.null(x$se)) {
      cat("SE:", round(x$se, digits), "\n")
      cat(sprintf("%.0f%% CI: [%.4f, %.4f]\n",
                  x$conf_level * 100, x$ci_lower, x$ci_upper))
    }
    cat("Naive estimate:", round(x$naive_estimate, digits), "\n")
  } else {
    cat("Results by threshold:\n")
    df <- data.frame(
      Threshold = x$threshold,
      Estimate = round(x$estimate, digits),
      Naive = round(x$naive_estimate, digits)
    )
    if (!is.null(x$se)) {
      df$SE <- round(x$se, digits)
      df$CI_Lower <- round(x$ci_lower, digits)
      df$CI_Upper <- round(x$ci_upper, digits)
    }
    print(df, row.names = FALSE)
  }

  invisible(x)
}


#' @export
print.tr_fpr <- function(x, digits = 4, ...) {
  cat("\nTransportable False Positive Rate Estimate\n")
  cat("==========================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Analysis:", x$analysis, "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N (source):", x$n_source, "\n")
  cat("N (target):", x$n_target, "\n\n")

  if (length(x$threshold) == 1) {
    cat("Threshold:", x$threshold, "\n")
    cat("Estimate:", round(x$estimate, digits), "\n")
    if (!is.null(x$se)) {
      cat("SE:", round(x$se, digits), "\n")
      cat(sprintf("%.0f%% CI: [%.4f, %.4f]\n",
                  x$conf_level * 100, x$ci_lower, x$ci_upper))
    }
    cat("Naive estimate:", round(x$naive_estimate, digits), "\n")
  } else {
    cat("Results by threshold:\n")
    df <- data.frame(
      Threshold = x$threshold,
      Estimate = round(x$estimate, digits),
      Naive = round(x$naive_estimate, digits)
    )
    if (!is.null(x$se)) {
      df$SE <- round(x$se, digits)
      df$CI_Lower <- round(x$ci_lower, digits)
      df$CI_Upper <- round(x$ci_upper, digits)
    }
    print(df, row.names = FALSE)
  }

  invisible(x)
}
