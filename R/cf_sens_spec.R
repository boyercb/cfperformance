#' Estimate Counterfactual Sensitivity
#'
#' Estimates the sensitivity (true positive rate) of a binary classifier at
#' one or more thresholds under a hypothetical intervention where treatment is
#' set to a specific level.
#'
#' @inheritParams cf_mse
#' @param threshold Numeric vector of classification thresholds. Predictions
#'   above this value are classified as positive. Can be a single value or
#'   a vector for computing sensitivity at multiple thresholds simultaneously.
#'   Default is 0.5.
#'
#' @return An object of class `c("cf_sensitivity", "cf_performance")` containing:
#'   \item{estimate}{Point estimate(s) of counterfactual sensitivity}
#'   \item{se}{Standard error(s) (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound(s)}
#'   \item{ci_upper}{Upper confidence interval bound(s)}
#'   \item{threshold}{Threshold value(s) used}
#'   \item{estimator}{Estimator used}
#'   \item{naive_estimate}{Naive sensitivity for comparison}
#'   \item{n_obs}{Number of observations}
#'   \item{treatment_level}{Counterfactual treatment level}
#'
#' @details
#' Sensitivity (also known as true positive rate or recall) is defined as:
#' \deqn{Sensitivity(c) = P(\hat{Y} > c | Y^{(a)} = 1)}
#'
#' where \eqn{Y^{(a)}} is the potential outcome under treatment level \eqn{a}.
#'
#' The function implements three estimators following Coston et al. (2020):
#'
#' **Conditional Loss (CL) Estimator**: Weights by the predicted probability of
#' being a case under the counterfactual:
#' \deqn{\hat{\psi}_{sens,cl} = \frac{\sum_i I(\hat{h}(X_i) > c) \hat{m}(X_i)}{\sum_i \hat{m}(X_i)}}
#' where \eqn{\hat{m}(X) = P(Y=1|X, A=a)}.
#'
#' **IPW Estimator**: Weights by the inverse probability of treatment:
#' \deqn{\hat{\psi}_{sens,ipw} = \frac{\sum_i I(\hat{h}(X_i) > c, Y_i=1, A_i=a) / \hat{e}(X_i)}{\sum_i I(Y_i=1, A_i=a) / \hat{e}(X_i)}}
#' where \eqn{\hat{e}(X) = P(A=a|X)}.
#'
#' **Doubly Robust (DR) Estimator**: Combines CL and IPW for protection against
#' misspecification of either model.
#'
#' @references
#' Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
#' "Counterfactual risk assessments, evaluation, and fairness."
#' *Proceedings of the 2020 Conference on Fairness, Accountability, and
#' Transparency*, 582-593.
#'
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @seealso [cf_specificity()], [cf_fpr()], [cf_auc()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 1000
#' x <- rnorm(n)
#' a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate counterfactual sensitivity at default threshold (0.5)
#' result <- cf_sensitivity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"
#' )
#' print(result)
#'
#' # Estimate at multiple thresholds
#' result_multi <- cf_sensitivity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   threshold = c(0.3, 0.5, 0.7),
#'   se_method = "none"
#' )
#' print(result_multi)
cf_sensitivity <- function(predictions,
                           outcomes,
                           treatment,
                           covariates,
                           threshold = 0.5,
                           treatment_level = 0,
                           estimator = c("dr", "cl", "ipw", "naive"),
                           propensity_model = NULL,
                           outcome_model = NULL,
                           se_method = c("none", "bootstrap", "influence"),
                           n_boot = 200,
                           conf_level = 0.95,
                           cross_fit = FALSE,
                           n_folds = 5,
                           parallel = FALSE,
                           ncores = NULL,
                           ps_trim = NULL,
                           ...) {

  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  # Influence function SE requires cross-fitting
  if (se_method == "influence" && !cross_fit) {
    stop("Influence function standard errors require cross_fit = TRUE")
  }

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  # Validate inputs
  .validate_inputs(predictions, outcomes, treatment, covariates)

  if (!all(outcomes %in% c(0, 1))) {
    stop("Sensitivity requires binary outcomes (0/1)")
  }

  n <- length(outcomes)
  n_thresholds <- length(threshold)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Cross-fitting implementation for DR estimator
  if (cross_fit && estimator == "dr") {
    # Detect if ml_learners are provided
    use_ml_propensity <- is_ml_learner(propensity_model)
    use_ml_outcome <- is_ml_learner(outcome_model)

    # Cross-fit nuisance models
    cf_nuisance <- .cross_fit_nuisance_sens_spec(
      treatment = treatment,
      outcomes = outcomes,
      covariates = covariates,
      treatment_level = treatment_level,
      K = n_folds,
      propensity_learner = if (use_ml_propensity) propensity_model else NULL,
      outcome_learner = if (use_ml_outcome) outcome_model else NULL,
      ps_trim_spec = ps_trim_spec
    )

    # Compute point estimates and optionally SE using cross-fitted nuisance
    if (se_method == "influence") {
      # Compute estimate and SE together via influence function
      results <- lapply(threshold, function(c) {
        .compute_cf_sensitivity_crossfit(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m,
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
        .compute_cf_sensitivity_crossfit(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m,
          return_se = FALSE
        )
      })
    }

    nuisance <- list(
      propensity = propensity_model,
      outcome = outcome_model,
      cross_fitted = TRUE,
      n_folds = n_folds,
      ps = cf_nuisance$ps,
      m = cf_nuisance$m
    )

    # Bootstrap SE for cross-fitted estimates (if requested)
    if (se_method == "bootstrap") {
      boot_result <- .bootstrap_cf_sens_spec_crossfit(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = threshold,
        treatment_level = treatment_level,
        metric = "sensitivity",
        n_folds = n_folds,
        n_boot = n_boot,
        conf_level = conf_level,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        parallel = parallel,
        ncores = ncores
      )
      se <- boot_result$se
      ci_lower <- boot_result$ci_lower
      ci_upper <- boot_result$ci_upper
    }

  } else {
    # Non-cross-fitted estimation (original code path)

    # Warn if cross_fit requested for non-DR estimator
    if (cross_fit && estimator != "dr") {
      warning("Cross-fitting is only implemented for estimator='dr'. Proceeding without cross-fitting.")
    }

    # Fit nuisance models if needed
    if (estimator != "naive") {
      nuisance <- .fit_nuisance_models(
        treatment = treatment,
        outcomes = outcomes,
        covariates = covariates,
        treatment_level = treatment_level,
        propensity_model = propensity_model,
        outcome_model = outcome_model
      )
    } else {
      nuisance <- list(propensity = NULL, outcome = NULL)
    }

    # Compute point estimates for all thresholds
    estimate <- sapply(threshold, function(c) {
      .compute_cf_sensitivity(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = c,
        treatment_level = treatment_level,
        estimator = estimator,
        propensity_model = nuisance$propensity,
        outcome_model = nuisance$outcome
      )
    })

    # Bootstrap standard errors
    if (se_method == "bootstrap") {
      boot_result <- .bootstrap_cf_sens_spec(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = threshold,
        treatment_level = treatment_level,
        estimator = estimator,
        metric = "sensitivity",
        n_boot = n_boot,
        conf_level = conf_level,
        parallel = parallel,
        ncores = ncores
      )
      se <- boot_result$se
      ci_lower <- boot_result$ci_lower
      ci_upper <- boot_result$ci_upper
    }
  }

  # Naive estimate (always compute for comparison)
  naive_estimate <- sapply(threshold, function(c) {
    .compute_cf_sens_spec_naive(predictions, outcomes, c, "sensitivity")
  })

  # Build result object
  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    threshold = threshold,
    estimator = estimator,
    naive_estimate = naive_estimate,
    n_obs = n,
    treatment_level = treatment_level,
    conf_level = conf_level,
    metric = "sensitivity"
  )

  class(result) <- c("cf_sensitivity", "cf_performance")
  return(result)
}


#' Estimate Counterfactual Specificity
#'
#' Estimates the specificity (true negative rate) of a binary classifier at
#' one or more thresholds under a hypothetical intervention where treatment is
#' set to a specific level.
#'
#' @inheritParams cf_sensitivity
#'
#' @return An object of class `c("cf_specificity", "cf_performance")` containing:
#'   \item{estimate}{Point estimate(s) of counterfactual specificity}
#'   \item{se}{Standard error(s) (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound(s)}
#'   \item{ci_upper}{Upper confidence interval bound(s)}
#'   \item{threshold}{Threshold value(s) used}
#'   \item{estimator}{Estimator used}
#'   \item{naive_estimate}{Naive specificity for comparison}
#'   \item{n_obs}{Number of observations}
#'   \item{treatment_level}{Counterfactual treatment level}
#'
#' @details
#' Specificity (also known as true negative rate) is defined as:
#' \deqn{Specificity(c) = P(\hat{Y} \leq c | Y^{(a)} = 0)}
#'
#' where \eqn{Y^{(a)}} is the potential outcome under treatment level \eqn{a}.
#' The estimators mirror those for sensitivity (see [cf_sensitivity()]).
#'
#' @references
#' Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
#' "Counterfactual risk assessments, evaluation, and fairness."
#' *Proceedings of the 2020 Conference on Fairness, Accountability, and
#' Transparency*, 582-593.
#'
#' @seealso [cf_sensitivity()], [cf_fpr()], [cf_auc()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 1000
#' x <- rnorm(n)
#' a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate counterfactual specificity at default threshold (0.5)
#' result <- cf_specificity(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   estimator = "dr",
#'   se_method = "none"
#' )
#' print(result)
cf_specificity <- function(predictions,
                           outcomes,
                           treatment,
                           covariates,
                           threshold = 0.5,
                           treatment_level = 0,
                           estimator = c("dr", "cl", "ipw", "naive"),
                           propensity_model = NULL,
                           outcome_model = NULL,
                           se_method = c("none", "bootstrap", "influence"),
                           n_boot = 200,
                           conf_level = 0.95,
                           cross_fit = FALSE,
                           n_folds = 5,
                           parallel = FALSE,
                           ncores = NULL,
                           ps_trim = NULL,
                           ...) {

  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  # Influence function SE requires cross-fitting
  if (se_method == "influence" && !cross_fit) {
    stop("Influence function standard errors require cross_fit = TRUE")
  }

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  # Validate inputs
  .validate_inputs(predictions, outcomes, treatment, covariates)

  if (!all(outcomes %in% c(0, 1))) {
    stop("Specificity requires binary outcomes (0/1)")
  }

  n <- length(outcomes)
  n_thresholds <- length(threshold)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Cross-fitting implementation for DR estimator
  if (cross_fit && estimator == "dr") {
    # Detect if ml_learners are provided
    use_ml_propensity <- is_ml_learner(propensity_model)
    use_ml_outcome <- is_ml_learner(outcome_model)

    # Cross-fit nuisance models
    cf_nuisance <- .cross_fit_nuisance_sens_spec(
      treatment = treatment,
      outcomes = outcomes,
      covariates = covariates,
      treatment_level = treatment_level,
      K = n_folds,
      propensity_learner = if (use_ml_propensity) propensity_model else NULL,
      outcome_learner = if (use_ml_outcome) outcome_model else NULL,
      ps_trim_spec = ps_trim_spec
    )

    # Compute point estimates and optionally SE using cross-fitted nuisance
    if (se_method == "influence") {
      # Compute estimate and SE together via influence function
      results <- lapply(threshold, function(c) {
        .compute_cf_specificity_crossfit(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m,
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
        .compute_cf_specificity_crossfit(
          predictions = predictions,
          outcomes = outcomes,
          treatment = treatment,
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m,
          return_se = FALSE
        )
      })
    }

    nuisance <- list(
      propensity = propensity_model,
      outcome = outcome_model,
      cross_fitted = TRUE,
      n_folds = n_folds,
      ps = cf_nuisance$ps,
      m = cf_nuisance$m
    )

    # Bootstrap SE for cross-fitted estimates (if requested)
    if (se_method == "bootstrap") {
      boot_result <- .bootstrap_cf_sens_spec_crossfit(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = threshold,
        treatment_level = treatment_level,
        metric = "specificity",
        n_folds = n_folds,
        n_boot = n_boot,
        conf_level = conf_level,
        propensity_learner = if (use_ml_propensity) propensity_model else NULL,
        outcome_learner = if (use_ml_outcome) outcome_model else NULL,
        parallel = parallel,
        ncores = ncores
      )
      se <- boot_result$se
      ci_lower <- boot_result$ci_lower
      ci_upper <- boot_result$ci_upper
    }

  } else {
    # Non-cross-fitted estimation (original code path)

    # Warn if cross_fit requested for non-DR estimator
    if (cross_fit && estimator != "dr") {
      warning("Cross-fitting is only implemented for estimator='dr'. Proceeding without cross-fitting.")
    }

    # Fit nuisance models if needed
    if (estimator != "naive") {
      nuisance <- .fit_nuisance_models(
        treatment = treatment,
        outcomes = outcomes,
        covariates = covariates,
        treatment_level = treatment_level,
        propensity_model = propensity_model,
        outcome_model = outcome_model
      )
    } else {
      nuisance <- list(propensity = NULL, outcome = NULL)
    }

    # Compute point estimates for all thresholds
    estimate <- sapply(threshold, function(c) {
      .compute_cf_specificity(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = c,
        treatment_level = treatment_level,
        estimator = estimator,
        propensity_model = nuisance$propensity,
        outcome_model = nuisance$outcome
      )
    })

    # Bootstrap standard errors
    if (se_method == "bootstrap") {
      boot_result <- .bootstrap_cf_sens_spec(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        covariates = covariates,
        threshold = threshold,
        treatment_level = treatment_level,
        estimator = estimator,
        metric = "specificity",
        n_boot = n_boot,
        conf_level = conf_level,
        parallel = parallel,
        ncores = ncores
      )
      se <- boot_result$se
      ci_lower <- boot_result$ci_lower
      ci_upper <- boot_result$ci_upper
    }
  }

  # Naive estimate (always compute for comparison)
  naive_estimate <- sapply(threshold, function(c) {
    .compute_cf_sens_spec_naive(predictions, outcomes, c, "specificity")
  })

  # Build result object
  result <- list(
    estimate = estimate,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    threshold = threshold,
    estimator = estimator,
    naive_estimate = naive_estimate,
    n_obs = n,
    treatment_level = treatment_level,
    conf_level = conf_level,
    metric = "specificity"
  )

  class(result) <- c("cf_specificity", "cf_performance")
  return(result)
}


#' Estimate Counterfactual False Positive Rate
#'
#' Estimates the false positive rate (FPR) of a binary classifier, which is
#' the complement of specificity (FPR = 1 - specificity).
#'
#' @inheritParams cf_specificity
#'
#' @return An object of class `c("cf_fpr", "cf_performance")` containing:
#'   \item{estimate}{Point estimate(s) of counterfactual FPR}
#'   \item{se}{Standard error(s) (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound(s)}
#'   \item{ci_upper}{Upper confidence interval bound(s)}
#'   \item{threshold}{Threshold value(s) used}
#'   \item{estimator}{Estimator used}
#'   \item{naive_estimate}{Naive FPR for comparison}
#'
#' @details
#' False positive rate is defined as:
#' \deqn{FPR(c) = P(\hat{Y} > c | Y^{(a)} = 0) = 1 - Specificity(c)}
#'
#' This function is provided as a convenience for ROC curve construction,
#' where FPR is typically plotted on the x-axis.
#'
#' @seealso [cf_sensitivity()], [cf_specificity()], [cf_tpr()]
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' result <- cf_fpr(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   se_method = "none"
#' )
#' print(result)
cf_fpr <- function(predictions,
                   outcomes,
                   treatment,
                   covariates,
                   threshold = 0.5,
                   treatment_level = 0,
                   estimator = c("dr", "cl", "ipw", "naive"),
                   propensity_model = NULL,
                   outcome_model = NULL,
                   se_method = c("none", "bootstrap", "influence"),
                   n_boot = 200,
                   conf_level = 0.95,
                   cross_fit = FALSE,
                   n_folds = 5,
                   parallel = FALSE,
                   ncores = NULL,
                   ps_trim = NULL,
                   ...) {

  # Call specificity and transform
  spec_result <- cf_specificity(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    covariates = covariates,
    threshold = threshold,
    treatment_level = treatment_level,
    estimator = estimator,
    propensity_model = propensity_model,
    outcome_model = outcome_model,
    se_method = se_method,
    n_boot = n_boot,
    conf_level = conf_level,
    cross_fit = cross_fit,
    n_folds = n_folds,
    parallel = parallel,
    ncores = ncores,
    ps_trim = ps_trim,
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
  class(spec_result) <- c("cf_fpr", "cf_performance")

  return(spec_result)
}


#' @rdname cf_sensitivity
#' @export
cf_tpr <- cf_sensitivity


#' @rdname cf_specificity
#' @export
cf_tnr <- cf_specificity


# ==============================================================================
# Cross-fitting for counterfactual sensitivity/specificity
# ==============================================================================

#' Cross-fit nuisance models for counterfactual sens/spec
#'
#' Implements K-fold cross-fitting for nuisance model estimation.
#'
#' @inheritParams cf_sensitivity
#' @param K Number of folds for cross-fitting.
#' @param propensity_learner Optional ml_learner for propensity model.
#' @param outcome_learner Optional ml_learner for outcome model.
#'
#' @return List containing cross-fitted nuisance function predictions.
#'
#' @keywords internal
.cross_fit_nuisance_sens_spec <- function(treatment, outcomes, covariates,
                                           treatment_level, K = 5,
                                           propensity_learner = NULL,
                                           outcome_learner = NULL,
                                           ps_trim_spec = NULL) {

  n <- length(outcomes)

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Create fold assignments
  folds <- sample(rep(1:K, length.out = n))

  # Initialize output vectors
  ps_cf <- numeric(n)    # Cross-fitted propensity scores P(A=a|X)
  m_cf <- numeric(n)     # Cross-fitted outcome probabilities P(Y=1|X,A=a)

  for (k in 1:K) {
    train_idx <- which(folds != k)
    val_idx <- which(folds == k)

    # --- Propensity model: P(A=1|X) ---
    ps_data <- data.frame(A = treatment[train_idx], covariates[train_idx, , drop = FALSE])

    if (!is.null(propensity_learner) && is_ml_learner(propensity_learner)) {
      ps_model <- .fit_ml_learner(propensity_learner, A ~ .,
                                   data = ps_data, family = "binomial")
      ps_pred <- .predict_ml_learner(ps_model, covariates[val_idx, , drop = FALSE])
    } else {
      ps_model <- glm(A ~ ., data = ps_data, family = binomial())
      ps_pred <- predict(ps_model, newdata = covariates[val_idx, , drop = FALSE],
                         type = "response")
    }

    # Convert to P(A=treatment_level|X)
    if (treatment_level == 0) {
      ps_pred <- 1 - ps_pred
    }
    ps_cf[val_idx] <- .trim_propensity(ps_pred, ps_trim_spec$method, ps_trim_spec$bounds)

    # --- Outcome model: E[Y | X, A=a] ---
    train_a <- train_idx[treatment[train_idx] == treatment_level]
    om_data <- data.frame(Y = outcomes[train_a], covariates[train_a, , drop = FALSE])

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
    m_cf[val_idx] <- pmin(pmax(m_pred, 0), 1)
  }

  list(
    ps = ps_cf,
    m = m_cf,
    folds = folds
  )
}


#' Compute DR sensitivity with cross-fitted nuisance
#'
#' Computes the doubly robust sensitivity estimator using cross-fitted
#' nuisance functions and returns both point estimate and influence
#' function-based standard error.
#'
#' @inheritParams cf_sensitivity
#' @param ps Cross-fitted propensity scores.
#' @param m_hat Cross-fitted outcome probabilities.
#' @param return_se Logical; if TRUE, returns list with estimate and SE.
#'
#' @return If return_se=FALSE, numeric estimate. If TRUE, list with estimate and se.
#' @keywords internal
.compute_cf_sensitivity_crossfit <- function(predictions, outcomes, treatment,
                                              threshold, treatment_level,
                                              ps, m_hat, return_se = FALSE) {

  n <- length(outcomes)

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_pos <- as.numeric(predictions > threshold)
  I_case <- as.numeric(outcomes == 1)

  # IPW weights
  ipw_weight <- I_a / ps

  # Numerator: DR for P(pred > c, Y=1) = E[I(pred>c) * Y^a]
  # DR estimator: E[m(X)*I(pred>c)] + E[(I_a/ps)*(I(pred>c)*Y - m(X)*I(pred>c))]
  num_cl <- mean(I_pos * m_hat)
  num_ipw <- mean(ipw_weight * I_pos * I_case)
  num_aug <- mean(ipw_weight * I_pos * m_hat)
  mu_1 <- num_cl + num_ipw - num_aug

  # Denominator: DR for P(Y=1) = E[Y^a]
  denom_cl <- mean(m_hat)
  denom_ipw <- mean(ipw_weight * I_case)
  denom_aug <- mean(ipw_weight * m_hat)
  mu_0 <- denom_cl + denom_ipw - denom_aug

  if (mu_0 <= 0) {
    if (return_se) return(list(estimate = NA_real_, se = NA_real_))
    return(NA_real_)
  }

  estimate <- mu_1 / mu_0

  if (!return_se) return(estimate)

  # Influence function for ratio estimator: phi = (phi_1 - psi * phi_0) / mu_0
  # phi_1: influence function for numerator (DR for joint probability)
  # phi_0: influence function for denominator (DR for marginal)

  # phi_1_i = I(pred>c)*m(X) + (I_a/ps)*(I(pred>c)*Y - I(pred>c)*m(X)) - mu_1
  phi_1 <- I_pos * m_hat + ipw_weight * (I_pos * I_case - I_pos * m_hat) - mu_1

  # phi_0_i = m(X) + (I_a/ps)*(Y - m(X)) - mu_0
  phi_0 <- m_hat + ipw_weight * (I_case - m_hat) - mu_0

  # Influence function for ratio (delta method)
  phi <- (phi_1 - estimate * phi_0) / mu_0

  se <- sqrt(var(phi) / n)

  list(estimate = estimate, se = se)
}


#' Compute DR specificity with cross-fitted nuisance
#'
#' Computes the doubly robust specificity estimator using cross-fitted
#' nuisance functions and returns both point estimate and influence
#' function-based standard error.
#'
#' @inheritParams cf_specificity
#' @param ps Cross-fitted propensity scores.
#' @param m_hat Cross-fitted outcome probabilities.
#' @param return_se Logical; if TRUE, returns list with estimate and SE.
#'
#' @return If return_se=FALSE, numeric estimate. If TRUE, list with estimate and se.
#' @keywords internal
.compute_cf_specificity_crossfit <- function(predictions, outcomes, treatment,
                                              threshold, treatment_level,
                                              ps, m_hat, return_se = FALSE) {

  n <- length(outcomes)

  # For specificity, we need P(Y=0|X) = 1 - m_hat
  m_hat_neg <- 1 - m_hat

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_neg <- as.numeric(predictions <= threshold)
  I_control <- as.numeric(outcomes == 0)

  # IPW weights
  ipw_weight <- I_a / ps

  # Numerator: DR for P(pred <= c, Y=0) = E[I(pred<=c) * (1-Y^a)]
  num_cl <- mean(I_neg * m_hat_neg)
  num_ipw <- mean(ipw_weight * I_neg * I_control)
  num_aug <- mean(ipw_weight * I_neg * m_hat_neg)
  mu_1 <- num_cl + num_ipw - num_aug

  # Denominator: DR for P(Y=0) = E[1-Y^a]
  denom_cl <- mean(m_hat_neg)
  denom_ipw <- mean(ipw_weight * I_control)
  denom_aug <- mean(ipw_weight * m_hat_neg)
  mu_0 <- denom_cl + denom_ipw - denom_aug

  if (mu_0 <= 0) {
    if (return_se) return(list(estimate = NA_real_, se = NA_real_))
    return(NA_real_)
  }

  estimate <- mu_1 / mu_0

  if (!return_se) return(estimate)

  # Influence function for ratio estimator
  # phi_1: influence function for numerator
  phi_1 <- I_neg * m_hat_neg + ipw_weight * (I_neg * I_control - I_neg * m_hat_neg) - mu_1

  # phi_0: influence function for denominator
  phi_0 <- m_hat_neg + ipw_weight * (I_control - m_hat_neg) - mu_0

  # Influence function for ratio (delta method)
  phi <- (phi_1 - estimate * phi_0) / mu_0

  se <- sqrt(var(phi) / n)

  list(estimate = estimate, se = se)
}


#' Bootstrap for cross-fitted sens/spec
#' @keywords internal
.bootstrap_cf_sens_spec_crossfit <- function(predictions, outcomes, treatment,
                                              covariates, threshold,
                                              treatment_level, metric,
                                              n_folds, n_boot, conf_level,
                                              propensity_learner = NULL,
                                              outcome_learner = NULL,
                                              parallel, ncores,
                                              ps_trim_spec = NULL) {

  n <- length(outcomes)
  n_thresholds <- length(threshold)
  alpha <- 1 - conf_level

  boot_fn <- function(idx) {
    # Cross-fit on bootstrap sample
    cf_nuisance <- tryCatch({
      .cross_fit_nuisance_sens_spec(
        treatment = treatment[idx],
        outcomes = outcomes[idx],
        covariates = covariates[idx, , drop = FALSE],
        treatment_level = treatment_level,
        K = n_folds,
        propensity_learner = propensity_learner,
        outcome_learner = outcome_learner,
        ps_trim_spec = ps_trim_spec
      )
    }, error = function(e) NULL)

    if (is.null(cf_nuisance)) return(rep(NA_real_, n_thresholds))

    sapply(threshold, function(c) {
      if (metric == "sensitivity") {
        .compute_cf_sensitivity_crossfit(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m
        )
      } else {
        .compute_cf_specificity_crossfit(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          threshold = c,
          treatment_level = treatment_level,
          ps = cf_nuisance$ps,
          m_hat = cf_nuisance$m
        )
      }
    })
  }

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    boot_results <- parallel::mclapply(
      1:n_boot,
      function(b) {
        set.seed(b)
        idx <- sample(n, n, replace = TRUE)
        boot_fn(idx)
      },
      mc.cores = ncores
    )
    boot_matrix <- do.call(rbind, boot_results)
  } else {
    boot_matrix <- matrix(NA_real_, nrow = n_boot, ncol = n_thresholds)
    for (b in 1:n_boot) {
      idx <- sample(n, n, replace = TRUE)
      boot_matrix[b, ] <- boot_fn(idx)
    }
  }

  # Compute SEs and CIs
  se <- apply(boot_matrix, 2, sd, na.rm = TRUE)
  ci_lower <- apply(boot_matrix, 2, quantile, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- apply(boot_matrix, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  list(
    se = se,
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper)
  )
}


# ==============================================================================
# Internal Functions: Computation
# ==============================================================================

#' Compute counterfactual sensitivity
#' @noRd
.compute_cf_sensitivity <- function(predictions, outcomes, treatment, covariates,
                                     threshold, treatment_level, estimator,
                                     propensity_model, outcome_model) {

  if (estimator == "naive") {
    return(.compute_cf_sens_spec_naive(predictions, outcomes, threshold, "sensitivity"))
  }

  # Refit nuisance models if NULL (e.g., in bootstrap)
  if (is.null(propensity_model) || is.null(outcome_model)) {
    nuisance <- .fit_nuisance_models(
      treatment = treatment,
      outcomes = outcomes,
      covariates = covariates,
      treatment_level = treatment_level,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
    propensity_model <- nuisance$propensity
    outcome_model <- nuisance$outcome
  }

  # Convert covariates to data frame
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  n <- length(outcomes)

  # Get propensity scores P(A=a|X)
  ps <- predict(propensity_model, newdata = covariates, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }
  # Clip extreme propensity scores
  ps <- pmax(ps, 0.01)
  ps <- pmin(ps, 0.99)

  # Get outcome probabilities m_hat = P(Y=1|X, A=a)
  m_hat <- predict(outcome_model, newdata = covariates, type = "response")

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_pos <- as.numeric(predictions > threshold)
  I_case <- as.numeric(outcomes == 1)

  if (estimator == "cl") {
    # Conditional loss estimator
    # Sens = E[I(pred > c) * m(X)] / E[m(X)]
    num <- mean(I_pos * m_hat)
    denom <- mean(m_hat)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "ipw") {
    # IPW estimator
    # Sens = E[I(pred > c) * I(Y=1) * I(A=a) / e(X)] / E[I(Y=1) * I(A=a) / e(X)]
    ipw_weight <- I_a / ps

    num <- mean(I_pos * I_case * ipw_weight)
    denom <- mean(I_case * ipw_weight)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "dr") {
    # Doubly robust estimator
    ipw_weight <- I_a / ps

    # Numerator: DR for P(pred > c, Y=1)
    num_cl <- mean(I_pos * m_hat)
    num_ipw <- mean(ipw_weight * I_pos * I_case)
    num_aug <- mean(ipw_weight * I_pos * m_hat)

    # Denominator: DR for P(Y=1)
    denom_cl <- mean(m_hat)
    denom_ipw <- mean(ipw_weight * I_case)
    denom_aug <- mean(ipw_weight * m_hat)

    num <- num_cl + num_ipw - num_aug
    denom <- denom_cl + denom_ipw - denom_aug

    if (denom <= 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


#' Compute counterfactual specificity
#' @noRd
.compute_cf_specificity <- function(predictions, outcomes, treatment, covariates,
                                     threshold, treatment_level, estimator,
                                     propensity_model, outcome_model) {

  if (estimator == "naive") {
    return(.compute_cf_sens_spec_naive(predictions, outcomes, threshold, "specificity"))
  }

  # Refit nuisance models if NULL (e.g., in bootstrap)
  if (is.null(propensity_model) || is.null(outcome_model)) {
    nuisance <- .fit_nuisance_models(
      treatment = treatment,
      outcomes = outcomes,
      covariates = covariates,
      treatment_level = treatment_level,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
    propensity_model <- nuisance$propensity
    outcome_model <- nuisance$outcome
  }

  # Convert covariates to data frame
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  n <- length(outcomes)

  # Get propensity scores P(A=a|X)
  ps <- predict(propensity_model, newdata = covariates, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }
  ps <- pmax(ps, 0.01)
  ps <- pmin(ps, 0.99)

  # Get outcome probabilities m_hat = P(Y=1|X, A=a)
  m_hat <- predict(outcome_model, newdata = covariates, type = "response")
  # For specificity, we need P(Y=0|X) = 1 - m_hat
  m_hat_neg <- 1 - m_hat

  # Indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_neg <- as.numeric(predictions <= threshold)
  I_control <- as.numeric(outcomes == 0)

  if (estimator == "cl") {
    # Conditional loss estimator
    # Spec = E[I(pred <= c) * (1 - m(X))] / E[1 - m(X)]
    num <- mean(I_neg * m_hat_neg)
    denom <- mean(m_hat_neg)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "ipw") {
    # IPW estimator
    ipw_weight <- I_a / ps

    num <- mean(I_neg * I_control * ipw_weight)
    denom <- mean(I_control * ipw_weight)

    if (denom == 0) return(NA_real_)
    return(num / denom)

  } else if (estimator == "dr") {
    # Doubly robust estimator
    ipw_weight <- I_a / ps

    # Numerator: DR for P(pred <= c, Y=0)
    num_cl <- mean(I_neg * m_hat_neg)
    num_ipw <- mean(ipw_weight * I_neg * I_control)
    num_aug <- mean(ipw_weight * I_neg * m_hat_neg)

    # Denominator: DR for P(Y=0)
    denom_cl <- mean(m_hat_neg)
    denom_ipw <- mean(ipw_weight * I_control)
    denom_aug <- mean(ipw_weight * m_hat_neg)

    num <- num_cl + num_ipw - num_aug
    denom <- denom_cl + denom_ipw - denom_aug

    if (denom <= 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


#' Compute naive sensitivity/specificity
#' @noRd
.compute_cf_sens_spec_naive <- function(predictions, outcomes, threshold, metric) {
  if (metric == "sensitivity") {
    cases <- outcomes == 1
    if (sum(cases) == 0) return(NA_real_)
    return(mean(predictions[cases] > threshold))
  } else {
    controls <- outcomes == 0
    if (sum(controls) == 0) return(NA_real_)
    return(mean(predictions[controls] <= threshold))
  }
}


# ==============================================================================
# Internal Functions: Bootstrap
# ==============================================================================

#' Bootstrap for counterfactual sensitivity/specificity
#' @noRd
.bootstrap_cf_sens_spec <- function(predictions, outcomes, treatment, covariates,
                                     threshold, treatment_level, estimator,
                                     metric, n_boot, conf_level, parallel, ncores) {

  n <- length(outcomes)
  n_thresholds <- length(threshold)
  alpha <- 1 - conf_level

  # Bootstrap function for single replicate
  boot_fn <- function(idx) {
    sapply(threshold, function(c) {
      if (metric == "sensitivity") {
        .compute_cf_sensitivity(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          covariates = covariates[idx, , drop = FALSE],
          threshold = c,
          treatment_level = treatment_level,
          estimator = estimator,
          propensity_model = NULL,  # Refit in each bootstrap
          outcome_model = NULL
        )
      } else {
        .compute_cf_specificity(
          predictions = predictions[idx],
          outcomes = outcomes[idx],
          treatment = treatment[idx],
          covariates = covariates[idx, , drop = FALSE],
          threshold = c,
          treatment_level = treatment_level,
          estimator = estimator,
          propensity_model = NULL,
          outcome_model = NULL
        )
      }
    })
  }

  # Generate bootstrap samples
  boot_indices <- lapply(seq_len(n_boot), function(b) {
    sample(n, n, replace = TRUE)
  })

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
print.cf_sensitivity <- function(x, digits = 4, ...) {
  cat("\nCounterfactual Sensitivity Estimate\n")
  cat("====================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N:", x$n_obs, "\n\n")

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

  cat("\n")
  invisible(x)
}


#' @export
print.cf_specificity <- function(x, digits = 4, ...) {
  cat("\nCounterfactual Specificity Estimate\n")
  cat("====================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N:", x$n_obs, "\n\n")

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

  cat("\n")
  invisible(x)
}


#' @export
print.cf_fpr <- function(x, digits = 4, ...) {
  cat("\nCounterfactual False Positive Rate Estimate\n")
  cat("============================================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N:", x$n_obs, "\n\n")

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

  cat("\n")
  invisible(x)
}
