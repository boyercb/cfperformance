#' Estimate Transportable Area Under the ROC Curve in the Target Population
#'
#' Estimates the area under the receiver operating characteristic curve (AUC)
#' of a prediction model in a target population using data transported from
#' a source population. Supports both **counterfactual** (under hypothetical
#' intervention) and **factual** (observational) prediction model
#' transportability.
#'
#' @inheritParams tr_mse
#'
#' @return An object of class `c("tr_auc", "tr_performance")` containing:
#'   \item{estimate}{Point estimate of transportable AUC}
#'   \item{se}{Standard error (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound}
#'   \item{ci_upper}{Upper confidence interval bound}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{naive_estimate}{Naive AUC for comparison}
#'   \item{n_target}{Number of target observations}
#'   \item{n_source}{Number of source observations}
#'   \item{treatment_level}{Treatment level (NULL for factual mode)}
#'
#' @details
#' This function implements estimators for transporting prediction model
#' AUC from a source population to a target population. It supports two modes:
#'
#' ## Counterfactual Mode (treatment provided)
#' When `treatment` is specified, estimates the counterfactual AUC under a
#' hypothetical intervention. The AUC represents the probability that a randomly
#' selected case (Y^a = 1) has a higher predicted risk than a randomly selected
#' non-case (Y^a = 0) in the target population. This requires:
#' - Selection model: P(S=0|X)
#' - Propensity model in source: P(A=a|X, S=1)
#' - Outcome model: E[Y|X, A, S=1]
#'
#' ## Factual Mode (treatment = NULL)
#' When `treatment` is `NULL`, estimates the AUC of observed outcomes in the
#' target population. This is appropriate for factual prediction model
#' transportability without causal interpretation. Only requires:
#' - Selection model: P(S=0|X)
#' - Outcome model: E[Y|X, S=1]
#'
#' ## Analysis Types
#' **Transportability Analysis** (`analysis = "transport"`): Uses outcome data
#' from the source population to estimate AUC in the target population.
#'
#' **Joint Analysis** (`analysis = "joint"`): Pools source and target data to
#' estimate AUC. More efficient when both populations have outcome data.
#'
#' For observational analysis (single population), use [cf_auc()] instead.
#'
#' @references
#' Steingrimsson, J. A., et al. (2023). "Transporting a Prediction Model for
#' Use in a New Target Population." *American Journal of Epidemiology*,
#' 192(2), 296-304. \doi{10.1093/aje/kwac128}
#'
#' Li, S., et al. (2023). "Efficient estimation of the expected prediction
#' error under covariate shift." *Biometrics*, 79(1), 295-307.
#' \doi{10.1111/biom.13583}
#'
#' Voter, S. R., et al. (2025). "Transportability of machine learning-based
#' counterfactual prediction models with application to CASS."
#' *Diagnostic and Prognostic Research*, 9(4).
#' \doi{10.1186/s41512-025-00201-y}
#'
#' Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
#' "Estimating the area under the ROC curve when transporting a prediction
#' model to a target population." *Biometrics*.
#'
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @seealso [cf_auc()], [tr_mse()]
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
#' # Predictions from some model
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate transportable AUC
#' result <- tr_auc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   treatment_level = 0,
#'   analysis = "transport",
#'   estimator = "dr",
#'   se_method = "none"  # Skip SE for speed
#' )
#' print(result)
#'
#' @importFrom stats as.formula
tr_auc <- function(predictions,
                   outcomes,
                   treatment = NULL,
                   source,
                   covariates,
                   treatment_level = NULL,
                   analysis = c("transport", "joint"),
                   estimator = c("dr", "om", "ipw", "naive"),
                   selection_model = NULL,
                   propensity_model = NULL,
                   outcome_model = NULL,
                   se_method = c("bootstrap", "influence", "none"),
                   n_boot = 500,
                   conf_level = 0.95,
                   boot_ci_type = c("percentile", "normal", "basic"),
                   stratified_boot = TRUE,
                   cross_fit = FALSE,
                   n_folds = 5,
                   parallel = FALSE,
                   ncores = NULL,
                   ps_trim = NULL,
                   ...) {

  # Input validation
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)
  boot_ci_type <- match.arg(boot_ci_type)

  # Parse propensity score trimming specification
  ps_trim_spec <- .parse_ps_trim(ps_trim)

  .validate_transport_inputs(predictions, outcomes, treatment, source, covariates, treatment_level)
  
  # Determine mode: factual (no treatment) or counterfactual
  factual_mode <- is.null(treatment)
  
  # Default treatment_level to 1 for backward compatibility if treatment is provided
  if (!factual_mode && is.null(treatment_level)) {
    treatment_level <- 1
  }

  # Check binary outcomes
  if (!all(outcomes %in% c(0, 1))) {
    stop("AUC requires binary outcomes (0/1)")
  }

  n <- length(outcomes)
  n_target <- sum(source == 0)
  n_source <- sum(source == 1)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Cross-fitting implementation for DR estimator
  if (cross_fit && estimator == "dr") {
    # Determine if we need SE from cross-fit function
    compute_cf_se <- (se_method == "influence")

    # Use cross-fitting for nuisance model estimation
    cf_result <- .compute_tr_auc_crossfit(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      analysis = analysis,
      K = n_folds,
      selection_learner = selection_model,
      propensity_learner = propensity_model,
      outcome_learner = outcome_model,
      parallel = parallel,
      return_se = compute_cf_se,
      ps_trim_spec = ps_trim_spec,
      ...
    )

    estimate <- cf_result$estimate

    # Store nuisance info for result object
    nuisance <- list(
      selection = selection_model,
      propensity = propensity_model,
      outcome = outcome_model,
      cross_fitted = TRUE,
      n_folds = n_folds,
      pi_s0 = cf_result$pi_s0,
      ps = cf_result$ps,
      q = cf_result$q
    )

    # Compute naive estimate for comparison
    naive_estimate <- .compute_tr_auc_naive(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      treatment_level = treatment_level,
      analysis = analysis
    )

    # Compute standard errors
    if (se_method == "influence") {
      # Use influence function SE from cross-fitting (Li et al. 2022)
      se <- cf_result$se
      z_alpha <- qnorm(1 - (1 - conf_level) / 2)
      ci_lower <- estimate - z_alpha * se
      ci_upper <- estimate + z_alpha * se
    } else if (se_method == "bootstrap") {
      boot_result <- .bootstrap_tr_auc_crossfit(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        n_folds = n_folds,
        n_boot = n_boot,
        conf_level = conf_level,
        boot_ci_type = boot_ci_type,
        stratified = stratified_boot,
        selection_learner = selection_model,
        propensity_learner = propensity_model,
        outcome_learner = outcome_model,
        parallel = parallel,
        ncores = ncores,
        ps_trim_spec = ps_trim_spec,
        point_estimate = estimate,
        ...
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

    # Fit nuisance models if not provided
    if (estimator != "naive") {
      nuisance <- .fit_transport_nuisance_auc(
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
    } else {
      nuisance <- list(selection = NULL, propensity = NULL, outcome = NULL)
    }

    # Compute point estimate
    estimate <- .compute_tr_auc(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = estimator,
      selection_model = nuisance$selection,
      propensity_model = nuisance$propensity,
      outcome_model = nuisance$outcome,
      ps_trim_spec = ps_trim_spec
    )

    # Compute naive estimate for comparison
    naive_estimate <- .compute_tr_auc_naive(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      treatment_level = treatment_level,
      analysis = analysis
    )

    # Compute standard errors
    if (se_method == "bootstrap") {
      boot_result <- .bootstrap_tr_auc(
        predictions = predictions,
        outcomes = outcomes,
        treatment = treatment,
        source = source,
        covariates = covariates,
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = estimator,
        n_boot = n_boot,
        conf_level = conf_level,
        boot_ci_type = boot_ci_type,
        stratified = stratified_boot,
        parallel = parallel,
        ncores = ncores,
        ps_trim_spec = ps_trim_spec,
        selection_model = nuisance$selection,
        propensity_model = nuisance$propensity,
        outcome_model = nuisance$outcome,
        point_estimate = estimate,
        ...
      )
      se <- boot_result$se
      ci_lower <- boot_result$ci_lower
      ci_upper <- boot_result$ci_upper
    } else if (se_method == "influence") {
      # Placeholder for influence function SE
      warning("Influence function SE not yet implemented for tr_auc; use bootstrap")
      se <- NA_real_
      ci_lower <- NA_real_
      ci_upper <- NA_real_
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
    analysis = analysis,
    metric = "auc",
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

  class(result) <- c("tr_auc", "tr_performance")
  return(result)
}


# ==============================================================================
# Cross-fitting for transportability AUC
# ==============================================================================

#' Cross-fit nuisance models for transportability AUC
#'
#' Implements K-fold cross-fitting for nuisance model estimation in
#' transportability AUC estimation.
#'
#' @inheritParams tr_auc
#' @param K Number of folds for cross-fitting.
#' @param selection_learner Optional ml_learner for selection model.
#' @param propensity_learner Optional ml_learner for propensity model.
#' @param outcome_learner Optional ml_learner for outcome model.
#' @param parallel Logical for parallel processing.
#' @param ... Additional arguments.
#'
#' @return List containing cross-fitted nuisance function predictions.
#'
#' @keywords internal
.cross_fit_transport_nuisance_auc <- function(treatment, outcomes, source, covariates,
                                               treatment_level, analysis,
                                               K = 5,
                                               selection_learner = NULL,
                                               propensity_learner = NULL,
                                               outcome_learner = NULL,
                                               parallel = FALSE,
                                               ps_trim_spec = NULL,
                                               ...) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)

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
  q_cf <- numeric(n)      # Cross-fitted outcome probabilities

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

    # --- Outcome model: E[Y | X, A=a, S] ---
    if (analysis == "transport") {
      # Model E[Y | X, A=a, S=1] using source data with treatment_level
      train_s1_a <- train_idx[source[train_idx] == 1 & treatment[train_idx] == treatment_level]
      om_data <- data.frame(Y = outcomes[train_s1_a], covariates[train_s1_a, , drop = FALSE])
    } else {
      # Model E[Y | X, A=a] using all data with treatment_level
      train_a <- train_idx[treatment[train_idx] == treatment_level]
      om_data <- data.frame(Y = outcomes[train_a], covariates[train_a, , drop = FALSE])
    }

    if (!is.null(outcome_learner) && is_ml_learner(outcome_learner)) {
      om_model <- .fit_ml_learner(outcome_learner, Y ~ .,
                                   data = om_data, family = "binomial")
      q_pred <- .predict_ml_learner(om_model, covariates[val_idx, , drop = FALSE])
    } else {
      om_model <- glm(Y ~ ., data = om_data, family = binomial())
      q_pred <- predict(om_model, newdata = covariates[val_idx, , drop = FALSE],
                        type = "response")
    }
    # Clip to [0, 1] bounds (not propensity trimming - this is an outcome model)
    q_cf[val_idx] <- pmin(pmax(q_pred, 0), 1)
  }

  list(
    pi_s0 = pi_s0_cf,
    ps = ps_cf,
    q = q_cf,
    folds = folds
  )
}


#' Compute transportable AUC with cross-fitting
#'
#' Computes the DR transportable AUC estimator using cross-fitted nuisance functions.
#'
#' @inheritParams tr_auc
#' @param K Number of folds for cross-fitting.
#' @param selection_learner Optional ml_learner for selection model.
#' @param propensity_learner Optional ml_learner for propensity model.
#' @param outcome_learner Optional ml_learner for outcome model.
#' @param parallel Logical for parallel processing.
#' @param ... Additional arguments.
#'
#' @return List with estimate and cross-fitted nuisance values.
#'
#' @keywords internal
.compute_tr_auc_crossfit <- function(predictions, outcomes, treatment, source,
                                      covariates, treatment_level, analysis,
                                      K = 5,
                                      selection_learner = NULL,
                                      propensity_learner = NULL,
                                      outcome_learner = NULL,
                                      parallel = FALSE,
                                      return_se = FALSE,
                                      ps_trim_spec = NULL,
                                      ...) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)
  n0 <- sum(source == 0)

  # Get cross-fitted nuisance functions
  cf_nuisance <- .cross_fit_transport_nuisance_auc(
    treatment = treatment,
    outcomes = outcomes,
    source = source,
    covariates = covariates,
    treatment_level = treatment_level,
    analysis = analysis,
    K = K,
    selection_learner = selection_learner,
    propensity_learner = propensity_learner,
    outcome_learner = outcome_learner,
    parallel = parallel,
    ps_trim_spec = ps_trim_spec,
    ...
  )

  pi_s0 <- cf_nuisance$pi_s0
  pi_s1 <- 1 - pi_s0
  ps <- cf_nuisance$ps
  q_hat <- cf_nuisance$q

  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)
  I_s1 <- as.numeric(source == 1)

  # Concordance indicator matrix
  ind_f <- outer(predictions, predictions, ">")

  # Compute DR AUC using cross-fitted nuisance functions
  if (analysis == "transport") {
    result <- .tr_auc_transport_cf(
      predictions = predictions,
      outcomes = outcomes,
      I_a = I_a,
      I_s0 = I_s0,
      I_s1 = I_s1,
      pi_s0 = pi_s0,
      pi_s1 = pi_s1,
      ps = ps,
      q_hat = q_hat,
      ind_f = ind_f,
      n0 = n0,
      return_se = return_se
    )
  } else {
    result <- .tr_auc_joint_cf(
      predictions = predictions,
      outcomes = outcomes,
      I_a = I_a,
      I_s0 = I_s0,
      pi_s0 = pi_s0,
      ps = ps,
      q_hat = q_hat,
      ind_f = ind_f,
      n0 = n0,
      return_se = return_se
    )
  }

  # Handle backward compatibility: result may be scalar or list
  if (is.list(result)) {
    estimate <- result$estimate
    se <- result$se
  } else {
    estimate <- result
    se <- NULL
  }

  out <- list(
    estimate = estimate,
    pi_s0 = pi_s0,
    ps = ps,
    q = q_hat,
    folds = cf_nuisance$folds
  )
  if (return_se) {
    out$se <- se
  }
  out
}


#' Transport AUC DR estimator with cross-fitted nuisance
#' @noRd
.tr_auc_transport_cf <- function(predictions, outcomes, I_a, I_s0, I_s1,
                                  pi_s0, pi_s1, ps, q_hat, ind_f, n0,
                                  return_se = FALSE) {

  n <- length(outcomes)

  # OM component for target
  w_case <- I_s0 * q_hat
  w_control <- I_s0 * (1 - q_hat)
  mat_om0 <- outer(w_case, w_control, "*")
  mat_om1 <- mat_om0 * ind_f

  # IPW component for source
  ipw_weight <- I_s1 * I_a * pi_s0 / (pi_s1 * ps)
  ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
  ipw_weight[!is.finite(ipw_weight)] <- 0

  # Augmentation: IPW applied to (Y - q_hat)
  aug_case <- ipw_weight * (outcomes - q_hat)
  aug_control <- ipw_weight * ((1 - outcomes) - (1 - q_hat))

  # Cross terms for DR
  mat_aug0_case <- outer(aug_case, w_control, "*")
  mat_aug0_control <- outer(w_case, aug_control, "*")
  mat_aug1_case <- mat_aug0_case * ind_f
  mat_aug1_control <- mat_aug0_control * ind_f

  # Combine
  num <- sum(mat_om1) + sum(mat_aug1_case) + sum(mat_aug1_control)
  denom <- sum(mat_om0) + sum(mat_aug0_case) + sum(mat_aug0_control)

  if (denom == 0 || !is.finite(denom)) {
    if (return_se) {
      return(list(estimate = NA_real_, se = NA_real_))
    }
    return(NA_real_)
  }

  estimate <- num / denom

  if (!return_se) {
    return(estimate)
  }

  # Influence function-based SE using Hajek projection (Li et al. 2022)
  # For the DR AUC estimator, we use the DeLong-like approach with DR weights
  #
  # The key insight is that for each "effective case" i and "effective control" j,

  # we compute the weighted placement values and then use their variances.
  #
  # Effective case weight: w_i = q_hat_i * I_S0_i + (Y_i - q_hat_i) * ipw_weight_i
  # Effective control weight: w_j = (1-q_hat_j) * I_S0_j + ((1-Y_j) - (1-q_hat_j)) * ipw_weight_j

  # Combined effective weights for case and control
  eff_case <- w_case + aug_case  # q_hat*I_s0 + (Y - q_hat)*ipw
  eff_control <- w_control + aug_control  # (1-q_hat)*I_s0 + ((1-Y) - (1-q_hat))*ipw

  # For variance: compute weighted placement values
  # V10_i: weighted proportion of controls with lower prediction than case i
  # V01_j: weighted proportion of cases with higher prediction than control j

  # Row-wise computation: for each i, V10_i = sum_j w_j * I(pred_i > pred_j) / sum_j w_j
  # Column-wise computation: for each j, V01_j = sum_i w_i * I(pred_i > pred_j) / sum_i w_i

  # First compute the normalizing constants
  sum_eff_case <- sum(eff_case)
  sum_eff_control <- sum(eff_control)

  if (sum_eff_case == 0 || sum_eff_control == 0) {
    return(list(estimate = estimate, se = NA_real_))
  }

  # Compute weighted placement values for each observation
  # V10: for observation i contributing as a case
  V10 <- numeric(n)
  V01 <- numeric(n)

  for (i in seq_len(n)) {
    if (eff_case[i] > 0) {
      # Weighted proportion of controls with lower prediction
      V10[i] <- sum(eff_control * ind_f[i, ]) / sum_eff_control
    }
    if (eff_control[i] > 0) {
      # Weighted proportion of cases with higher prediction
      V01[i] <- sum(eff_case * ind_f[, i]) / sum_eff_case
    }
  }

  # DeLong-like variance: weighted variance of placement values
  # Weight by effective contribution
  cases_idx <- which(eff_case > 0)
  controls_idx <- which(eff_control > 0)

  n1_eff <- length(cases_idx)
  n0_eff <- length(controls_idx)

  if (n1_eff <= 1 || n0_eff <= 1) {
    return(list(estimate = estimate, se = NA_real_))
  }

  # Weighted variance of V10 among effective cases
  V10_cases <- V10[cases_idx]
  w_cases <- eff_case[cases_idx]
  w_cases_norm <- w_cases / sum(w_cases)
  mean_V10 <- sum(w_cases_norm * V10_cases)
  S10 <- sum(w_cases_norm * (V10_cases - mean_V10)^2) * n1_eff / (n1_eff - 1)

  # Weighted variance of V01 among effective controls
  V01_controls <- V01[controls_idx]
  w_controls <- eff_control[controls_idx]
  w_controls_norm <- w_controls / sum(w_controls)
  mean_V01 <- sum(w_controls_norm * V01_controls)
  S01 <- sum(w_controls_norm * (V01_controls - mean_V01)^2) * n0_eff / (n0_eff - 1)

  # DeLong variance estimate
  se <- sqrt(S10 / n1_eff + S01 / n0_eff)

  list(estimate = estimate, se = se)
}


#' Joint AUC DR estimator with cross-fitted nuisance
#' @noRd
.tr_auc_joint_cf <- function(predictions, outcomes, I_a, I_s0,
                              pi_s0, ps, q_hat, ind_f, n0,
                              return_se = FALSE) {

  n <- length(outcomes)

  # OM component for target
  w_case <- I_s0 * q_hat
  w_control <- I_s0 * (1 - q_hat)
  mat_om0 <- outer(w_case, w_control, "*")
  mat_om1 <- mat_om0 * ind_f

  # IPW component for all with treatment
  ipw_weight <- I_a * pi_s0 / ps
  ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
  ipw_weight[!is.finite(ipw_weight)] <- 0

  # Augmentation
  aug_case <- ipw_weight * (outcomes - q_hat)
  aug_control <- ipw_weight * ((1 - outcomes) - (1 - q_hat))

  # Cross terms for DR
  mat_aug0_case <- outer(aug_case, w_control, "*")
  mat_aug0_control <- outer(w_case, aug_control, "*")
  mat_aug1_case <- mat_aug0_case * ind_f
  mat_aug1_control <- mat_aug0_control * ind_f

  # Combine
  num <- sum(mat_om1) + sum(mat_aug1_case) + sum(mat_aug1_control)
  denom <- sum(mat_om0) + sum(mat_aug0_case) + sum(mat_aug0_control)

  if (denom == 0 || !is.finite(denom)) {
    if (return_se) {
      return(list(estimate = NA_real_, se = NA_real_))
    }
    return(NA_real_)
  }

  estimate <- num / denom

  if (!return_se) {
    return(estimate)
  }

  # Influence function-based SE using Hajek projection (Li et al. 2022)
  # Same approach as transport version

  # Combined effective weights
  eff_case <- w_case + aug_case
  eff_control <- w_control + aug_control

  sum_eff_case <- sum(eff_case)
  sum_eff_control <- sum(eff_control)

  if (sum_eff_case == 0 || sum_eff_control == 0) {
    return(list(estimate = estimate, se = NA_real_))
  }

  # Compute weighted placement values
  V10 <- numeric(n)
  V01 <- numeric(n)

  for (i in seq_len(n)) {
    if (eff_case[i] > 0) {
      V10[i] <- sum(eff_control * ind_f[i, ]) / sum_eff_control
    }
    if (eff_control[i] > 0) {
      V01[i] <- sum(eff_case * ind_f[, i]) / sum_eff_case
    }
  }

  cases_idx <- which(eff_case > 0)
  controls_idx <- which(eff_control > 0)

  n1_eff <- length(cases_idx)
  n0_eff <- length(controls_idx)

  if (n1_eff <= 1 || n0_eff <= 1) {
    return(list(estimate = estimate, se = NA_real_))
  }

  # Weighted variance of V10 among effective cases
  V10_cases <- V10[cases_idx]
  w_cases <- eff_case[cases_idx]
  w_cases_norm <- w_cases / sum(w_cases)
  mean_V10 <- sum(w_cases_norm * V10_cases)
  S10 <- sum(w_cases_norm * (V10_cases - mean_V10)^2) * n1_eff / (n1_eff - 1)

  # Weighted variance of V01 among effective controls
  V01_controls <- V01[controls_idx]
  w_controls <- eff_control[controls_idx]
  w_controls_norm <- w_controls / sum(w_controls)
  mean_V01 <- sum(w_controls_norm * V01_controls)
  S01 <- sum(w_controls_norm * (V01_controls - mean_V01)^2) * n0_eff / (n0_eff - 1)

  # DeLong variance estimate
  se <- sqrt(S10 / n1_eff + S01 / n0_eff)

  list(estimate = estimate, se = se)
}


# ==============================================================================
# Internal Functions: Nuisance Model Fitting for AUC
# ==============================================================================

#' Fit nuisance models for transportability AUC estimation
#' @noRd
.fit_transport_nuisance_auc <- function(treatment, outcomes, source, covariates,
                                        treatment_level, analysis,
                                        selection_model, propensity_model,
                                        outcome_model) {

  # Determine if we're in factual mode (no treatment)
  factual_mode <- is.null(treatment)

  # Create data frame - handle factual mode (no treatment column)
  covariates_df <- as.data.frame(covariates)
  covariate_names <- names(covariates_df)

  if (factual_mode) {
    df <- cbind(
      Y = outcomes,
      S = source,
      covariates_df
    )
  } else {
    df <- cbind(
      Y = outcomes,
      A = treatment,
      S = source,
      covariates_df
    )
  }

  # Selection model: P(S=0|X) - target probability
  if (is.null(selection_model)) {
    selection_formula <- as.formula(
      paste("I(S == 0) ~", paste(covariate_names, collapse = " + "))
    )
    selection_model <- glm(selection_formula, data = df, family = binomial())
  }

  # Propensity model: P(A=1|X,S) - only for counterfactual mode
  if (is.null(propensity_model) && !factual_mode) {
    if (analysis == "transport") {
      # For transport: P(A=1|X, S=1) - propensity in source
      propensity_formula <- as.formula(
        paste("A ~", paste(covariate_names, collapse = " + "))
      )
      propensity_model <- glm(propensity_formula,
                              data = df[df$S == 1, ],
                              family = binomial())
    } else {
      # For joint: P(A=1|X) - marginal propensity
      propensity_formula <- as.formula(
        paste("A ~", paste(covariate_names, collapse = " + "))
      )
      propensity_model <- glm(propensity_formula, data = df, family = binomial())
    }
  }

  # Outcome model: E[Y|X, A=a, S] for AUC (probability of outcome)
  if (is.null(outcome_model)) {
    outcome_formula <- as.formula(
      paste("Y ~", paste(covariate_names, collapse = " + "))
    )

    if (factual_mode) {
      # Factual mode: use all source observations
      if (analysis == "transport") {
        subset_data <- df[df$S == 1, ]
      } else {
        subset_data <- df  # All data
      }
    } else {
      # Counterfactual mode: use treatment-matched observations
      if (analysis == "transport") {
        # Fit on source data with treatment level
        subset_data <- df[df$S == 1 & df$A == treatment_level, ]
      } else {
        # Fit on all data with treatment level (joint)
        subset_data <- df[df$A == treatment_level, ]
      }
    }

    if (nrow(subset_data) < 5) {
      warning("Very few observations for outcome model fitting")
    }
    outcome_model <- glm(outcome_formula, data = subset_data, family = binomial())
    # Store the full data for predictions
    attr(outcome_model, "full_data") <- df
  }

  list(
    selection = selection_model,
    propensity = propensity_model,
    outcome = outcome_model
  )
}


# ==============================================================================
# Internal Functions: AUC Computation
# ==============================================================================

#' Compute transportable AUC
#' @noRd
.compute_tr_auc <- function(predictions, outcomes, treatment, source, covariates,
                            treatment_level, analysis, estimator,
                            selection_model, propensity_model, outcome_model,
                            ps_trim_spec = NULL) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  # Dispatch to factual mode if treatment is NULL
  if (is.null(treatment)) {
    return(.compute_tr_auc_factual(
      predictions = predictions,
      outcomes = outcomes,
      source = source,
      covariates = covariates,
      analysis = analysis,
      estimator = estimator,
      selection_model = selection_model,
      outcome_model = outcome_model,
      ps_trim_spec = ps_trim_spec
    ))
  }

  if (estimator == "naive") {
    return(.compute_tr_auc_naive(predictions, outcomes, treatment, source,
                                  treatment_level, analysis))
  }

  if (analysis == "transport") {
    return(.tr_auc_transport(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      ps_trim_spec = ps_trim_spec
    ))
  } else {
    return(.tr_auc_joint(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      source = source,
      covariates = covariates,
      treatment_level = treatment_level,
      estimator = estimator,
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model,
      ps_trim_spec = ps_trim_spec
    ))
  }
}


#' Compute naive AUC for transportability setting
#' @noRd
.compute_tr_auc_naive <- function(predictions, outcomes, treatment, source,
                                   treatment_level, analysis) {
  # For transport: use source data with treatment level
  # For joint: use all data with treatment level
  # Handle factual mode (treatment is NULL)
  if (is.null(treatment)) {
    if (analysis == "transport") {
      idx <- source == 1
    } else {
      idx <- rep(TRUE, length(source))
    }
  } else {
    if (analysis == "transport") {
      idx <- source == 1 & treatment == treatment_level
    } else {
      idx <- treatment == treatment_level
    }
  }

  preds_subset <- predictions[idx]
  outcomes_subset <- outcomes[idx]

  n1 <- sum(outcomes_subset == 1)
  n0 <- sum(outcomes_subset == 0)

  if (n1 == 0 || n0 == 0) {
    warning("AUC undefined: no cases or no non-cases in relevant subset")
    return(NA_real_)
  }

  # Wilcoxon-Mann-Whitney statistic
  r <- rank(c(preds_subset[outcomes_subset == 1], preds_subset[outcomes_subset == 0]))
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


#' Compute transportable AUC for factual prediction models
#'
#' Implements estimators for transporting factual prediction models
#' (without treatment/intervention) to a target population.
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of binary outcomes (only available in source).
#' @param source Numeric vector indicating source (1) or target (0) population.
#' @param covariates Matrix or data frame of covariates.
#' @param analysis Character: "transport" or "joint".
#' @param estimator Character: "naive", "ipw", "om", or "dr".
#' @param selection_model Fitted selection model or ML learner.
#' @param outcome_model Fitted outcome model or ML learner (for om/dr).
#' @param ps_trim_spec Propensity score trimming specification.
#' @return Numeric scalar: estimated AUC in target population.
#' @noRd
.compute_tr_auc_factual <- function(predictions, outcomes, source, covariates,
                                         analysis, estimator, selection_model,
                                         outcome_model, ps_trim_spec = NULL) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)
  n0 <- sum(source == 0)

  I_s1 <- source == 1
  I_s0 <- source == 0

  if (estimator == "naive") {
    # Use source data only
    if (analysis == "transport") {
      idx <- I_s1
    } else {
      idx <- rep(TRUE, n)
    }

    preds_subset <- predictions[idx]
    outcomes_subset <- outcomes[idx]

    n1 <- sum(outcomes_subset == 1)
    n0_cases <- sum(outcomes_subset == 0)

    if (n1 == 0 || n0_cases == 0) {
      warning("AUC undefined: no cases or no non-cases in relevant subset")
      return(NA_real_)
    }

    # Wilcoxon-Mann-Whitney statistic
    r <- rank(c(preds_subset[outcomes_subset == 1], preds_subset[outcomes_subset == 0]))
    return((sum(r[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0_cases))
  }

  # Get selection probabilities
  # Note: selection model predicts P(S=0|X) (target probability)
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  if (is_ml_learner(selection_model)) {
    p_s0 <- .predict_ml_learner(selection_model, covariates)
  } else {
    p_s0 <- predict(selection_model, newdata = covariates, type = "response")
  }

  # Trim selection probabilities
  if (!is.null(ps_trim_spec)) {
    p_s0 <- .trim_propensity(p_s0, ps_trim_spec$method, ps_trim_spec$bounds)
  }

  # P(S=1|X) = 1 - P(S=0|X)
  p_s1 <- 1 - p_s0

  # Inverse-odds weights: w(X) = P(S=0|X) / P(S=1|X)
  iow <- p_s0 / p_s1

  # Source data only
  preds_s1 <- predictions[I_s1]
  outcomes_s1 <- outcomes[I_s1]
  iow_s1 <- iow[I_s1]

  # Cases and non-cases in source
  cases_idx <- which(outcomes_s1 == 1)
  noncases_idx <- which(outcomes_s1 == 0)

  n1 <- length(cases_idx)
  n0_cases <- length(noncases_idx)

  if (n1 == 0 || n0_cases == 0) {
    warning("AUC undefined: no cases or no non-cases in source")
    return(NA_real_)
  }

  if (estimator == "ipw") {
    # IPW AUC: weight concordant pairs by inverse-odds weights
    # AUC = sum_{i: Y=1, j: Y=0} w_i * w_j * I(pred_i > pred_j) / sum_{i: Y=1, j: Y=0} w_i * w_j

    numerator <- 0
    denominator <- 0

    for (i in cases_idx) {
      for (j in noncases_idx) {
        w_ij <- iow_s1[i] * iow_s1[j]
        # Concordance indicator (with 0.5 for ties)
        conc <- ifelse(preds_s1[i] > preds_s1[j], 1,
                       ifelse(preds_s1[i] == preds_s1[j], 0.5, 0))
        numerator <- numerator + w_ij * conc
        denominator <- denominator + w_ij
      }
    }

    return(numerator / denominator)

  } else if (estimator == "om") {
    # Outcome model estimator for AUC
    # This requires modeling P(Y=1|X) in source and using it for target
    # AUC = E_{target}[P(Y=1|X) > P(Y=1|X')] for random X from cases, X' from non-cases

    if (is_ml_learner(outcome_model)) {
      p_y1 <- .predict_ml_learner(outcome_model, covariates)
    } else {
      p_y1 <- predict(outcome_model, newdata = covariates, type = "response")
    }

    # Use target population predictions
    preds_target <- predictions[I_s0]
    p_y1_target <- p_y1[I_s0]

    # Estimate prevalence-weighted AUC in target
    # For OM, we use expected case/non-case probabilities
    # This is an approximation - proper OM for AUC is complex
    # Simple approach: use ranking on E[Y|X]

    n_target <- sum(I_s0)
    if (n_target < 2) {
      return(NA_real_)
    }

    # Compute all pairwise comparisons weighted by predicted probabilities
    numerator <- 0
    denominator <- 0

    for (i in 1:n_target) {
      for (j in 1:n_target) {
        if (i != j) {
          # Weight by probability of i being case and j being non-case
          w_ij <- p_y1_target[i] * (1 - p_y1_target[j])
          conc <- ifelse(preds_target[i] > preds_target[j], 1,
                         ifelse(preds_target[i] == preds_target[j], 0.5, 0))
          numerator <- numerator + w_ij * conc
          denominator <- denominator + w_ij
        }
      }
    }

    return(numerator / denominator)

  } else if (estimator == "dr") {
    # Doubly robust estimator for factual AUC transportability
    # Combines outcome model (OM) and inverse probability weighting (IPW)
    # DR = OM + IPW - Augmentation (to remove double counting)

    if (is_ml_learner(outcome_model)) {
      p_y1 <- .predict_ml_learner(outcome_model, covariates)
    } else {
      p_y1 <- predict(outcome_model, newdata = covariates, type = "response")
    }

    # Concordance indicator matrix for all observations
    ind_f <- outer(predictions, predictions, ">")

    # OM component: Expected concordance in target using outcome model
    # Weight target observations by predicted case/non-case probabilities
    p_y1_target <- p_y1[I_s0]
    w_case_om <- rep(0, n)
    w_control_om <- rep(0, n)
    w_case_om[I_s0] <- p_y1[I_s0]
    w_control_om[I_s0] <- 1 - p_y1[I_s0]

    mat_om0 <- outer(w_case_om, w_control_om, "*")
    mat_om1 <- mat_om0 * ind_f

    # IPW component: Weighted concordance from source
    # For factual mode: weight by inverse-odds only (no propensity)
    # Clip extreme weights for stability
    iow_clipped <- pmin(iow, quantile(iow[iow > 0 & is.finite(iow)], 0.99, na.rm = TRUE))
    iow_clipped[!is.finite(iow_clipped)] <- 0

    w_case_ipw <- rep(0, n)
    w_control_ipw <- rep(0, n)
    w_case_ipw[I_s1] <- iow_clipped[I_s1] * (outcomes[I_s1] == 1)
    w_control_ipw[I_s1] <- iow_clipped[I_s1] * (outcomes[I_s1] == 0)

    mat_ipw0 <- outer(w_case_ipw, w_control_ipw, "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    # Augmentation term: Remove double-counting between OM and IPW
    # This is weighted by IPW weights but uses OM predictions
    w_aug_case <- rep(0, n)
    w_aug_control <- rep(0, n)
    w_aug_case[I_s1] <- iow_clipped[I_s1] * p_y1[I_s1]
    w_aug_control[I_s1] <- iow_clipped[I_s1] * (1 - p_y1[I_s1])

    mat_aug0 <- outer(w_aug_case, w_aug_control, "*")
    mat_aug1 <- mat_aug0 * ind_f

    # DR estimator: OM + IPW - Augmentation
    num <- sum(mat_om1) + sum(mat_ipw1) - sum(mat_aug1)
    denom <- sum(mat_om0) + sum(mat_ipw0) - sum(mat_aug0)

    if (denom == 0 || !is.finite(denom)) return(NA_real_)
    return(num / denom)
  }
}


#' Transport AUC estimators (use source outcomes for target)
#' @noRd
.tr_auc_transport <- function(predictions, outcomes, treatment, source, covariates,
                               treatment_level, estimator, selection_model,
                               propensity_model, outcome_model,
                               ps_trim_spec = NULL) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)
  n0 <- sum(source == 0)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X) - needed for IPW and DR
  if (estimator %in% c("ipw", "dr")) {
    pi_s0 <- predict(selection_model, newdata = df, type = "response")
    pi_s1 <- 1 - pi_s0
  }

  # Get propensity scores P(A=a|X, S=1) - needed for IPW and DR
  if (estimator %in% c("ipw", "dr")) {
    ps_s1 <- predict(propensity_model, newdata = df, type = "response")
    if (treatment_level == 0) {
      ps_s1 <- 1 - ps_s1
    }
  }

  # Get outcome probabilities q_hat = E[Y|X, A=a, S=1] - needed for OM and DR
  if (estimator %in% c("om", "dr")) {
    q_hat <- predict(outcome_model, newdata = df, type = "response")
  }

  # Treatment indicator
  I_a <- as.numeric(treatment == treatment_level)
  I_s1 <- as.numeric(source == 1)
  I_s0 <- as.numeric(source == 0)

  # Concordance indicator matrix (predictions i > predictions j)
  ind_f <- outer(predictions, predictions, ">")

  if (estimator == "om") {
    # Outcome model estimator: weight by target selection and outcome probabilities
    # Numerator: P(Y=1|X_i, A=a, S=1) * P(Y=0|X_j, A=a, S=1) * I(pred_i > pred_j) for target
    # Use outcome model to weight pairs in target population

    # For target individuals, weight by outcome probabilities
    w_case <- I_s0 * q_hat
    w_control <- I_s0 * (1 - q_hat)

    # Construct weighted concordance
    mat_om0 <- outer(w_case, w_control, "*")
    mat_om1 <- mat_om0 * ind_f

    denom <- sum(mat_om0)
    if (denom == 0) return(NA_real_)
    return(sum(mat_om1) / denom)

  } else if (estimator == "ipw") {
    # IPW estimator: weight source observations by selection/propensity weights
    # Weight = I(S=1, A=a) * P(S=0|X) / (P(S=1|X) * P(A=a|X, S=1))

    ipw_weight <- I_s1 * I_a * pi_s0 / (pi_s1 * ps_s1)
    # Clip extreme weights
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    # Weighted concordance among source observations
    mat_ipw0 <- outer(ipw_weight * (outcomes == 1), ipw_weight * (outcomes == 0), "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    denom <- sum(mat_ipw0)
    if (denom == 0) return(NA_real_)
    return(sum(mat_ipw1) / denom)

  } else if (estimator == "dr") {
    # Doubly robust estimator combines OM and IPW

    # OM component for target
    w_case_om <- I_s0 * q_hat
    w_control_om <- I_s0 * (1 - q_hat)
    mat_om0 <- outer(w_case_om, w_control_om, "*")
    mat_om1 <- mat_om0 * ind_f

    # IPW component for source
    ipw_weight <- I_s1 * I_a * pi_s0 / (pi_s1 * ps_s1)
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    mat_ipw0 <- outer(ipw_weight * (outcomes == 1), ipw_weight * (outcomes == 0), "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    # Augmentation term (source weighted by IPW and outcome model)
    # This corrects for outcome model bias
    aug_case <- ipw_weight * q_hat
    aug_control <- ipw_weight * (1 - q_hat)
    mat_aug0 <- outer(aug_case, aug_control, "*")
    mat_aug1 <- mat_aug0 * ind_f

    # DR = OM + IPW - Augmentation
    num <- sum(mat_om1) + sum(mat_ipw1) - sum(mat_aug1)
    denom <- sum(mat_om0) + sum(mat_ipw0) - sum(mat_aug0)

    if (denom == 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


#' Joint AUC estimators (pool source and target)
#' @noRd
.tr_auc_joint <- function(predictions, outcomes, treatment, source, covariates,
                           treatment_level, estimator, selection_model,
                           propensity_model, outcome_model,
                           ps_trim_spec = NULL) {

  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)
  n0 <- sum(source == 0)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X) - needed for IPW and DR
  if (estimator %in% c("ipw", "dr")) {
    pi_s0 <- predict(selection_model, newdata = df, type = "response")
  }

  # Get propensity scores P(A=a|X) - marginal, needed for IPW and DR
  if (estimator %in% c("ipw", "dr")) {
    ps <- predict(propensity_model, newdata = df, type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
  }

  # Get outcome probabilities q_hat = E[Y|X, A=a] (joint model) - needed for OM and DR
  if (estimator %in% c("om", "dr")) {
    q_hat <- predict(outcome_model, newdata = df, type = "response")
  }

  # Treatment indicator
  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)

  # Concordance indicator matrix
  ind_f <- outer(predictions, predictions, ">")

  if (estimator == "om") {
    # Outcome model estimator for joint
    # Weight by target selection and joint outcome probabilities

    w_case <- I_s0 * q_hat
    w_control <- I_s0 * (1 - q_hat)

    mat_om0 <- outer(w_case, w_control, "*")
    mat_om1 <- mat_om0 * ind_f

    denom <- sum(mat_om0)
    if (denom == 0) return(NA_real_)
    return(sum(mat_om1) / denom)

  } else if (estimator == "ipw") {
    # IPW estimator: weight all observations with treatment level
    # Weight = I(A=a) * P(S=0|X) / P(A=a|X)

    ipw_weight <- I_a * pi_s0 / ps
    # Clip extreme weights
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    mat_ipw0 <- outer(ipw_weight * (outcomes == 1), ipw_weight * (outcomes == 0), "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    denom <- sum(mat_ipw0)
    if (denom == 0) return(NA_real_)
    return(sum(mat_ipw1) / denom)

  } else if (estimator == "dr") {
    # Doubly robust estimator for joint

    # OM component for target
    w_case_om <- I_s0 * q_hat
    w_control_om <- I_s0 * (1 - q_hat)
    mat_om0 <- outer(w_case_om, w_control_om, "*")
    mat_om1 <- mat_om0 * ind_f

    # IPW component for all with treatment level
    ipw_weight <- I_a * pi_s0 / ps
    ipw_weight <- pmin(ipw_weight, quantile(ipw_weight[ipw_weight > 0], 0.99, na.rm = TRUE))
    ipw_weight[!is.finite(ipw_weight)] <- 0

    mat_ipw0 <- outer(ipw_weight * (outcomes == 1), ipw_weight * (outcomes == 0), "*")
    mat_ipw1 <- mat_ipw0 * ind_f

    # Augmentation term
    aug_case <- ipw_weight * q_hat
    aug_control <- ipw_weight * (1 - q_hat)
    mat_aug0 <- outer(aug_case, aug_control, "*")
    mat_aug1 <- mat_aug0 * ind_f

    # DR = OM + IPW - Augmentation
    num <- sum(mat_om1) + sum(mat_ipw1) - sum(mat_aug1)
    denom <- sum(mat_om0) + sum(mat_ipw0) - sum(mat_aug0)

    if (denom == 0) return(NA_real_)
    return(num / denom)
  }

  NA_real_
}


# ==============================================================================
# Internal Functions: Bootstrap
# ==============================================================================

#' Bootstrap standard errors for transportability AUC
#' @noRd
.bootstrap_tr_auc <- function(predictions, outcomes, treatment, source, covariates,
                               treatment_level, analysis, estimator,
                               n_boot, conf_level, 
                               boot_ci_type = c("percentile", "normal", "basic"),
                               stratified = TRUE,
                               parallel, ncores, ps_trim_spec = NULL,
                               selection_model = NULL,
                               propensity_model = NULL,
                               outcome_model = NULL,
                               point_estimate = NULL, ...) {

  boot_ci_type <- match.arg(boot_ci_type)
  
  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)

  # Indices for source and target
  idx_source <- which(source == 1)
  idx_target <- which(source == 0)
  n_source <- length(idx_source)
  n_target <- length(idx_target)

  # Function to generate bootstrap indices
  get_boot_idx <- function() {
    if (stratified) {
      # Stratified bootstrap: sample separately from source and target
      boot_source <- sample(idx_source, n_source, replace = TRUE)
      boot_target <- sample(idx_target, n_target, replace = TRUE)
      c(boot_source, boot_target)
    } else {
      # Simple bootstrap
      sample(n, n, replace = TRUE)
    }
  }

  boot_fn <- function(idx) {
    tryCatch({
      # Determine if we're in factual mode (no treatment)
      factual_mode <- is.null(treatment)
      
      # For propensity and outcome models, subset appropriately
      I_s1 <- source[idx] == 1
      if (!factual_mode) {
        I_a <- treatment[idx] == treatment_level
      }

      # Prepare bootstrap data for selection model (no A or Y columns)
      # Selection model: P(S=0|X) - only needs covariates
      boot_data_sel <- data.frame(
        S0 = as.integer(source[idx] == 0),
        covariates[idx, , drop = FALSE]
      )

      if (factual_mode) {
        # Factual mode: no propensity model, outcome model uses all source
        boot_data_ps <- NULL
        if (analysis == "transport") {
          boot_data_om <- data.frame(Y = outcomes[idx], covariates[idx, , drop = FALSE])[I_s1, , drop = FALSE]
        } else {
          boot_data_om <- data.frame(Y = outcomes[idx], covariates[idx, , drop = FALSE])
        }
      } else if (analysis == "transport") {
        boot_data_ps <- data.frame(A = treatment[idx], covariates[idx, , drop = FALSE])[I_s1, , drop = FALSE]
        boot_data_om <- data.frame(Y = outcomes[idx], covariates[idx, , drop = FALSE])[I_s1 & I_a, , drop = FALSE]
      } else {
        boot_data_ps <- data.frame(A = treatment[idx], covariates[idx, , drop = FALSE])
        boot_data_om <- data.frame(Y = outcomes[idx], covariates[idx, , drop = FALSE])[I_a, , drop = FALSE]
      }

      # Refit nuisance models preserving user's formula if provided
      # Fit models based on estimator requirements
      sel_model_b <- NULL
      ps_model_b <- NULL
      out_model_b <- NULL

      # Selection model needed for all estimators except naive
      if (estimator != "naive") {
        if (!is.null(selection_model) && inherits(selection_model, "glm")) {
          # User provided a fitted model - use update() to preserve formula
          sel_model_b <- tryCatch(
            update(selection_model, data = boot_data_sel),
            error = function(e) glm(S0 ~ ., data = boot_data_sel, family = binomial())
          )
        } else {
          sel_model_b <- glm(S0 ~ ., data = boot_data_sel, family = binomial())
        }
      }

      # Propensity model needed for IPW and DR (only in counterfactual mode)
      if (estimator %in% c("ipw", "dr") && !factual_mode) {
        if (!is.null(propensity_model) && inherits(propensity_model, "glm")) {
          ps_model_b <- tryCatch(
            update(propensity_model, data = boot_data_ps),
            error = function(e) glm(A ~ ., data = boot_data_ps, family = binomial())
          )
        } else {
          ps_model_b <- glm(A ~ ., data = boot_data_ps, family = binomial())
        }
      }

      # Outcome model needed for OM and DR (for AUC, outcomes are always binary)
      if (estimator %in% c("om", "dr") && nrow(boot_data_om) > 0) {
        if (!is.null(outcome_model) && inherits(outcome_model, "glm")) {
          out_model_b <- tryCatch(
            update(outcome_model, data = boot_data_om),
            error = function(e) glm(Y ~ ., data = boot_data_om, family = binomial())
          )
        } else {
          out_model_b <- glm(Y ~ ., data = boot_data_om, family = binomial())
        }
      }

      nuisance_boot <- list(selection = sel_model_b, propensity = ps_model_b, outcome = out_model_b)

      .compute_tr_auc(
        predictions = predictions[idx],
        outcomes = outcomes[idx],
        treatment = if (factual_mode) NULL else treatment[idx],
        source = source[idx],
        covariates = covariates[idx, , drop = FALSE],
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = estimator,
        selection_model = nuisance_boot$selection,
        propensity_model = nuisance_boot$propensity,
        outcome_model = nuisance_boot$outcome,
        ps_trim_spec = ps_trim_spec
      )
    }, error = function(e) NA_real_)
  }

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    boot_estimates <- parallel::mclapply(
      1:n_boot,
      function(b) {
        set.seed(b)
        idx <- get_boot_idx()
        boot_fn(idx)
      },
      mc.cores = ncores
    )
    boot_estimates <- unlist(boot_estimates)
  } else {
    boot_estimates <- numeric(n_boot)
    for (b in 1:n_boot) {
      idx <- get_boot_idx()
      boot_estimates[b] <- boot_fn(idx)
    }
  }

  # Remove NA values
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]

  if (length(boot_estimates) < n_boot * 0.5) {
    warning("More than 50% of bootstrap samples failed")
  }

  se <- sd(boot_estimates, na.rm = TRUE)
  alpha <- 1 - conf_level
  
  # Compute confidence intervals based on specified method
  if (boot_ci_type == "percentile") {
    ci_lower <- quantile(boot_estimates, alpha / 2, na.rm = TRUE)
    ci_upper <- quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)
  } else if (boot_ci_type == "normal") {
    z <- qnorm(1 - alpha / 2)
    ci_lower <- point_estimate - z * se
    ci_upper <- point_estimate + z * se
  } else if (boot_ci_type == "basic") {
    # Basic bootstrap: 2*theta_hat - quantile(1-alpha/2), 2*theta_hat - quantile(alpha/2)
    ci_lower <- 2 * point_estimate - quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)
    ci_upper <- 2 * point_estimate - quantile(boot_estimates, alpha / 2, na.rm = TRUE)
  }

  list(
    se = se,
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper),
    boot_estimates = boot_estimates
  )
}


#' Bootstrap for cross-fitted transportable AUC
#'
#' @inheritParams .bootstrap_tr_auc
#' @param n_folds Number of folds for cross-fitting.
#' @param selection_learner Optional ml_learner for selection model.
#' @param propensity_learner Optional ml_learner for propensity model.
#' @param outcome_learner Optional ml_learner for outcome model.
#'
#' @return List with se, ci_lower, ci_upper, boot_estimates.
#' @keywords internal
.bootstrap_tr_auc_crossfit <- function(predictions, outcomes, treatment, source, covariates,
                                        treatment_level, analysis, n_folds,
                                        n_boot, conf_level,
                                        boot_ci_type = c("percentile", "normal", "basic"),
                                        stratified = TRUE,
                                        selection_learner = NULL,
                                        propensity_learner = NULL,
                                        outcome_learner = NULL,
                                        parallel, ncores, ps_trim_spec = NULL,
                                        point_estimate = NULL, ...) {

  boot_ci_type <- match.arg(boot_ci_type)
  
  # Default ps_trim_spec if not provided
  if (is.null(ps_trim_spec)) {
    ps_trim_spec <- .parse_ps_trim(NULL)
  }

  n <- length(outcomes)

  # Indices for source and target
  idx_source <- which(source == 1)
  idx_target <- which(source == 0)
  n_source <- length(idx_source)
  n_target <- length(idx_target)

  # Function to generate bootstrap indices
  get_boot_idx <- function() {
    if (stratified) {
      boot_source <- sample(idx_source, n_source, replace = TRUE)
      boot_target <- sample(idx_target, n_target, replace = TRUE)
      c(boot_source, boot_target)
    } else {
      sample(n, n, replace = TRUE)
    }
  }

  boot_fn <- function(idx) {
    tryCatch({
      .compute_tr_auc_crossfit(
        predictions = predictions[idx],
        outcomes = outcomes[idx],
        treatment = treatment[idx],
        source = source[idx],
        covariates = covariates[idx, , drop = FALSE],
        treatment_level = treatment_level,
        analysis = analysis,
        K = n_folds,
        selection_learner = selection_learner,
        propensity_learner = propensity_learner,
        outcome_learner = outcome_learner,
        parallel = FALSE,
        ps_trim_spec = ps_trim_spec,
        ...
      )$estimate
    }, error = function(e) NA_real_)
  }

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores() - 1)
    }
    boot_estimates <- parallel::mclapply(
      1:n_boot,
      function(b) {
        set.seed(b)
        idx <- get_boot_idx()
        boot_fn(idx)
      },
      mc.cores = ncores
    )
    boot_estimates <- unlist(boot_estimates)
  } else {
    boot_estimates <- numeric(n_boot)
    for (b in 1:n_boot) {
      idx <- get_boot_idx()
      boot_estimates[b] <- boot_fn(idx)
    }
  }

  # Remove NA values
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]

  if (length(boot_estimates) < n_boot * 0.5) {
    warning("More than 50% of bootstrap samples failed")
  }

  se <- sd(boot_estimates, na.rm = TRUE)
  alpha <- 1 - conf_level
  
  # Compute confidence intervals based on specified method
  if (boot_ci_type == "percentile") {
    ci_lower <- quantile(boot_estimates, alpha / 2, na.rm = TRUE)
    ci_upper <- quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)
  } else if (boot_ci_type == "normal") {
    z <- qnorm(1 - alpha / 2)
    ci_lower <- point_estimate - z * se
    ci_upper <- point_estimate + z * se
  } else if (boot_ci_type == "basic") {
    # Basic bootstrap: 2*theta_hat - quantile(1-alpha/2), 2*theta_hat - quantile(alpha/2)
    ci_lower <- 2 * point_estimate - quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)
    ci_upper <- 2 * point_estimate - quantile(boot_estimates, alpha / 2, na.rm = TRUE)
  }

  list(
    se = se,
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper),
    boot_estimates = boot_estimates
  )
}
