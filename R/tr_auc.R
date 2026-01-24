#' Estimate (Counterfactual) Area Under the ROC Curve in the Target Population
#'
#' Estimates the area under the receiver operating characteristic curve (AUC)
#' of a prediction model in a target population using data transported from
#' a source population (typically an RCT).
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
#'   \item{treatment_level}{Treatment level}
#'
#' @details
#' This function implements estimators for transporting prediction model
#' AUC from a source population (typically an RCT) to a target population.
#' The AUC is defined as the probability that a randomly selected case
#' has a higher predicted risk than a randomly selected non-case.
#'
#' **Transportability Analysis**: Uses outcome data from the source/RCT
#' population to estimate AUC in the target population. Requires:
#' - Selection model: P(S=0|X)
#' - Propensity model in source: P(A=1|X, S=1)
#' - Outcome model trained on source data: E\[Y|X, A, S=1\]
#'
#' **Joint Analysis**: Pools source and target data to estimate AUC
#' in the target population. More efficient when both populations have
#' outcome data.
#'
#' For observational analysis (single population), use [cf_auc()] instead.
#'
#' The estimators use U-statistic formulations that weight concordant pairs
#' by their estimated probabilities of occurring in the target population
#' under the counterfactual treatment.
#'
#' @references
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
                   treatment,
                   source,
                   covariates,
                   treatment_level = 0,
                   analysis = c("transport", "joint"),
                   estimator = c("dr", "om", "ipw", "naive"),
                   selection_model = NULL,
                   propensity_model = NULL,
                   outcome_model = NULL,
                   se_method = c("bootstrap", "influence", "none"),
                   n_boot = 500,
                   conf_level = 0.95,
                   stratified_boot = TRUE,
                   cross_fit = FALSE,
                   n_folds = 5,
                   parallel = FALSE,
                   ncores = NULL,
                   ...) {

  # Input validation
  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)
  se_method <- match.arg(se_method)

  .validate_transport_inputs(predictions, outcomes, treatment, source, covariates)

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
        stratified = stratified_boot,
        selection_learner = selection_model,
        propensity_learner = propensity_model,
        outcome_learner = outcome_model,
        parallel = parallel,
        ncores = ncores,
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
      outcome_model = nuisance$outcome
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
        stratified = stratified_boot,
        parallel = parallel,
        ncores = ncores,
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
                                               ...) {

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
    pi_s0_cf[val_idx] <- pmax(pmin(pi_s0_pred, 0.99), 0.01)

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
    ps_cf[val_idx] <- pmax(pmin(ps_pred, 0.975), 0.025)

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
    q_cf[val_idx] <- pmax(pmin(q_pred, 0.99), 0.01)
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
                                      ...) {

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

  # Propensity model: P(A=1|X,S) depends on analysis type
  if (is.null(propensity_model)) {
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

    if (analysis == "transport") {
      # Fit on source data with treatment level
      subset_data <- df[df$S == 1 & df$A == treatment_level, ]
      if (nrow(subset_data) < 5) {
        warning("Very few observations for outcome model fitting")
      }
      outcome_model <- glm(outcome_formula, data = subset_data, family = binomial())
    } else {
      # Fit on all data with treatment level (joint)
      subset_data <- df[df$A == treatment_level, ]
      outcome_model <- glm(outcome_formula, data = subset_data, family = binomial())
    }
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
                            selection_model, propensity_model, outcome_model) {

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
      outcome_model = outcome_model
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
      outcome_model = outcome_model
    ))
  }
}


#' Compute naive AUC for transportability setting
#' @noRd
.compute_tr_auc_naive <- function(predictions, outcomes, treatment, source,
                                   treatment_level, analysis) {
  # For transport: use source data with treatment level
  # For joint: use all data with treatment level
  if (analysis == "transport") {
    idx <- source == 1 & treatment == treatment_level
  } else {
    idx <- treatment == treatment_level
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


#' Transport AUC estimators (use source outcomes for target)
#' @noRd
.tr_auc_transport <- function(predictions, outcomes, treatment, source, covariates,
                               treatment_level, estimator, selection_model,
                               propensity_model, outcome_model) {

  n <- length(outcomes)
  n0 <- sum(source == 0)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X)
  pi_s0 <- predict(selection_model, newdata = df, type = "response")
  pi_s1 <- 1 - pi_s0

  # Get propensity scores P(A=a|X, S=1)
  ps_s1 <- predict(propensity_model, newdata = df, type = "response")
  if (treatment_level == 0) {
    ps_s1 <- 1 - ps_s1
  }

  # Get outcome probabilities q_hat = E[Y|X, A=a, S=1]
  q_hat <- predict(outcome_model, newdata = df, type = "response")

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
                           propensity_model, outcome_model) {

  n <- length(outcomes)
  n0 <- sum(source == 0)
  df <- cbind(Y = outcomes, A = treatment, S = source, as.data.frame(covariates))

  # Get selection probabilities P(S=0|X)
  pi_s0 <- predict(selection_model, newdata = df, type = "response")

  # Get propensity scores P(A=a|X) - marginal
  ps <- predict(propensity_model, newdata = df, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }

  # Get outcome probabilities q_hat = E[Y|X, A=a] (joint model)
  q_hat <- predict(outcome_model, newdata = df, type = "response")

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
                               n_boot, conf_level, stratified = TRUE,
                               parallel, ncores, ...) {

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
      # Re-fit nuisance models on bootstrap sample
      if (estimator != "naive") {
        nuisance_boot <- .fit_transport_nuisance_auc(
          treatment = treatment[idx],
          outcomes = outcomes[idx],
          source = source[idx],
          covariates = covariates[idx, , drop = FALSE],
          treatment_level = treatment_level,
          analysis = analysis,
          selection_model = NULL,
          propensity_model = NULL,
          outcome_model = NULL
        )
      } else {
        nuisance_boot <- list(selection = NULL, propensity = NULL, outcome = NULL)
      }

      .compute_tr_auc(
        predictions = predictions[idx],
        outcomes = outcomes[idx],
        treatment = treatment[idx],
        source = source[idx],
        covariates = covariates[idx, , drop = FALSE],
        treatment_level = treatment_level,
        analysis = analysis,
        estimator = estimator,
        selection_model = nuisance_boot$selection,
        propensity_model = nuisance_boot$propensity,
        outcome_model = nuisance_boot$outcome
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
  ci_lower <- quantile(boot_estimates, alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)

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
                                        n_boot, conf_level, stratified = TRUE,
                                        selection_learner = NULL,
                                        propensity_learner = NULL,
                                        outcome_learner = NULL,
                                        parallel, ncores, ...) {

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
  ci_lower <- quantile(boot_estimates, alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)

  list(
    se = se,
    ci_lower = unname(ci_lower),
    ci_upper = unname(ci_upper),
    boot_estimates = boot_estimates
  )
}
