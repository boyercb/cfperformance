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
