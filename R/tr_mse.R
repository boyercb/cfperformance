#' Estimate Transportable Mean Squared Error
#'
#' Estimates the mean squared error (MSE) of a prediction model in a target
#' population using data transported from a source population (typically an RCT).
#'
#' @param predictions Numeric vector of model predictions.
#' @param outcomes Numeric vector of observed outcomes.
#' @param treatment Numeric vector of treatment indicators (0/1).
#' @param source Numeric vector of population indicators (1=source/RCT, 0=target).
#' @param covariates A matrix or data frame of baseline covariates.
#' @param treatment_level The treatment level of interest (default: 0).
#' @param analysis Character string specifying the type of analysis:
#'   - `"transport"`: Use source outcomes for target estimation (default)
#'   - `"joint"`: Pool source and target data
#' @param estimator Character string specifying the estimator:
#'   - `"naive"`: Naive estimator (biased)
#'   - `"om"`: Outcome model estimator
#'   - `"ipw"`: Inverse probability weighting estimator
#'   - `"dr"`: Doubly robust estimator (default)
#' @param selection_model Optional fitted selection model for P(S=0|X). If NULL,
#'   a logistic regression model is fit using the covariates.
#' @param propensity_model Optional fitted propensity score model for P(A=1|X,S=1).
#'   If NULL, a logistic regression model is fit using source data.
#' @param outcome_model Optional fitted outcome model for E\[L(Y,g)|X,A,S\].
#'   If NULL, a regression model is fit using the relevant data.
#' @param se_method Method for standard error estimation:
#'   - `"bootstrap"`: Bootstrap standard errors (default)
#'   - `"influence"`: Influence function-based standard errors
#'   - `"none"`: No standard error estimation
#' @param n_boot Number of bootstrap replications (default: 500).
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param stratified_boot Logical indicating whether to use stratified bootstrap
#'   that preserves the source/target ratio (default: TRUE). Recommended for
#'   transportability analysis.
#' @param cross_fit Logical indicating whether to use cross-fitting for
#'   nuisance model estimation (default: FALSE).
#' @param n_folds Number of folds for cross-fitting (default: 5).
#' @param parallel Logical indicating whether to use parallel processing
#'   for bootstrap (default: FALSE).
#' @param ncores Number of cores for parallel processing (default: NULL,
#'   which uses all available cores minus one).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class `c("tr_mse", "tr_performance")` containing:
#'   \item{estimate}{Point estimate of transportable MSE}
#'   \item{se}{Standard error (if computed)}
#'   \item{ci_lower}{Lower confidence interval bound}
#'   \item{ci_upper}{Upper confidence interval bound}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{naive_estimate}{Naive MSE for comparison}
#'   \item{n_target}{Number of target observations}
#'   \item{n_source}{Number of source observations}
#'   \item{treatment_level}{Treatment level}
#'
#' @details
#' This function implements estimators for transporting prediction model
#' performance from a source population (typically an RCT) to a target
#' population.
#'
#' **Transportability Analysis**: Uses outcome data from the source/RCT
#' population to estimate performance in the target population. Requires:
#' - Selection model: P(S=0|X)
#' - Propensity model in source: P(A=1|X, S=1)
#' - Outcome model trained on source data
#'
#' **Joint Analysis**: Pools source and target data to estimate performance
#' in the target population. More efficient when both populations have
#' outcome data.
#'
#' For observational analysis (single population), use [cf_mse()] instead.
#'
#' @references
#' Voter, S. R., et al. (2025). "Transportability of machine learning-based
#' counterfactual prediction models with application to CASS."
#' *Diagnostic and Prognostic Research*, 9(4).
#' \doi{10.1186/s41512-025-00201-y}
#'
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @seealso [cf_mse()], [tr_auc()]
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
#' # Outcome
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' # Predictions from some model
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Estimate transportable MSE
#' result <- tr_mse(
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
tr_mse <- function(predictions,
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

  n <- length(outcomes)
  n_target <- sum(source == 0)
  n_source <- sum(source == 1)

  # Initialize SE variables
  se <- NULL
  ci_lower <- NULL
  ci_upper <- NULL

  # Fit nuisance models if not provided
  if (estimator != "naive") {
    nuisance <- .fit_transport_nuisance(
      treatment = treatment,
      outcomes = outcomes,
      source = source,
      covariates = covariates,
      predictions = predictions,
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
  estimate <- .compute_tr_mse(
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

  # Compute naive estimate for comparison (among target with treatment_level)
  naive_estimate <- .compute_tr_mse(
    predictions = predictions,
    outcomes = outcomes,
    treatment = treatment,
    source = source,
    covariates = covariates,
    treatment_level = treatment_level,
    analysis = analysis,
    estimator = "naive",
    selection_model = NULL,
    propensity_model = NULL,
    outcome_model = NULL
  )

  # Compute standard errors
  if (se_method == "bootstrap") {
    boot_result <- .bootstrap_tr_mse(
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
    se <- .influence_se_tr_mse(
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
    z <- qnorm(1 - (1 - conf_level) / 2)
    ci_lower <- estimate - z * se
    ci_upper <- estimate + z * se
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
    metric = "mse",
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

  class(result) <- c("tr_mse", "tr_performance")
  return(result)
}


# =============================================================================
# Internal functions for computing transportable MSE
# =============================================================================

#' Compute transportable MSE estimate
#' @noRd
.compute_tr_mse <- function(predictions, outcomes, treatment, source, covariates,
                            treatment_level, analysis, estimator,
                            selection_model, propensity_model, outcome_model) {

  n <- length(outcomes)
  n0 <- sum(source == 0)  # Target population size
  loss <- (outcomes - predictions)^2

  # Indicators
  I_s0 <- source == 0  # Target population
  I_s1 <- source == 1  # Source population
  I_a <- treatment == treatment_level

  if (estimator == "naive") {
    # Naive: mean loss among target with treatment_level
    if (analysis == "transport") {
      # For transport, we don't have outcomes in target, use source
      return(mean(loss[I_s1 & I_a]))
    } else {
      # For joint, use all with treatment_level
      return(mean(loss[I_a]))
    }
  }

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  # Get selection scores P(S=0|X)
  if (!is.null(selection_model)) {
    p_s0 <- predict(selection_model, newdata = covariates, type = "response")
    p_s0 <- pmax(pmin(p_s0, 0.99), 0.01)
    p_s1 <- 1 - p_s0
  }

  # Get propensity scores
  if (!is.null(propensity_model)) {
    ps <- predict(propensity_model, newdata = covariates, type = "response")
    if (treatment_level == 0) {
      ps <- 1 - ps
    }
    ps <- pmax(pmin(ps, 0.99), 0.01)
  }

  # Get outcome model predictions (conditional expected loss)
  if (!is.null(outcome_model)) {
    h <- predict(outcome_model, newdata = covariates, type = "response")
  }

  # Compute estimator based on analysis type
  if (analysis == "transport") {
    return(.tr_mse_transport(
      loss = loss, I_s0 = I_s0, I_s1 = I_s1, I_a = I_a,
      p_s0 = p_s0, p_s1 = p_s1, ps = ps, h = h,
      estimator = estimator, n0 = n0
    ))
  } else {  # joint
    return(.tr_mse_joint(
      loss = loss, I_s0 = I_s0, I_a = I_a,
      p_s0 = p_s0, ps = ps, h = h,
      estimator = estimator, n0 = n0
    ))
  }
}


#' Transportability estimators (using source data for target estimation)
#' @noRd
.tr_mse_transport <- function(loss, I_s0, I_s1, I_a, p_s0, p_s1, ps, h,
                               estimator, n0) {

  if (estimator == "om") {
    # Outcome model estimator: E[h(X) | S=0]
    return(mean(h[I_s0]))

  } else if (estimator == "ipw") {
    # IPW estimator
    # psi = (1/n0) * sum_{S=1, A=a} [P(S=0|X) / (P(A=a|X,S=1) * P(S=1|X))] * L
    weights <- (p_s0 / (ps * p_s1))[I_s1 & I_a]
    return(sum(weights * loss[I_s1 & I_a]) / n0)

  } else if (estimator == "dr") {
    # Doubly robust estimator
    # DR = (1/n0) * [sum_{S=0} h(X) + sum_{S=1,A=a} w * (L - h)]
    # where w = P(S=0|X) / (P(S=1|X) * P(A=a|X,S=1))

    # First term: outcome model in target
    term1 <- sum(h[I_s0])

    # Second term: augmentation from source
    idx_s1_a <- I_s1 & I_a
    weights <- (p_s0 / (p_s1 * ps))[idx_s1_a]
    term2 <- sum(weights * (loss[idx_s1_a] - h[idx_s1_a]))

    return((term1 + term2) / n0)
  }
}


#' Joint estimators (pooling source and target data)
#' @noRd
.tr_mse_joint <- function(loss, I_s0, I_a, p_s0, ps, h, estimator, n0) {

  if (estimator == "om") {
    # Outcome model estimator: E[h(X) | S=0]
    return(mean(h[I_s0]))

  } else if (estimator == "ipw") {
    # IPW estimator
    # psi = (1/n0) * sum_{A=a} [P(S=0|X) / P(A=a|X)] * L
    weights <- (p_s0 / ps)[I_a]
    return(sum(weights * loss[I_a]) / n0)

  } else if (estimator == "dr") {
    # Doubly robust estimator
    # DR = (1/n0) * [sum_{S=0} h(X) + sum_{A=a} (P(S=0|X)/P(A=a|X)) * (L - h)]

    # First term: outcome model in target
    term1 <- sum(h[I_s0])

    # Second term: augmentation from all with treatment
    weights <- (p_s0 / ps)[I_a]
    term2 <- sum(weights * (loss[I_a] - h[I_a]))

    return((term1 + term2) / n0)
  }
}


# =============================================================================
# Nuisance model fitting for transportability
# =============================================================================

#' Fit nuisance models for transportability analysis
#' @noRd
.fit_transport_nuisance <- function(treatment, outcomes, source, covariates,
                                     predictions, treatment_level, analysis,
                                     selection_model, propensity_model,
                                     outcome_model) {

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  n <- length(outcomes)

  # Indicators
  I_s0 <- source == 0
  I_s1 <- source == 1
  I_a <- treatment == treatment_level

  # Compute loss for outcome modeling
  loss <- (outcomes - predictions)^2

  # --- Selection model: P(S=0|X) ---
  if (is.null(selection_model)) {
    sel_data <- cbind(S0 = as.numeric(I_s0), covariates)
    selection_model <- glm(S0 ~ ., data = sel_data, family = binomial())
  }

  # --- Propensity model ---
  if (is.null(propensity_model)) {
    if (analysis == "transport") {
      # P(A=1|X, S=1) - fit on source population only
      ps_data <- cbind(A = treatment, covariates)[I_s1, ]
      propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
    } else {
      # P(A=1|X) - fit on all data for joint analysis
      ps_data <- cbind(A = treatment, covariates)
      propensity_model <- glm(A ~ ., data = ps_data, family = binomial())
    }
  }

  # --- Outcome model: E[L(Y, g(X)) | X, A=a, S] ---
  if (is.null(outcome_model)) {
    if (analysis == "transport") {
      # Model E[L | X, A=a, S=1] using source data
      om_data <- cbind(L = loss, covariates)[I_s1 & I_a, ]
      outcome_model <- glm(L ~ ., data = om_data, family = gaussian())
    } else {
      # Model E[L | X, A=a] using all data with treatment_level
      om_data <- cbind(L = loss, covariates)[I_a, ]
      outcome_model <- glm(L ~ ., data = om_data, family = gaussian())
    }
  }

  list(
    selection = selection_model,
    propensity = propensity_model,
    outcome = outcome_model
  )
}


# =============================================================================
# Input validation for transportability
# =============================================================================

#' Validate inputs for transportability functions
#' @noRd
.validate_transport_inputs <- function(predictions, outcomes, treatment,
                                        source, covariates) {
  n <- length(outcomes)

  if (length(predictions) != n) {
    stop("predictions and outcomes must have the same length")
  }
  if (length(treatment) != n) {
    stop("treatment and outcomes must have the same length")
  }
  if (length(source) != n) {
    stop("source and outcomes must have the same length")
  }
  if (nrow(covariates) != n) {
    stop("covariates must have the same number of rows as outcomes")
  }
  if (!all(treatment %in% c(0, 1))) {
    stop("treatment must be binary (0/1)")
  }
  if (!all(source %in% c(0, 1))) {
    stop("source must be binary (0=target, 1=source)")
  }
  if (sum(source == 0) == 0) {
    stop("No target population observations (source=0)")
  }
  if (sum(source == 1) == 0) {
    stop("No source population observations (source=1)")
  }
  if (any(is.na(predictions)) || any(is.na(outcomes)) ||
      any(is.na(treatment)) || any(is.na(source))) {
    stop("Missing values not allowed in predictions, outcomes, treatment, or source")
  }
}


# =============================================================================
# Bootstrap for transportability MSE
# =============================================================================

#' Bootstrap standard errors for transportable MSE
#'
#' @param stratified Logical indicating whether to use stratified bootstrap
#'   that samples separately within source and target populations, preserving
#'   the original ratio.
#' @noRd
.bootstrap_tr_mse <- function(predictions, outcomes, treatment, source,
                               covariates, treatment_level, analysis,
                               estimator, n_boot, conf_level,
                               stratified = TRUE, parallel, ncores, ...) {

  n <- length(outcomes)

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

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

  # Bootstrap function - fits nuisance models and computes estimate
  boot_fn <- function(idx) {
    # Subset data
    pred_b <- predictions[idx]
    out_b <- outcomes[idx]
    trt_b <- treatment[idx]
    src_b <- source[idx]
    cov_b <- covariates[idx, , drop = FALSE]

    # Check we have both source and target in bootstrap sample
    if (sum(src_b == 0) < 5 || sum(src_b == 1) < 5) {
      return(NA_real_)
    }

    # Fit nuisance models on bootstrap sample
    nuisance_b <- .fit_transport_nuisance(
      treatment = trt_b,
      outcomes = out_b,
      source = src_b,
      covariates = cov_b,
      predictions = pred_b,
      treatment_level = treatment_level,
      analysis = analysis,
      selection_model = NULL,
      propensity_model = NULL,
      outcome_model = NULL
    )

    # Compute estimate
    .compute_tr_mse(
      predictions = pred_b,
      outcomes = out_b,
      treatment = trt_b,
      source = src_b,
      covariates = cov_b,
      treatment_level = treatment_level,
      analysis = analysis,
      estimator = estimator,
      selection_model = nuisance_b$selection,
      propensity_model = nuisance_b$propensity,
      outcome_model = nuisance_b$outcome
    )
  }

  # Run bootstrap
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    ncores <- ncores %||% (parallel::detectCores() - 1)
    boot_estimates <- parallel::mclapply(
      1:n_boot,
      function(b) {
        idx <- get_boot_idx()
        tryCatch(boot_fn(idx), error = function(e) NA_real_)
      },
      mc.cores = ncores
    )
    boot_estimates <- unlist(boot_estimates)
  } else {
    boot_estimates <- vapply(1:n_boot, function(b) {
      idx <- get_boot_idx()
      tryCatch(boot_fn(idx), error = function(e) NA_real_)
    }, numeric(1))
  }

  # Remove failed bootstrap samples
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]

  # Compute SE and CI
  se <- sd(boot_estimates)
  alpha <- 1 - conf_level
  ci_lower <- quantile(boot_estimates, alpha / 2)
  ci_upper <- quantile(boot_estimates, 1 - alpha / 2)

  list(
    se = se,
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    boot_estimates = boot_estimates
  )
}


# =============================================================================
# Influence function-based SE
# =============================================================================

#' Influence function-based SE for transportable MSE
#'
#' Computes standard errors using the efficient influence function for
#' the doubly robust transportability MSE estimator.
#'
#' @details
#' For the transport DR estimator, the efficient influence function is:
#' \deqn{\phi_i = I(S_i=0) h_a(X_i) + I(S_i=1, A_i=a) \frac{P(S=0|X_i)}{P(S=1|X_i) P(A=a|X_i,S=1)} (L_i - h_a(X_i)) - \psi}
#'
#' For the joint DR estimator:
#' \deqn{\phi_i = I(S_i=0) h_a(X_i) + I(A_i=a) \frac{P(S=0|X_i)}{P(A=a|X_i)} (L_i - h_a(X_i)) - \psi}
#'
#' The variance is estimated as Var(phi) / n.
#'
#' @noRd
.influence_se_tr_mse <- function(predictions, outcomes, treatment, source,
                                  covariates, treatment_level, analysis,
                                  estimator, selection_model, propensity_model,
                                  outcome_model) {

  n <- length(outcomes)
  n0 <- sum(source == 0)

  # Compute loss
  loss <- (outcomes - predictions)^2

  # Treatment and source indicators
  I_a <- as.numeric(treatment == treatment_level)
  I_s0 <- as.numeric(source == 0)
  I_s1 <- as.numeric(source == 1)

  # Convert covariates to data frame if needed
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }

  if (estimator == "naive") {
    # For naive estimator, influence is centered loss in relevant subset
    if (analysis == "transport") {
      # Use source with treatment level
      idx <- I_s1 == 1 & I_a == 1
      psi_hat <- mean(loss[idx])
      # Influence function for mean of subset
      phi <- rep(0, n)
      phi[idx] <- loss[idx] - psi_hat
      se <- sqrt(var(phi[idx]) / sum(idx))
    } else {
      # Joint: use all with treatment level
      idx <- I_a == 1
      psi_hat <- mean(loss[idx])
      phi <- rep(0, n)
      phi[idx] <- loss[idx] - psi_hat
      se <- sqrt(var(phi[idx]) / sum(idx))
    }
    return(se)
  }

  # Get selection probabilities P(S=0|X)
  pi_s0 <- predict(selection_model, newdata = covariates, type = "response")
  pi_s1 <- 1 - pi_s0
  # Truncate for stability
  pi_s0 <- pmax(pmin(pi_s0, 0.99), 0.01)
  pi_s1 <- pmax(pmin(pi_s1, 0.99), 0.01)

  # Get propensity scores
  ps <- predict(propensity_model, newdata = covariates, type = "response")
  if (treatment_level == 0) {
    ps <- 1 - ps
  }
  ps <- pmax(pmin(ps, 0.99), 0.01)

  # Get outcome model predictions h_a(X) = E[L|X, A=a]
  # For MSE: h = E[Y|X,A=a] - 2*pred*E[Y|X,A=a] + pred^2
  q_hat <- predict(outcome_model, newdata = covariates, type = "response")
  # For loss model: h = E[L|X, A=a, S] directly
  h_a <- q_hat

  if (analysis == "transport") {
    # Transport estimator influence function
    if (estimator == "om") {
      # OM: phi = I(S=0) * h_a(X) - psi, normalized
      psi_hat <- mean(I_s0 * h_a) / mean(I_s0)
      phi <- (I_s0 * h_a - psi_hat * I_s0)
      # Variance of mean of target
      se <- sqrt(var(phi[I_s0 == 1]) / n0)

    } else if (estimator == "ipw") {
      # IPW: phi = I(S=1, A=a) * [P(S=0|X)/(P(S=1|X)*P(A|X,S=1))] * L - psi
      ipw_weight <- I_s1 * I_a * pi_s0 / (pi_s1 * ps)
      psi_hat <- sum(ipw_weight * loss) / n0
      phi <- (ipw_weight * loss - psi_hat * I_s0) / mean(I_s0)
      se <- sqrt(var(phi) / n)

    } else if (estimator == "dr") {
      # DR: phi = I(S=0)*h_a(X) + I(S=1,A=a)*[P(S=0|X)/(P(S=1|X)*P(A|X,S=1))]*(L-h_a) - psi
      ipw_weight <- I_s1 * I_a * pi_s0 / (pi_s1 * ps)
      augmentation <- ipw_weight * (loss - h_a)
      psi_hat <- (sum(I_s0 * h_a) + sum(augmentation)) / n0

      # Influence function (scaled by n0/n for target population inference)
      phi <- (I_s0 * h_a + augmentation - psi_hat * I_s0)
      # Divide by n0/n to get proper scaling
      phi <- phi / mean(I_s0)

      se <- sqrt(var(phi) / n)
    }

  } else {
    # Joint estimator influence function
    if (estimator == "om") {
      # OM: phi = I(S=0) * h_a(X) - psi
      psi_hat <- mean(I_s0 * h_a) / mean(I_s0)
      phi <- (I_s0 * h_a - psi_hat * I_s0)
      se <- sqrt(var(phi[I_s0 == 1]) / n0)

    } else if (estimator == "ipw") {
      # IPW: phi = I(A=a) * [P(S=0|X)/P(A|X)] * L - psi
      ipw_weight <- I_a * pi_s0 / ps
      psi_hat <- sum(ipw_weight * loss) / n0
      phi <- (ipw_weight * loss - psi_hat * I_s0) / mean(I_s0)
      se <- sqrt(var(phi) / n)

    } else if (estimator == "dr") {
      # DR: phi = I(S=0)*h_a(X) + I(A=a)*[P(S=0|X)/P(A|X)]*(L-h_a) - psi
      ipw_weight <- I_a * pi_s0 / ps
      augmentation <- ipw_weight * (loss - h_a)
      psi_hat <- (sum(I_s0 * h_a) + sum(augmentation)) / n0

      # Influence function
      phi <- (I_s0 * h_a + augmentation - psi_hat * I_s0) / mean(I_s0)

      se <- sqrt(var(phi) / n)
    }
  }

  return(se)
}
