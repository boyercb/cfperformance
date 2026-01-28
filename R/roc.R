#' Compute Transportable ROC Curve
#'
#' Computes a receiver operating characteristic (ROC) curve in a target
#' population using data transported from a source population. Returns
#' sensitivity (TPR) and false positive rate (FPR) at multiple thresholds.
#' Supports both **counterfactual** and **factual** prediction model
#' transportability.
#'
#' @inheritParams tr_sensitivity
#' @param n_thresholds Integer specifying the number of thresholds to evaluate.
#'   Thresholds are evenly spaced between 0 and 1. Default is 201.
#' @param thresholds Optional numeric vector of specific thresholds to use.
#'   If provided, overrides `n_thresholds`.
#' @param include_naive Logical indicating whether to also compute the naive
#'   ROC curve for comparison. Default is TRUE.
#'
#' @return An object of class `c("tr_roc", "roc_curve")` containing:
#'   \item{thresholds}{Thresholds used}
#'   \item{sensitivity}{Sensitivity (TPR) at each threshold}
#'   \item{fpr}{False positive rate at each threshold}
#'   \item{specificity}{Specificity at each threshold}
#'   \item{naive_sensitivity}{Naive sensitivity (if include_naive=TRUE)}
#'   \item{naive_fpr}{Naive FPR (if include_naive=TRUE)}
#'   \item{auc}{Area under the ROC curve (computed via trapezoidal rule)}
#'   \item{naive_auc}{Naive AUC (if include_naive=TRUE)}
#'   \item{estimator}{Estimator used}
#'   \item{analysis}{Analysis type}
#'   \item{n_source}{Number of source observations}
#'   \item{n_target}{Number of target observations}
#'   \item{treatment_level}{Treatment level (NULL for factual mode)}
#'
#' @details
#' The ROC curve plots sensitivity (true positive rate) against the false
#' positive rate (1 - specificity) at various classification thresholds.
#'
#' ## Counterfactual Mode (treatment provided)
#' When `treatment` is specified, computes the ROC curve for counterfactual
#' outcomes under a hypothetical intervention.
#'
#' ## Factual Mode (treatment = NULL)
#' When `treatment` is `NULL`, computes the ROC curve for observed outcomes
#' in the target population using inverse-odds weighting based on the
#' selection model only.
#'
#' This function computes transportable sensitivity and FPR at multiple
#' thresholds using the estimators from [tr_sensitivity()] and [tr_fpr()].
#' The area under the curve (AUC) is computed using the trapezoidal rule on
#' the discrete threshold grid. For exact AUC estimation, use [tr_auc()]
#' which employs the Wilcoxon-Mann-Whitney statistic.
#'
#' For efficient computation, all thresholds are evaluated in a single pass
#' through the data, with nuisance models fitted only once.
#'
#' @references
#' Steingrimsson, J. A., et al. (2023). "Transporting a Prediction Model for
#' Use in a New Target Population." *American Journal of Epidemiology*,
#' 192(2), 296-304. \doi{10.1093/aje/kwac128}
#'
#' Steingrimsson, J. A., Wen, L., Voter, S., & Dahabreh, I. J. (2024).
#' "Interpretable meta-analysis of model or marker performance."
#' *arXiv preprint arXiv:2409.13458*.
#'
#' @seealso [tr_sensitivity()], [tr_specificity()], [tr_fpr()], [plot.tr_roc()], [tr_auc()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
#' a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
#' y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' # Compute transportable ROC curve
#' roc <- tr_roc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   source = s,
#'   covariates = data.frame(x = x),
#'   n_thresholds = 51
#' )
#' print(roc)
#'
#' # Plot the ROC curve
#' plot(roc)
tr_roc <- function(predictions,
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
                   n_thresholds = 201,
                   thresholds = NULL,
                   include_naive = TRUE,
                   ...) {

  analysis <- match.arg(analysis)
  estimator <- match.arg(estimator)

  # Validate inputs
  .validate_transport_inputs(predictions, outcomes, treatment, source, covariates, treatment_level)
  
  # Determine mode: factual (no treatment) or counterfactual
  factual_mode <- is.null(treatment)
  
  # Default treatment_level to 1 for backward compatibility if treatment is provided
  if (!factual_mode && is.null(treatment_level)) {
    treatment_level <- 1
  }

  if (!all(outcomes %in% c(0, 1))) {
    stop("ROC curve requires binary outcomes (0/1)")
  }

  # Set up thresholds

  if (is.null(thresholds)) {
    thresholds <- seq(0, 1, length.out = n_thresholds)
  }
  thresholds <- sort(unique(thresholds))

  n <- length(outcomes)
  n_source <- sum(source == 1)
  n_target <- sum(source == 0)

  # Fit nuisance models once
  if (estimator != "naive") {
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

  # Compute sensitivity at all thresholds
  sensitivity <- sapply(thresholds, function(c) {
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
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
  })

  # Compute specificity at all thresholds
  specificity <- sapply(thresholds, function(c) {
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
      selection_model = selection_model,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
  })

  # FPR = 1 - specificity
  fpr <- 1 - specificity

  # Compute AUC using trapezoidal rule
  # Sort by FPR for proper integration
  ord <- order(fpr)
  auc <- .compute_auc_trapezoid(fpr[ord], sensitivity[ord])

  # Naive ROC curve
  naive_sensitivity <- NULL
  naive_fpr <- NULL
  naive_auc <- NULL

  if (include_naive) {
    naive_sensitivity <- sapply(thresholds, function(c) {
      .compute_tr_sens_spec_naive(predictions, outcomes, treatment, source,
                                   c, treatment_level, analysis, "sensitivity")
    })

    naive_specificity <- sapply(thresholds, function(c) {
      .compute_tr_sens_spec_naive(predictions, outcomes, treatment, source,
                                   c, treatment_level, analysis, "specificity")
    })

    naive_fpr <- 1 - naive_specificity

    ord_naive <- order(naive_fpr)
    naive_auc <- .compute_auc_trapezoid(naive_fpr[ord_naive], naive_sensitivity[ord_naive])
  }

  result <- list(
    thresholds = thresholds,
    sensitivity = sensitivity,
    fpr = fpr,
    specificity = specificity,
    naive_sensitivity = naive_sensitivity,
    naive_fpr = naive_fpr,
    auc = auc,
    naive_auc = naive_auc,
    estimator = estimator,
    analysis = analysis,
    treatment_level = treatment_level,
    n_source = n_source,
    n_target = n_target,
    include_naive = include_naive
  )

  class(result) <- c("tr_roc", "roc_curve")
  return(result)
}


#' Compute Counterfactual ROC Curve
#'
#' Computes a receiver operating characteristic (ROC) curve under a
#' hypothetical intervention where treatment is set to a specific level.
#'
#' @inheritParams cf_sensitivity
#' @param n_thresholds Integer specifying the number of thresholds to evaluate.
#'   Thresholds are evenly spaced between 0 and 1. Default is 201.
#' @param thresholds Optional numeric vector of specific thresholds to use.
#'   If provided, overrides `n_thresholds`.
#' @param include_naive Logical indicating whether to also compute the naive
#'   ROC curve for comparison. Default is TRUE.
#'
#' @return An object of class `c("cf_roc", "roc_curve")` containing:
#'   \item{thresholds}{Thresholds used}
#'   \item{sensitivity}{Sensitivity (TPR) at each threshold}
#'   \item{fpr}{False positive rate at each threshold}
#'   \item{specificity}{Specificity at each threshold}
#'   \item{naive_sensitivity}{Naive sensitivity (if include_naive=TRUE)}
#'   \item{naive_fpr}{Naive FPR (if include_naive=TRUE)}
#'   \item{auc}{Area under the ROC curve (computed via trapezoidal rule)}
#'   \item{naive_auc}{Naive AUC (if include_naive=TRUE)}
#'   \item{estimator}{Estimator used}
#'   \item{n_obs}{Number of observations}
#'
#' @details
#' The ROC curve plots sensitivity (true positive rate) against the false
#' positive rate (1 - specificity) at various classification thresholds.
#'
#' This function computes counterfactual sensitivity and FPR at multiple
#' thresholds using the estimators from [cf_sensitivity()] and [cf_fpr()].
#' The area under the curve (AUC) is computed using the trapezoidal rule on
#' the discrete threshold grid. For exact AUC estimation, use [cf_auc()]
#' which employs the Wilcoxon-Mann-Whitney statistic.
#'
#' @references
#' Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
#' "Counterfactual risk assessments, evaluation, and fairness."
#' *Proceedings of the 2020 Conference on Fairness, Accountability, and
#' Transparency*, 582-593.
#'
#' @seealso [cf_sensitivity()], [cf_specificity()], [cf_fpr()], [plot.cf_roc()], [cf_auc()]
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
#' # Compute counterfactual ROC curve
#' roc <- cf_roc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   n_thresholds = 51
#' )
#' print(roc)
#'
#' # Plot the ROC curve
#' plot(roc)
cf_roc <- function(predictions,
                   outcomes,
                   treatment,
                   covariates,
                   treatment_level = 0,
                   estimator = c("dr", "om", "ipw", "naive"),
                   propensity_model = NULL,
                   outcome_model = NULL,
                   n_thresholds = 201,
                   thresholds = NULL,
                   include_naive = TRUE,
                   ...) {

  estimator <- match.arg(estimator)

  # Validate inputs
  .validate_inputs(predictions, outcomes, treatment, covariates)

  if (!all(outcomes %in% c(0, 1))) {
    stop("ROC curve requires binary outcomes (0/1)")
  }

  # Set up thresholds
  if (is.null(thresholds)) {
    thresholds <- seq(0, 1, length.out = n_thresholds)
  }
  thresholds <- sort(unique(thresholds))

  n <- length(outcomes)

  # Fit nuisance models once
  if (estimator != "naive") {
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

  # Compute sensitivity at all thresholds
  sensitivity <- sapply(thresholds, function(c) {
    .compute_cf_sensitivity(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      threshold = c,
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
  })

  # Compute specificity at all thresholds
  specificity <- sapply(thresholds, function(c) {
    .compute_cf_specificity(
      predictions = predictions,
      outcomes = outcomes,
      treatment = treatment,
      covariates = covariates,
      threshold = c,
      treatment_level = treatment_level,
      estimator = estimator,
      propensity_model = propensity_model,
      outcome_model = outcome_model
    )
  })

  # FPR = 1 - specificity
  fpr <- 1 - specificity

  # Compute AUC using trapezoidal rule
  ord <- order(fpr)
  auc <- .compute_auc_trapezoid(fpr[ord], sensitivity[ord])

  # Naive ROC curve
  naive_sensitivity <- NULL
  naive_fpr <- NULL
  naive_auc <- NULL

  if (include_naive) {
    naive_sensitivity <- sapply(thresholds, function(c) {
      .compute_cf_sens_spec_naive(predictions, outcomes, c, "sensitivity")
    })

    naive_specificity <- sapply(thresholds, function(c) {
      .compute_cf_sens_spec_naive(predictions, outcomes, c, "specificity")
    })

    naive_fpr <- 1 - naive_specificity

    ord_naive <- order(naive_fpr)
    naive_auc <- .compute_auc_trapezoid(naive_fpr[ord_naive], naive_sensitivity[ord_naive])
  }

  result <- list(
    thresholds = thresholds,
    sensitivity = sensitivity,
    fpr = fpr,
    specificity = specificity,
    naive_sensitivity = naive_sensitivity,
    naive_fpr = naive_fpr,
    auc = auc,
    naive_auc = naive_auc,
    estimator = estimator,
    treatment_level = treatment_level,
    n_obs = n,
    include_naive = include_naive
  )

  class(result) <- c("cf_roc", "roc_curve")
  return(result)
}


# ==============================================================================
# Internal Functions
# ==============================================================================

#' Compute AUC using trapezoidal rule
#' @noRd
.compute_auc_trapezoid <- function(fpr, tpr) {
  # Remove NA values
  valid <- !is.na(fpr) & !is.na(tpr)
  fpr <- fpr[valid]
  tpr <- tpr[valid]

  if (length(fpr) < 2) return(NA_real_)

  # Sort by FPR
  ord <- order(fpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]

  # Trapezoidal integration
  n <- length(fpr)
  auc <- sum((fpr[2:n] - fpr[1:(n-1)]) * (tpr[2:n] + tpr[1:(n-1)]) / 2)

  return(auc)
}


# ==============================================================================
# Print Methods
# ==============================================================================

#' @export
print.tr_roc <- function(x, digits = 4, ...) {
  # Determine if factual mode (no treatment) or counterfactual mode
  factual_mode <- is.null(x$treatment_level)
  mode_label <- if (factual_mode) "Factual " else ""

  cat("\n", mode_label, "Transportable ROC Curve\n", sep = "")
  cat("=======================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Analysis:", x$analysis, "\n")
  if (!factual_mode) {
    cat("Treatment level:", x$treatment_level, "\n")
  }
  cat("N (source):", x$n_source, "\n")
  cat("N (target):", x$n_target, "\n")
  cat("Thresholds evaluated:", length(x$thresholds), "\n\n")

  cat("AUC:", round(x$auc, digits), "\n")
  if (x$include_naive) {
    cat("Naive AUC:", round(x$naive_auc, digits), "\n")
  }

  cat("\nUse plot() to visualize the ROC curve.\n")
  cat("\n")
  invisible(x)
}


#' @export
print.cf_roc <- function(x, digits = 4, ...) {
  cat("\nCounterfactual ROC Curve\n")
  cat("========================\n\n")

  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N:", x$n_obs, "\n")
  cat("Thresholds evaluated:", length(x$thresholds), "\n\n")

  cat("AUC:", round(x$auc, digits), "\n")
  if (x$include_naive) {
    cat("Naive AUC:", round(x$naive_auc, digits), "\n")
  }

  cat("\nUse plot() to visualize the ROC curve.\n")
  cat("\n")
  invisible(x)
}


# ==============================================================================
# Plot Methods
# ==============================================================================

#' Plot ROC Curve
#'
#' Creates a plot of the ROC curve.
#'
#' @param x An object of class `tr_roc` or `cf_roc`.
#' @param add_diagonal Logical indicating whether to add a diagonal reference
#'   line (representing random classifier). Default is TRUE.
#' @param show_naive Logical indicating whether to show the naive ROC curve
#'   if available. Default is TRUE.
#' @param main Title for the plot.
#' @param col Color for the ROC curve. Default is "blue".
#' @param naive_col Color for the naive ROC curve. Default is "gray50".
#' @param lwd Line width. Default is 2.
#' @param ... Additional arguments passed to plot().
#'
#' @return Invisibly returns the input object.
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
#' roc <- cf_roc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   n_thresholds = 51
#' )
#' plot(roc)
plot.tr_roc <- function(x, add_diagonal = TRUE, show_naive = TRUE,
                        main = NULL, col = "blue", naive_col = "gray50",
                        lwd = 2, ...) {

  if (is.null(main)) {
    main <- "ROC Curve in Target Population"
  }

  # Sort by FPR for proper plotting
  ord <- order(x$fpr)

  plot(x$fpr[ord], x$sensitivity[ord],
       type = "l", col = col, lwd = lwd,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "Sensitivity (True Positive Rate)",
       main = main, ...)

  # Add naive curve if requested
  if (show_naive && x$include_naive && !is.null(x$naive_fpr)) {
    ord_naive <- order(x$naive_fpr)
    lines(x$naive_fpr[ord_naive], x$naive_sensitivity[ord_naive],
          col = naive_col, lwd = lwd, lty = 2)
  }

  # Add diagonal reference line
  if (add_diagonal) {
    abline(0, 1, col = "gray70", lty = 3)
  }

  # Add legend with estimator name
  est_label <- toupper(x$estimator)
  legend_labels <- paste0(est_label, " (AUC = ", round(x$auc, 4), ")")
  legend_cols <- col
  legend_lty <- 1

  if (show_naive && x$include_naive) {
    legend_labels <- c(legend_labels,
                       paste0("Naive (AUC = ", round(x$naive_auc, 4), ")"))
    legend_cols <- c(legend_cols, naive_col)
    legend_lty <- c(legend_lty, 2)
  }

  legend("bottomright", legend = legend_labels, col = legend_cols,
         lty = legend_lty, lwd = lwd, bty = "n")

  invisible(x)
}


#' @rdname plot.tr_roc
#' @export
plot.cf_roc <- function(x, add_diagonal = TRUE, show_naive = TRUE,
                        main = NULL, col = "blue", naive_col = "gray50",
                        lwd = 2, ...) {

  if (is.null(main)) {
    main <- "Counterfactual ROC Curve"
  }

  # Sort by FPR for proper plotting
  ord <- order(x$fpr)

  plot(x$fpr[ord], x$sensitivity[ord],
       type = "l", col = col, lwd = lwd,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "Sensitivity (True Positive Rate)",
       main = main, ...)

  # Add naive curve if requested
  if (show_naive && x$include_naive && !is.null(x$naive_fpr)) {
    ord_naive <- order(x$naive_fpr)
    lines(x$naive_fpr[ord_naive], x$naive_sensitivity[ord_naive],
          col = naive_col, lwd = lwd, lty = 2)
  }

  # Add diagonal reference line
  if (add_diagonal) {
    abline(0, 1, col = "gray70", lty = 3)
  }

  # Add legend with estimator name
  est_label <- toupper(x$estimator)
  legend_labels <- paste0(est_label, " (AUC = ", round(x$auc, 4), ")")
  legend_cols <- col
  legend_lty <- 1

  if (show_naive && x$include_naive) {
    legend_labels <- c(legend_labels,
                       paste0("Naive (AUC = ", round(x$naive_auc, 4), ")"))
    legend_cols <- c(legend_cols, naive_col)
    legend_lty <- c(legend_lty, 2)
  }

  legend("bottomright", legend = legend_labels, col = legend_cols,
         lty = legend_lty, lwd = lwd, bty = "n")

  invisible(x)
}


# ==============================================================================
# Conversion to Data Frame
# ==============================================================================

#' Convert ROC Curve to Data Frame
#'
#' Converts an ROC curve object to a data frame suitable for use with ggplot2.
#'
#' @param x An object of class `tr_roc` or `cf_roc`.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns:
#'   \item{threshold}{Classification threshold}
#'   \item{fpr}{False positive rate}
#'   \item{sensitivity}{Sensitivity (TPR)}
#'   \item{specificity}{Specificity}
#'   \item{type}{Either "adjusted" or "naive"}
#'
#' @method as.data.frame tr_roc
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' a <- rbinom(n, 1, 0.5)
#' y <- rbinom(n, 1, plogis(-1 + x))
#' pred <- plogis(-1 + 0.8 * x)
#'
#' roc <- cf_roc(
#'   predictions = pred,
#'   outcomes = y,
#'   treatment = a,
#'   covariates = data.frame(x = x),
#'   n_thresholds = 21
#' )
#' df <- as.data.frame(roc)
#' head(df)
as.data.frame.tr_roc <- function(x, ...) {
  df <- data.frame(
    threshold = x$thresholds,
    fpr = x$fpr,
    sensitivity = x$sensitivity,
    specificity = x$specificity,
    type = "adjusted"
  )

  if (x$include_naive && !is.null(x$naive_fpr)) {
    df_naive <- data.frame(
      threshold = x$thresholds,
      fpr = x$naive_fpr,
      sensitivity = x$naive_sensitivity,
      specificity = 1 - x$naive_fpr,
      type = "naive"
    )
    df <- rbind(df, df_naive)
  }

  return(df)
}


#' @rdname as.data.frame.tr_roc
#' @method as.data.frame cf_roc
#' @export
as.data.frame.cf_roc <- function(x, ...) {
  df <- data.frame(
    threshold = x$thresholds,
    fpr = x$fpr,
    sensitivity = x$sensitivity,
    specificity = x$specificity,
    type = "adjusted"
  )

  if (x$include_naive && !is.null(x$naive_fpr)) {
    df_naive <- data.frame(
      threshold = x$thresholds,
      fpr = x$naive_fpr,
      sensitivity = x$naive_sensitivity,
      specificity = 1 - x$naive_fpr,
      type = "naive"
    )
    df <- rbind(df, df_naive)
  }

  return(df)
}
