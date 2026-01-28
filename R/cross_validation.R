#' Cross-Validation with Counterfactual Performance Metrics
#'
#' Performs K-fold cross-validation to estimate out-of-sample counterfactual
#' model performance. This function trains and evaluates prediction models
#' while properly accounting for treatment effects.
#'
#' @param formula A formula specifying the prediction model (e.g., `Y ~ X1 + X2`).
#' @param data A data frame containing the variables in the formula plus
#'   `treatment` and any additional covariates for nuisance models.
#' @param treatment Character string naming the treatment variable in `data`.
#' @param treatment_level The counterfactual treatment level (default: 0).
#' @param nuisance_covariates Character vector of covariate names for nuisance
#'   models. If NULL, uses all predictors from the formula.
#' @param metric Character string specifying the performance metric:
#'   `"mse"` (default), `"auc"`, or `"both"`.
#' @param estimator Character string specifying the estimator:
#'   `"dr"` (default), `"cl"` (conditional loss, for MSE), `"om"` (outcome model, for AUC),
#'   `"ipw"`, or `"naive"`. Note: `"cl"` is automatically mapped to `"om"` for AUC metrics.
#' @param K Number of folds (default: 5).
#' @param repeats Number of times to repeat K-fold CV (default: 1).
#' @param stratify Logical indicating whether to stratify folds by outcome
#'   (default: TRUE for binary outcomes).
#' @param seed Random seed for reproducibility (default: NULL).
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class `cf_cv` containing:
#'   \item{results}{Data frame with fold-level performance estimates}
#'   \item{summary}{Summary statistics across folds}
#'   \item{metric}{Performance metric used}
#'   \item{estimator}{Estimator used}
#'   \item{K}{Number of folds}
#'   \item{repeats}{Number of repeats}
#'   \item{call}{The matched call}
#'
#' @details
#' Cross-validation for counterfactual prediction models requires special care:
#'
#' 1. **Nuisance model estimation**: Propensity and outcome models are re-fit
#'    in each training fold to avoid overfitting.
#'
#' 2. **Sample splitting**: The prediction model is trained on the training
#'    fold and evaluated on the test fold using counterfactual estimators.
#'
#' 3. **Stratification**: For binary outcomes, stratified sampling ensures
#'    each fold has similar outcome prevalence.
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
#' "Estimating and evaluating counterfactual prediction models."
#' *Statistics in Medicine*, 44(23-24), e70287. \doi{10.1002/sim.70287}
#'
#' @seealso [cf_mse()], [cf_auc()], [cf_compare()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 300
#' data <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
#' data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))
#'
#' # 5-fold cross-validation
#' cv_result <- cf_cv(
#'   formula = y ~ x1 + x2,
#'   data = data,
#'   treatment = "a",
#'   treatment_level = 0,
#'   metric = "mse",
#'   K = 5
#' )
#' print(cv_result)
cf_cv <- function(formula,
                  data,
                  treatment,
                  treatment_level = 0,
                  nuisance_covariates = NULL,
                  metric = c("mse", "auc", "both"),
                  estimator = c("dr", "cl", "om", "ipw", "naive"),
                  K = 5,
                  repeats = 1,
                  stratify = TRUE,
                  seed = NULL,
                  ...) {

  metric <- match.arg(metric)
  estimator <- match.arg(estimator)

  # Extract variables from formula
  outcome_var <- all.vars(formula)[1]
  predictor_vars <- all.vars(formula)[-1]

  # Validate inputs
  if (!outcome_var %in% names(data)) {
    stop("Outcome variable '", outcome_var, "' not found in data")
  }
  if (!treatment %in% names(data)) {
    stop("Treatment variable '", treatment, "' not found in data")
  }
  if (!all(predictor_vars %in% names(data))) {
    missing <- predictor_vars[!predictor_vars %in% names(data)]
    stop("Predictor variables not found in data: ", paste(missing, collapse = ", "))
  }

  # Set nuisance covariates
  if (is.null(nuisance_covariates)) {
    nuisance_covariates <- predictor_vars
  }

  # Check for binary outcome if AUC requested
  outcomes <- data[[outcome_var]]
  if (metric %in% c("auc", "both") && !all(outcomes %in% c(0, 1))) {
    stop("AUC requires binary outcome (0/1)")
  }

  n <- nrow(data)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Storage for results
  all_results <- list()

  for (rep in seq_len(repeats)) {
    # Create fold assignments
    if (stratify && all(outcomes %in% c(0, 1))) {
      folds <- .stratified_folds(outcomes, K)
    } else {
      folds <- sample(rep(1:K, length.out = n))
    }

    for (k in 1:K) {
      # Split data
      train_idx <- which(folds != k)
      test_idx <- which(folds == k)

      train_data <- data[train_idx, , drop = FALSE]
      test_data <- data[test_idx, , drop = FALSE]

      # Fit prediction model on training data
      pred_model <- glm(formula, data = train_data, family = binomial())

      # Get predictions on test data
      predictions <- predict(pred_model, newdata = test_data, type = "response")

      # Extract test data components
      test_outcomes <- test_data[[outcome_var]]
      test_treatment <- test_data[[treatment]]
      test_covariates <- test_data[, nuisance_covariates, drop = FALSE]

      # Compute performance metrics
      result_row <- data.frame(
        repeat_id = rep,
        fold = k,
        n_train = length(train_idx),
        n_test = length(test_idx)
      )

      if (metric %in% c("mse", "both")) {
        mse_result <- cf_mse(
          predictions = predictions,
          outcomes = test_outcomes,
          treatment = test_treatment,
          covariates = test_covariates,
          treatment_level = treatment_level,
          estimator = estimator,
          se_method = "none"
        )
        result_row$mse <- mse_result$estimate
        result_row$mse_naive <- mse_result$naive_estimate
      }

      if (metric %in% c("auc", "both")) {
        # Map "cl" to "om" for AUC (cf_auc uses "om" for outcome model estimator)
        auc_estimator <- if (estimator == "cl") "om" else estimator
        auc_result <- cf_auc(
          predictions = predictions,
          outcomes = test_outcomes,
          treatment = test_treatment,
          covariates = test_covariates,
          treatment_level = treatment_level,
          estimator = auc_estimator,
          se_method = "none"
        )
        result_row$auc <- auc_result$estimate
        result_row$auc_naive <- auc_result$naive_estimate
      }

      all_results[[length(all_results) + 1]] <- result_row
    }
  }

  # Combine results
  results_df <- do.call(rbind, all_results)

  # Compute summary statistics
  summary_stats <- .summarize_cv_results(results_df, metric)

  # Construct return object
  result <- list(
    results = results_df,
    summary = summary_stats,
    metric = metric,
    estimator = estimator,
    K = K,
    repeats = repeats,
    formula = formula,
    treatment = treatment,
    treatment_level = treatment_level,
    call = match.call()
  )

  class(result) <- "cf_cv"
  return(result)
}


#' Compare Multiple Prediction Models
#'
#' Compares the counterfactual performance of multiple prediction models
#' using cross-validation or a held-out test set.
#'
#' @param models A named list of model formulas or fitted model objects.
#' @param data A data frame containing all variables.
#' @param treatment Character string naming the treatment variable.
#' @param treatment_level The counterfactual treatment level (default: 0).
#' @param nuisance_covariates Character vector of covariate names for nuisance
#'   models. If NULL, inferred from model formulas.
#' @param metric Character string specifying the performance metric.
#' @param estimator Character string specifying the estimator.
#' @param method Comparison method: `"cv"` for cross-validation (default),
#'   or `"holdout"` for train/test split.
#' @param K Number of folds for CV (default: 5).
#' @param test_prop Proportion of data for test set if method = "holdout".
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to cf_cv or internal functions.
#'
#' @return An object of class `cf_compare` containing:
#'   \item{results}{Data frame with model performance comparisons}
#'   \item{best_model}{Name of the best performing model}
#'   \item{metric}{Performance metric used}
#'   \item{method}{Comparison method used}
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 300
#' data <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = rnorm(n)
#' )
#' data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
#' data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))
#'
#' # Compare models
#' models <- list(
#'   "Simple" = y ~ x1,
#'   "Full" = y ~ x1 + x2 + x3
#' )
#'
#' comparison <- cf_compare(
#'   models = models,
#'   data = data,
#'   treatment = "a",
#'   metric = "mse",
#'   K = 3
#' )
#' print(comparison)
cf_compare <- function(models,
                       data,
                       treatment,
                       treatment_level = 0,
                       nuisance_covariates = NULL,
                       metric = c("mse", "auc", "both"),
                       estimator = c("dr", "cl", "om", "ipw", "naive"),
                       method = c("cv", "holdout"),
                       K = 5,
                       test_prop = 0.2,
                       seed = NULL,
                       ...) {

  metric <- match.arg(metric)
  estimator <- match.arg(estimator)
  method <- match.arg(method)

  if (!is.list(models)) {
    stop("'models' must be a list of formulas or fitted models")
  }

  if (is.null(names(models))) {
    names(models) <- paste0("Model_", seq_along(models))
  }

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Infer nuisance covariates from all models if not specified
  if (is.null(nuisance_covariates)) {
    all_vars <- unique(unlist(lapply(models, function(m) {
      if (inherits(m, "formula")) {
        all.vars(m)[-1]
      } else {
        names(m$model)[-1]
      }
    })))
    nuisance_covariates <- all_vars
  }

  results_list <- list()

  for (model_name in names(models)) {
    model <- models[[model_name]]

    # Convert to formula if needed
    if (!inherits(model, "formula")) {
      if (inherits(model, "glm") || inherits(model, "lm")) {
        model <- formula(model)
      } else {
        stop("Model '", model_name, "' must be a formula or fitted glm/lm object")
      }
    }

    if (method == "cv") {
      # Cross-validation
      cv_result <- cf_cv(
        formula = model,
        data = data,
        treatment = treatment,
        treatment_level = treatment_level,
        nuisance_covariates = nuisance_covariates,
        metric = metric,
        estimator = estimator,
        K = K,
        seed = seed,
        ...
      )

      model_result <- data.frame(
        model = model_name,
        stringsAsFactors = FALSE
      )

      if (metric %in% c("mse", "both")) {
        model_result$mse_mean <- cv_result$summary$mse_mean
        model_result$mse_se <- cv_result$summary$mse_se
        model_result$mse_naive_mean <- cv_result$summary$mse_naive_mean
      }

      if (metric %in% c("auc", "both")) {
        model_result$auc_mean <- cv_result$summary$auc_mean
        model_result$auc_se <- cv_result$summary$auc_se
        model_result$auc_naive_mean <- cv_result$summary$auc_naive_mean
      }

    } else {
      # Holdout method
      n <- nrow(data)
      test_idx <- sample(n, floor(n * test_prop))
      train_idx <- setdiff(1:n, test_idx)

      train_data <- data[train_idx, , drop = FALSE]
      test_data <- data[test_idx, , drop = FALSE]

      # Fit model
      pred_model <- glm(model, data = train_data, family = binomial())
      predictions <- predict(pred_model, newdata = test_data, type = "response")

      # Get outcome variable name
      outcome_var <- all.vars(model)[1]

      model_result <- data.frame(
        model = model_name,
        stringsAsFactors = FALSE
      )

      test_covariates <- test_data[, nuisance_covariates, drop = FALSE]

      if (metric %in% c("mse", "both")) {
        mse_result <- cf_mse(
          predictions = predictions,
          outcomes = test_data[[outcome_var]],
          treatment = test_data[[treatment]],
          covariates = test_covariates,
          treatment_level = treatment_level,
          estimator = estimator,
          se_method = "influence"
        )
        model_result$mse_mean <- mse_result$estimate
        model_result$mse_se <- mse_result$se
        model_result$mse_naive_mean <- mse_result$naive_estimate
      }

      if (metric %in% c("auc", "both")) {
        # Map "cl" to "om" for AUC (cf_auc uses "om" for outcome model estimator)
        auc_estimator <- if (estimator == "cl") "om" else estimator
        auc_result <- cf_auc(
          predictions = predictions,
          outcomes = test_data[[outcome_var]],
          treatment = test_data[[treatment]],
          covariates = test_covariates,
          treatment_level = treatment_level,
          estimator = auc_estimator,
          se_method = "influence"
        )
        model_result$auc_mean <- auc_result$estimate
        model_result$auc_se <- auc_result$se
        model_result$auc_naive_mean <- auc_result$naive_estimate
      }
    }

    results_list[[model_name]] <- model_result
  }

  # Combine results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL


  # Determine best model
  if (metric == "mse") {
    best_idx <- which.min(results_df$mse_mean)
  } else if (metric == "auc") {
    best_idx <- which.max(results_df$auc_mean)
  } else {
    # For "both", prefer MSE
    best_idx <- which.min(results_df$mse_mean)
  }
  best_model <- results_df$model[best_idx]

  result <- list(
    results = results_df,
    best_model = best_model,
    metric = metric,
    estimator = estimator,
    method = method,
    K = if (method == "cv") K else NULL,
    test_prop = if (method == "holdout") test_prop else NULL,
    call = match.call()
  )

  class(result) <- "cf_compare"
  return(result)
}


# ============================================================================
# Internal helper functions
# ============================================================================

#' Create Stratified Folds
#'
#' @param y Binary outcome vector
#' @param K Number of folds
#'
#' @return Vector of fold assignments
#' @keywords internal
.stratified_folds <- function(y, K) {
  n <- length(y)
  folds <- integer(n)

  # Split by outcome class
  idx_1 <- which(y == 1)
  idx_0 <- which(y == 0)

  # Assign folds within each class
  folds[idx_1] <- sample(rep(1:K, length.out = length(idx_1)))
  folds[idx_0] <- sample(rep(1:K, length.out = length(idx_0)))

  return(folds)
}


#' Summarize CV Results
#'
#' @param results_df Data frame of fold-level results
#' @param metric Performance metric
#'
#' @return Data frame with summary statistics
#' @keywords internal
.summarize_cv_results <- function(results_df, metric) {
  summary_stats <- list()

  if ("mse" %in% names(results_df)) {
    summary_stats$mse_mean <- mean(results_df$mse, na.rm = TRUE)
    summary_stats$mse_se <- sd(results_df$mse, na.rm = TRUE) / sqrt(nrow(results_df))
    summary_stats$mse_min <- min(results_df$mse, na.rm = TRUE)
    summary_stats$mse_max <- max(results_df$mse, na.rm = TRUE)
    summary_stats$mse_naive_mean <- mean(results_df$mse_naive, na.rm = TRUE)
  }

  if ("auc" %in% names(results_df)) {
    summary_stats$auc_mean <- mean(results_df$auc, na.rm = TRUE)
    summary_stats$auc_se <- sd(results_df$auc, na.rm = TRUE) / sqrt(nrow(results_df))
    summary_stats$auc_min <- min(results_df$auc, na.rm = TRUE)
    summary_stats$auc_max <- max(results_df$auc, na.rm = TRUE)
    summary_stats$auc_naive_mean <- mean(results_df$auc_naive, na.rm = TRUE)
  }

  as.data.frame(summary_stats)
}


# ============================================================================
# Print and Summary Methods
# ============================================================================

#' Print Method for cf_cv Objects
#'
#' @param x A `cf_cv` object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#' @export
print.cf_cv <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Counterfactual Cross-Validation Results\n")
  cat(paste(rep("-", 45), collapse = ""), "\n")

  cat("Folds:", x$K)
  if (x$repeats > 1) {
    cat(" (repeated", x$repeats, "times)")
  }
  cat("\n")
  cat("Estimator:", x$estimator, "\n")
  cat("Treatment level:", x$treatment_level, "\n\n")

  if ("mse_mean" %in% names(x$summary)) {
    cat("MSE:\n")
    cat("  Counterfactual:", round(x$summary$mse_mean, digits),
        "(SE:", round(x$summary$mse_se, digits), ")\n")
    cat("  Naive:         ", round(x$summary$mse_naive_mean, digits), "\n")
  }

  if ("auc_mean" %in% names(x$summary)) {
    if ("mse_mean" %in% names(x$summary)) cat("\n")
    cat("AUC:\n")
    cat("  Counterfactual:", round(x$summary$auc_mean, digits),
        "(SE:", round(x$summary$auc_se, digits), ")\n")
    cat("  Naive:         ", round(x$summary$auc_naive_mean, digits), "\n")
  }

  cat("\n")
  invisible(x)
}


#' Summary Method for cf_cv Objects
#'
#' @param object A `cf_cv` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a summary data frame.
#' @export
summary.cf_cv <- function(object, ...) {
  cat("\n")
  cat("Summary: Counterfactual Cross-Validation\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Settings:\n")
  cat("  Formula:", deparse(object$formula), "\n")
  cat("  K-folds:", object$K, "\n")
  cat("  Repeats:", object$repeats, "\n")
  cat("  Estimator:", object$estimator, "\n")
  cat("  Treatment level:", object$treatment_level, "\n\n")

  cat("Fold-level Results:\n")
  print(object$results, digits = 4, row.names = FALSE)
  cat("\n")

  cat("Summary Statistics:\n")
  print(object$summary, digits = 4, row.names = FALSE)
  cat("\n")

  invisible(object$summary)
}


#' Print Method for cf_compare Objects
#'
#' @param x A `cf_compare` object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#' @export
print.cf_compare <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Counterfactual Model Comparison\n")
  cat(paste(rep("-", 45), collapse = ""), "\n")

  cat("Method:", x$method)
  if (x$method == "cv") {
    cat(" (K =", x$K, ")\n")
  } else {
    cat(" (test_prop =", x$test_prop, ")\n")
  }
  cat("Estimator:", x$estimator, "\n\n")

  # Format results table
  results_print <- x$results

  if ("mse_mean" %in% names(results_print)) {
    results_print$mse_mean <- round(results_print$mse_mean, digits)
    if ("mse_se" %in% names(results_print)) {
      results_print$mse_se <- round(results_print$mse_se, digits)
    }
  }

  if ("auc_mean" %in% names(results_print)) {
    results_print$auc_mean <- round(results_print$auc_mean, digits)
    if ("auc_se" %in% names(results_print)) {
      results_print$auc_se <- round(results_print$auc_se, digits)
    }
  }

  print(results_print, row.names = FALSE)
  cat("\nBest model:", x$best_model, "\n\n")

  invisible(x)
}


#' Summary Method for cf_compare Objects
#'
#' @param object A `cf_compare` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the results data frame.
#' @export
summary.cf_compare <- function(object, ...) {
  cat("\n")
  cat("Summary: Counterfactual Model Comparison\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Settings:\n")
  cat("  Method:", object$method, "\n")
  cat("  Metric:", object$metric, "\n")
  cat("  Estimator:", object$estimator, "\n\n")

  cat("Model Performance:\n")
  print(object$results, digits = 4, row.names = FALSE)
  cat("\n")

  cat("Best Model: ", object$best_model, "\n")

  if (object$metric %in% c("mse", "both")) {
    best_mse <- object$results$mse_mean[object$results$model == object$best_model]
    cat("  MSE:", round(best_mse, 4), "\n")
  }

  if (object$metric %in% c("auc", "both")) {
    best_auc <- object$results$auc_mean[object$results$model == object$best_model]
    cat("  AUC:", round(best_auc, 4), "\n")
  }

  cat("\n")
  invisible(object$results)
}
