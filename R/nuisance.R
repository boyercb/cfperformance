#' Fit Nuisance Models for Counterfactual Performance Estimation
#'
#' Fits propensity score and outcome models for use in counterfactual
#' performance estimation.
#'
#' @param propensity_formula Formula for the propensity score model.
#' @param outcome_formula Formula for the outcome model.
#' @param data Data frame containing the variables.
#' @param treatment_level The counterfactual treatment level (default: 0).
#' @param propensity_method Method for propensity score estimation:
#'   - `"glm"`: Logistic regression (default)
#'   - `"gam"`: Generalized additive model (requires mgcv)
#' @param outcome_method Method for outcome model estimation:
#'   - `"glm"`: Logistic regression (default)
#'   - `"gam"`: Generalized additive model (requires mgcv)
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return An object of class `cf_nuisance` containing:
#'   \item{propensity}{Fitted propensity score model}
#'   \item{outcome}{Fitted outcome model}
#'   \item{treatment_level}{Counterfactual treatment level}
#'
#' @details
#' This function provides a convenient way to fit nuisance models that can
#' then be passed to [cf_mse()], [cf_auc()], or [cf_calibration()].
#'
#' For the outcome model, only observations with the counterfactual treatment
#' level are used for fitting, as required for the conditional loss estimator.
#'
#' @seealso [cf_mse()], [cf_auc()], [cf_calibration()]
#'
#' @export
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' df <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' df$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * df$x1))
#' df$y <- rbinom(n, 1, plogis(-1 + df$x1 + 0.5 * df$x2 - 0.5 * df$a))
#'
#' # Fit nuisance models
#' nuisance <- fit_nuisance(
#'   propensity_formula = a ~ x1 + x2,
#'   outcome_formula = y ~ x1 + x2,
#'   data = df,
#'   treatment_level = 0
#' )
#'
#' print(nuisance)
fit_nuisance <- function(propensity_formula,
                         outcome_formula,
                         data,
                         treatment_level = 0,
                         propensity_method = c("glm", "gam"),
                         outcome_method = c("glm", "gam"),
                         ...) {

  propensity_method <- match.arg(propensity_method)
  outcome_method <- match.arg(outcome_method)

  # Get treatment variable name from formula
  treatment_var <- as.character(propensity_formula[[2]])
  outcome_var <- as.character(outcome_formula[[2]])

  # Fit propensity score model
  if (propensity_method == "glm") {
    propensity_model <- glm(propensity_formula, data = data, family = binomial())
  } else if (propensity_method == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' required for GAM. Install it or use method='glm'.")
    }
    propensity_model <- mgcv::gam(propensity_formula, data = data,
                                   family = binomial())
  }

  # Subset data for outcome model (only those with counterfactual treatment)
  subset_idx <- data[[treatment_var]] == treatment_level
  outcome_data <- data[subset_idx, ]

  # Fit outcome model
  if (outcome_method == "glm") {
    outcome_model <- glm(outcome_formula, data = outcome_data, family = binomial())
  } else if (outcome_method == "gam") {
    outcome_model <- mgcv::gam(outcome_formula, data = outcome_data,
                                family = binomial())
  }

  # Store full data reference for prediction
  attr(outcome_model, "full_data") <- data

  result <- list(
    propensity = propensity_model,
    outcome = outcome_model,
    treatment_level = treatment_level,
    treatment_var = treatment_var,
    outcome_var = outcome_var,
    propensity_method = propensity_method,
    outcome_method = outcome_method
  )

  class(result) <- "cf_nuisance"
  return(result)
}


#' Print Method for cf_nuisance Objects
#'
#' @param x A `cf_nuisance` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.cf_nuisance <- function(x, ...) {
  cat("\n")
  cat("Nuisance Models for Counterfactual Performance\n")
  cat(paste(rep("-", 45), collapse = ""), "\n")
  cat("Treatment level:", x$treatment_level, "\n\n")

  cat("Propensity Score Model:\n")
  cat("  Method:", x$propensity_method, "\n")
  cat("  Formula:", deparse(formula(x$propensity)), "\n\n")

  cat("Outcome Model:\n")
  cat("  Method:", x$outcome_method, "\n")
  cat("  Formula:", deparse(formula(x$outcome)), "\n")
  cat("  Fitted on n =", nobs(x$outcome), "observations with A =",
      x$treatment_level, "\n")

  cat("\n")
  invisible(x)
}
