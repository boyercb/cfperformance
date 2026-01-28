#' Print Method for tr_performance Objects
#'
#' @param x A `tr_performance` object.
#' @param digits Number of digits to print (default: 4).
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.tr_performance <- function(x, digits = 4, ...) {

  # Determine if factual mode (no treatment) or counterfactual mode
  factual_mode <- is.null(x$treatment_level)
  mode_label <- if (factual_mode) "Factual" else "Counterfactual"

  cat("\n")
  cat(mode_label, "Transportable", toupper(x$metric), "Estimation\n")
  cat(paste(rep("-", 45), collapse = ""), "\n")

  cat("Analysis:", x$analysis, "\n")
  cat("Estimator:", x$estimator, "\n")
  if (!factual_mode) {
    cat("Treatment level:", x$treatment_level, "\n")
  }
  cat("N target:", x$n_target, " | N source:", x$n_source, "\n")
  cat("\n")

  if (x$metric == "calibration") {
    cat("Calibration Metrics:\n")
    cat("  ICI (Integrated Calibration Index):", round(x$ici, digits), "\n")
    cat("  E50 (Median absolute error):", round(x$e50, digits), "\n")
    cat("  E90 (90th percentile error):", round(x$e90, digits), "\n")
    cat("  Emax (Maximum error):", round(x$emax, digits), "\n")
  } else {
    cat("Estimate:", round(x$estimate, digits))

    if (!is.null(x$se) && !is.na(x$se)) {
      cat(" (SE:", round(x$se, digits), ")")
    }
    cat("\n")

    if (!is.null(x$ci_lower) && !is.na(x$ci_lower)) {
      cat(sprintf("%d%% CI: [%s, %s]\n",
                  round(x$conf_level * 100),
                  round(x$ci_lower, digits),
                  round(x$ci_upper, digits)))
    }

    if (!is.null(x$naive_estimate)) {
      cat("\nNaive estimate:", round(x$naive_estimate, digits), "\n")
    }
  }

  cat("\n")
  invisible(x)
}


#' Summary Method for tr_performance Objects
#'
#' @param object A `tr_performance` object.
#' @param ... Additional arguments (ignored).
#'
#' @return A summary list with key statistics.
#'
#' @export
summary.tr_performance <- function(object, ...) {
  # Determine if factual mode (no treatment) or counterfactual mode
  factual_mode <- is.null(object$treatment_level)
  mode_label <- if (factual_mode) "Factual" else "Counterfactual"

  cat("\n")
  cat("Summary:", mode_label, "Transportable", toupper(object$metric), "Estimation\n")
  cat(paste(rep("=", 55), collapse = ""), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Settings:\n")
  cat("  Mode:", mode_label, "\n")
  cat("  Analysis type:", object$analysis, "\n")
  cat("  Estimator:", object$estimator, "\n")
  if (!factual_mode) {
    cat("  Treatment level:", object$treatment_level, "\n")
  }
  cat("  Target sample size:", object$n_target, "\n")
  cat("  Source sample size:", object$n_source, "\n")
  cat("\n")

  if (object$metric == "calibration") {
    cat("Calibration Statistics:\n")
    cat("  Smoother:", object$smoother, "\n")
    cat("  ICI:", round(object$ici, 4), "\n")
    cat("  E50:", round(object$e50, 4), "\n")
    cat("  E90:", round(object$e90, 4), "\n")
    cat("  Emax:", round(object$emax, 4), "\n")
  } else {
    cat("Results:\n")

    results_df <- data.frame(
      Estimator = c("Transportable", "Naive"),
      Estimate = c(object$estimate, object$naive_estimate),
      SE = c(object$se, NA),
      CI_lower = c(object$ci_lower, NA),
      CI_upper = c(object$ci_upper, NA)
    )
    print(results_df, digits = 4, row.names = FALSE)

    if (!is.null(object$naive_estimate)) {
      diff <- object$estimate - object$naive_estimate
      cat("\nDifference (Transportable - Naive):", round(diff, 4), "\n")
    }
  }

  cat("\n")
  invisible(object)
}


#' Coefficient Method for tr_performance Objects
#'
#' @param object A `tr_performance` object.
#' @param ... Additional arguments (ignored).
#'
#' @return The point estimate.
#'
#' @export
coef.tr_performance <- function(object, ...) {
  if (object$metric == "calibration") {
    return(c(ici = object$ici, e50 = object$e50,
             e90 = object$e90, emax = object$emax))
  }
  object$estimate
}


#' Confidence Interval Method for tr_performance Objects
#'
#' @param object A `tr_performance` object.
#' @param parm Parameter (ignored, only one parameter).
#' @param level Confidence level (default uses object's conf_level).
#' @param ... Additional arguments (ignored).
#'
#' @return A matrix with confidence interval bounds.
#'
#' @export
confint.tr_performance <- function(object, parm = NULL, level = NULL, ...) {
  # Use object's confidence level if not specified
  level <- level %||% object$conf_level %||% 0.95

  if (object$metric == "calibration") {
    # Calibration has multiple metrics
    if (is.null(object$ci_lower)) {
      warning("Confidence intervals not computed. Set se_method = 'bootstrap'.")
      return(NULL)
    }

    ci <- matrix(
      c(object$ci_lower$ici, object$ci_upper$ici,
        object$ci_lower$e50, object$ci_upper$e50,
        object$ci_lower$e90, object$ci_upper$e90,
        object$ci_lower$emax, object$ci_upper$emax),
      nrow = 4, ncol = 2, byrow = TRUE
    )
    colnames(ci) <- c(
      paste0((1 - level) / 2 * 100, "%"),
      paste0((1 - (1 - level) / 2) * 100, "%")
    )
    rownames(ci) <- c("ici", "e50", "e90", "emax")
    return(ci)
  }

  if (is.null(object$ci_lower) || is.na(object$ci_lower)) {
    warning("Confidence intervals not computed. Set se_method != 'none'.")
    return(NULL)
  }

  ci <- matrix(c(object$ci_lower, object$ci_upper), nrow = 1)
  colnames(ci) <- c(
    paste0((1 - level) / 2 * 100, "%"),
    paste0((1 - (1 - level) / 2) * 100, "%")
  )
  rownames(ci) <- paste0("tr_", object$metric)
  ci
}
