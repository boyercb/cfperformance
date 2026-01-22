#' Print Method for cf_performance Objects
#'
#' @param x A `cf_performance` object.
#' @param digits Number of digits to print (default: 4).
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.cf_performance <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Counterfactual", toupper(x$metric), "Estimation\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  cat("Estimator:", x$estimator, "\n")
  cat("Treatment level:", x$treatment_level, "\n")
  cat("N observations:", x$n_obs, "\n")
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


#' Summary Method for cf_performance Objects
#'
#' @param object A `cf_performance` object.
#' @param ... Additional arguments (ignored).
#'
#' @return A summary list with key statistics.
#'
#' @export
summary.cf_performance <- function(object, ...) {
  cat("\n")
  cat("Summary: Counterfactual", toupper(object$metric), "Estimation\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Settings:\n")
  cat("  Estimator:", object$estimator, "\n")
  cat("  Treatment level:", object$treatment_level, "\n")
  cat("  Sample size:", object$n_obs, "\n")
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
      Estimator = c("Counterfactual", "Naive"),
      Estimate = c(object$estimate, object$naive_estimate),
      SE = c(object$se, NA),
      CI_lower = c(object$ci_lower, NA),
      CI_upper = c(object$ci_upper, NA)
    )
    print(results_df, digits = 4, row.names = FALSE)

    if (!is.null(object$naive_estimate)) {
      diff <- object$estimate - object$naive_estimate
      cat("\nDifference (Counterfactual - Naive):", round(diff, 4), "\n")
    }
  }

  cat("\n")
  invisible(object)
}


#' Coefficient Method for cf_performance Objects
#'
#' @param object A `cf_performance` object.
#' @param ... Additional arguments (ignored).
#'
#' @return The point estimate.
#'
#' @export
coef.cf_performance <- function(object, ...) {
  if (object$metric == "calibration") {
    return(c(ici = object$ici, e50 = object$e50,
             e90 = object$e90, emax = object$emax))
  }
  object$estimate
}


#' Confidence Interval Method for cf_performance Objects
#'
#' @param object A `cf_performance` object.
#' @param parm Parameter (ignored, only one parameter).
#' @param level Confidence level (default uses object's conf_level).
#' @param ... Additional arguments (ignored).
#'
#' @return A matrix with confidence interval bounds.
#'
#' @export
confint.cf_performance <- function(object, parm = NULL, level = NULL, ...) {
  if (object$metric == "calibration") {
    warning("Confidence intervals not available for calibration")
    return(NULL)
  }

  if (is.null(object$ci_lower) || is.na(object$ci_lower)) {
    warning("Confidence intervals not computed. Set se_method != 'none'.")
    return(NULL)
  }

  # Use object's confidence level if not specified
  level <- level %||% object$conf_level %||% 0.95

  ci <- matrix(c(object$ci_lower, object$ci_upper), nrow = 1)
  colnames(ci) <- c(
    paste0((1 - level) / 2 * 100, "%"),
    paste0((1 - (1 - level) / 2) * 100, "%")
  )
  rownames(ci) <- object$metric
  ci
}
