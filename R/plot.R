#' Plot Method for cf_calibration Objects
#'
#' Creates a calibration plot showing predicted vs observed probabilities
#' under the counterfactual intervention.
#'
#' @param x A `cf_calibration` object.
#' @param add_histogram Logical; add histogram of predictions (default: TRUE).
#' @param add_rug Logical; add rug plot (default: TRUE).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object (if ggplot2 available) or base R plot.
#'
#' @export
plot.cf_calibration <- function(x, add_histogram = TRUE, add_rug = TRUE, ...) {

  estimator_label <- toupper(x$estimator)
  treatment_label <- x$treatment_level
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R plot
    plot(x$predicted, x$observed,
         type = "l", lwd = 2,
         xlab = "Predicted probability",
         ylab = sprintf("Probability under A = %s (%s)", treatment_label, estimator_label),
         main = "Counterfactual Calibration Curve",
         xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, lty = 2, col = "gray")
    return(invisible(NULL))
  }

  # ggplot2 version
  df <- data.frame(
    predicted = x$predicted,
    observed = x$observed
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_line(linewidth = 1.2, color = "#2E86AB") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "gray50") +
    ggplot2::labs(
      x = "Predicted probability",
      y = sprintf("Probability under A = %s (%s)", treatment_label, estimator_label),
      title = "Counterfactual Calibration Curve",
      subtitle = sprintf("ICI = %.3f, Emax = %.3f", x$ici, x$emax)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw()

  return(p)
}


#' Plot Method for cf_auc Objects
#'
#' Creates an ROC curve plot.
#'
#' @param x A `cf_auc` object.
#' @param ... Additional arguments (ignored).
#'
#' @return A plot object.
#'
#' @export
plot.cf_auc <- function(x, ...) {
  # Note: Full ROC curve plotting would require storing more information

  # This is a placeholder that could be expanded

  cat("ROC curve plotting requires additional data storage.\n")
  cat("Counterfactual AUC:", round(x$estimate, 3), "\n")
  cat("Naive AUC:", round(x$naive_estimate, 3), "\n")

  invisible(NULL)
}
