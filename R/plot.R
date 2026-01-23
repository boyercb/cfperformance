#' Plot Method for cf_calibration Objects
#'
#' Creates a calibration plot showing predicted vs observed probabilities
#' under the counterfactual intervention, with an optional histogram showing
#' the distribution of predicted probabilities and optional confidence bands.
#'
#' @param x A `cf_calibration` object.
#' @param add_reference Logical; add 45-degree reference line (default: TRUE).
#' @param show_metrics Logical; show calibration metrics on plot (default: TRUE).
#' @param add_histogram Logical; add histogram of predictions below the
#'   calibration curve (default: TRUE).
#' @param add_rug Logical; add rug plot to show individual predictions
#'   (default: FALSE).
#' @param add_ci Logical; add bootstrap confidence bands if available
#'   (default: TRUE if bootstrap was run).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot object (if ggplot2 available) or base R plot.
#'
#' @export
plot.cf_calibration <- function(x, add_reference = TRUE, show_metrics = TRUE,
                                 add_histogram = TRUE, add_rug = FALSE, 
                                 add_ci = !is.null(x$boot_curves), ...) {

  estimator_label <- toupper(x$estimator)
  treatment_label <- x$treatment_level
  has_ci <- !is.null(x$boot_curves) && add_ci
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R plot - simple version without histogram
    if (add_histogram) {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
      par(mar = c(0, 4, 3, 2))
    }
    
    plot(x$predicted, x$observed,
         type = "l", lwd = 2, col = "#2E86AB",
         xlab = if (add_histogram) "" else "Predicted probability",
         ylab = sprintf("Probability under A = %s (%s)", treatment_label, estimator_label),
         main = "Counterfactual Calibration Curve",
         xlim = c(0, 1), ylim = c(0, 1),
         xaxt = if (add_histogram) "n" else "s")
    if (add_reference) {
      abline(0, 1, lty = 2, col = "gray")
    }
    
    # Add CI band if available
    if (has_ci) {
      polygon(c(x$boot_curves$predicted, rev(x$boot_curves$predicted)),
              c(x$boot_curves$ci_lower, rev(x$boot_curves$ci_upper)),
              col = rgb(46/255, 134/255, 171/255, 0.2), border = NA)
    }
    
    if (add_histogram && !is.null(x$predictions_raw)) {
      par(mar = c(4, 4, 0, 2))
      hist(x$predictions_raw, breaks = 30, col = "#2E86AB", border = "white",
           main = "", xlab = "Predicted probability", xlim = c(0, 1))
    }
    
    return(invisible(NULL))
  }

  # ggplot2 version
  df_curve <- data.frame(
    predicted = x$predicted,
    observed = x$observed
  )

  # Build subtitle
  subtitle_text <- NULL
  if (show_metrics) {
    subtitle_text <- sprintf("ICI = %.3f, Emax = %.3f%s", x$ici, x$emax,
                             if (has_ci) sprintf(" (%d%% CI)", round(x$conf_level * 100)) else "")
  }

  # Main calibration curve
  p_cal <- ggplot2::ggplot(df_curve, ggplot2::aes(x = .data$predicted, y = .data$observed))
  
  # Add confidence band first (so it's behind the line)
  if (has_ci) {
    df_ci <- x$boot_curves
    p_cal <- p_cal +
      ggplot2::geom_ribbon(
        data = df_ci,
        ggplot2::aes(x = .data$predicted, ymin = .data$ci_lower, ymax = .data$ci_upper),
        fill = "#2E86AB", alpha = 0.2, inherit.aes = FALSE
      )
  }
  
  p_cal <- p_cal +
    ggplot2::geom_line(linewidth = 1.2, color = "#2E86AB") +
    ggplot2::labs(
      y = sprintf("Probability under A = %s (%s)", treatment_label, estimator_label),
      title = "Counterfactual Calibration Curve",
      subtitle = subtitle_text
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw()
  
  if (add_reference) {
    p_cal <- p_cal + ggplot2::geom_abline(slope = 1, intercept = 0,
                                           linetype = "dashed", color = "gray50")
  }
  
  # Add rug if requested
  if (add_rug && !is.null(x$predictions_raw)) {
    df_raw <- data.frame(predictions_raw = x$predictions_raw)
    p_cal <- p_cal + 
      ggplot2::geom_rug(data = df_raw, 
                        ggplot2::aes(x = .data$predictions_raw, y = NULL),
                        alpha = 0.3, color = "#2E86AB")
  }
  
  # Add histogram subplot if requested
  if (add_histogram && !is.null(x$predictions_raw)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      # Fall back to just the calibration plot with rug
      message("Install 'patchwork' package for histogram subplot. Showing rug plot instead.")
      df_raw <- data.frame(predictions_raw = x$predictions_raw)
      p_cal <- p_cal + 
        ggplot2::geom_rug(data = df_raw, 
                          ggplot2::aes(x = .data$predictions_raw, y = NULL),
                          alpha = 0.3, color = "#2E86AB") +
        ggplot2::labs(x = "Predicted probability")
      return(p_cal)
    }
    
    # Remove x-axis label from calibration plot
    p_cal <- p_cal + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    # Create histogram
    df_raw <- data.frame(predictions_raw = x$predictions_raw)
    p_hist <- ggplot2::ggplot(df_raw, ggplot2::aes(x = .data$predictions_raw)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2E86AB", color = "white", alpha = 0.8) +
      ggplot2::labs(x = "Predicted probability", y = "Count") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 5.5, b = 5.5, l = 5.5))
    
    # Combine plots
    p_combined <- patchwork::wrap_plots(p_cal, p_hist, ncol = 1, heights = c(3, 1))
    return(p_combined)
  }
  
  # No histogram - add x-axis label
  p_cal <- p_cal + ggplot2::labs(x = "Predicted probability")
  return(p_cal)
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
