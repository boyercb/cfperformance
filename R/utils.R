# Internal utility functions for cfperformance

# Check if a package is available
.check_package <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(purpose)) {
      msg <- paste0(msg, " for ", purpose)
    }
    msg <- paste0(msg, ". Please install it.")
    stop(msg, call. = FALSE)
  }
}


# Truncate extreme weights
.truncate_weights <- function(weights, lower = 0.01, upper = 0.99) {
  lower_q <- quantile(weights, lower, na.rm = TRUE)
  upper_q <- quantile(weights, upper, na.rm = TRUE)
  pmin(pmax(weights, lower_q), upper_q)
}


# Stabilize weights
.stabilize_weights <- function(weights, treatment, treatment_level) {
  # Multiply by marginal probability
  p_marginal <- mean(treatment == treatment_level)
  weights * p_marginal
}


# Safe division (avoid Inf)
.safe_divide <- function(num, denom, default = 0) {
  result <- num / denom
  result[!is.finite(result)] <- default
  result
}


# Format confidence interval for printing
.format_ci <- function(lower, upper, digits = 3) {
  sprintf("[%.*f, %.*f]", digits, lower, digits, upper)
}


# Create progress bar (if progressr available)
.make_progress <- function(n, ...) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps = n, ...)
  } else {
    function(...) invisible(NULL)
  }
}
