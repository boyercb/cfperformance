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


#' Trim propensity scores to avoid extreme values
#'
#' Internal function for trimming propensity scores using either absolute bounds
#' or quantile-based trimming.
#'
#' @param ps Numeric vector of propensity scores (probabilities).
#' @param method Character string specifying trimming method:
#'   - `"absolute"`: Trim to fixed probability bounds (default)
#'   - `"quantile"`: Trim based on quantiles of the propensity distribution
#'   - `"none"`: No trimming
#' @param bounds Numeric vector of length 2 specifying trimming bounds.
#'   For `method = "absolute"`: probability bounds (default: c(0.01, 0.99)).
#'   For `method = "quantile"`: quantile bounds (default: c(0.01, 0.99)).
#' @param weights Optional numeric vector of weights for quantile calculation.
#'   Only used when `method = "quantile"`.
#'
#' @return Trimmed propensity scores.
#'
#' @details
#' For `method = "absolute"`, propensity scores are clipped to `[bounds[1], bounds[2]]`.
#' This is the most common approach and ensures that no propensity score falls below

#' or above the specified thresholds.
#'
#' For `method = "quantile"`, the bounds are computed as the `bounds[1]` and `bounds[2]`
#' quantiles of the propensity score distribution, and then scores are clipped to
#' these data-dependent bounds. This can be useful when the propensity score
#' distribution varies across datasets.
#'
#' @keywords internal
.trim_propensity <- function(ps, method = "absolute", bounds = c(0.01, 0.99),
                              weights = NULL) {
  # Handle NULL or "none" method

if (is.null(method) || identical(method, "none")) {
    return(ps)
  }

  method <- match.arg(method, c("absolute", "quantile", "none"))

  if (method == "none") {
    return(ps)
  }

  # Validate bounds
  if (length(bounds) != 2 || !is.numeric(bounds)) {
    stop("bounds must be a numeric vector of length 2", call. = FALSE)
  }
  if (bounds[1] >= bounds[2]) {
    stop("bounds[1] must be less than bounds[2]", call. = FALSE)
  }

  if (method == "absolute") {
    # Clip to fixed probability bounds
    if (bounds[1] < 0 || bounds[2] > 1) {
      stop("For method='absolute', bounds must be in [0, 1]", call. = FALSE)
    }
    ps_trimmed <- pmax(pmin(ps, bounds[2]), bounds[1])

  } else if (method == "quantile") {
    # Clip to data-dependent quantile bounds
    if (bounds[1] < 0 || bounds[2] > 1) {
      stop("For method='quantile', bounds must be quantile levels in [0, 1]", call. = FALSE)
    }
    if (!is.null(weights)) {
      # Weighted quantiles (approximation using sorted values)
      ord <- order(ps)
      ps_sorted <- ps[ord]
      w_sorted <- weights[ord]
      cum_w <- cumsum(w_sorted) / sum(w_sorted)
      lower_q <- ps_sorted[which.min(abs(cum_w - bounds[1]))]
      upper_q <- ps_sorted[which.min(abs(cum_w - bounds[2]))]
    } else {
      lower_q <- quantile(ps, bounds[1], na.rm = TRUE)
      upper_q <- quantile(ps, bounds[2], na.rm = TRUE)
    }
    ps_trimmed <- pmax(pmin(ps, upper_q), lower_q)
  }

  return(ps_trimmed)
}


#' Parse ps_trim argument into method and bounds
#'
#' Internal function to parse the `ps_trim` argument into method and bounds.
#'
#' @param ps_trim Either:
#'   - A list with `method` and `bounds` elements
#'   - A numeric vector of length 2 (interpreted as absolute bounds)
#'   - A single number (interpreted as symmetric absolute bounds: c(x, 1-x))
#'   - "none" or NULL for no trimming
#'   - "quantile" for quantile-based trimming with default bounds
#'
#' @return A list with `method` and `bounds` elements.
#'
#' @keywords internal
.parse_ps_trim <- function(ps_trim) {
  # Default: absolute trimming at c(0.01, 0.99)
  default <- list(method = "absolute", bounds = c(0.01, 0.99))

  if (is.null(ps_trim)) {
    return(default)
  }

  # Handle "none" string
  if (is.character(ps_trim)) {
    if (ps_trim == "none") {
      return(list(method = "none", bounds = NULL))
    } else if (ps_trim == "quantile") {
      return(list(method = "quantile", bounds = c(0.01, 0.99)))
    } else if (ps_trim == "absolute") {
      return(default)
    } else {
      stop("If ps_trim is a string, it must be 'none', 'absolute', or 'quantile'",
           call. = FALSE)
    }
  }

  # Handle list input
  if (is.list(ps_trim)) {
    method <- ps_trim$method %||% "absolute"
    bounds <- ps_trim$bounds %||% c(0.01, 0.99)
    return(list(method = method, bounds = bounds))
  }

  # Handle numeric input
  if (is.numeric(ps_trim)) {
    if (length(ps_trim) == 2) {
      # Two values: treat as absolute bounds
      return(list(method = "absolute", bounds = ps_trim))
    } else if (length(ps_trim) == 1) {
      # Single value: symmetric bounds c(x, 1-x)
      return(list(method = "absolute", bounds = c(ps_trim, 1 - ps_trim)))
    } else {
      stop("If ps_trim is numeric, it must have length 1 or 2", call. = FALSE)
    }
  }

  stop("ps_trim must be NULL, 'none', 'quantile', a list, or a numeric vector",
       call. = FALSE)
}


# Truncate extreme weights (legacy function - kept for compatibility)
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
