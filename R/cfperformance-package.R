#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats predict coef confint glm binomial gaussian weighted.mean
#' @importFrom stats model.matrix reformulate quantile sd var qnorm approx
#' @importFrom stats loess median formula nobs plogis
#' @importFrom graphics abline hist layout par polygon lines legend
#' @importFrom grDevices rgb
#' @importFrom rlang .data abort warn inform
## usethis namespace: end
NULL

#' cfperformance: Counterfactual Prediction Model Performance
#'
#' Provides methods for estimating model performance measures (MSE, AUC,
#' calibration) under hypothetical/counterfactual interventions.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{cf_mse}}: Estimate counterfactual mean squared error
#'   \item \code{\link{cf_auc}}: Estimate counterfactual AUC
#'   \item \code{\link{cf_calibration}}: Estimate counterfactual calibration curve
#'   \item \code{\link{fit_nuisance}}: Fit nuisance models for estimation
#' }
#'
#' @section Estimators:
#' Three main estimators are implemented:
#' \itemize{
#'   \item \strong{Conditional Loss (CL)}: Relies on correct specification of
#'     the outcome model \eqn{h_a(X) = E[L(Y, \mu(X^*)) | X, A=a]}
#'   \item \strong{Inverse Probability Weighting (IPW)}: Relies on correct 
#'     specification of the propensity score \eqn{e_a(X) = Pr[A=a | X]}
#'   \item \strong{Doubly Robust (DR)}: Consistent if either the outcome model
#'     or propensity score is correctly specified
#' }
#'
#' @references
#' Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025). 
#' "Estimating and evaluating counterfactual prediction models." 
#' Statistics in Medicine. \doi{10.1002/sim.70287}
#'
#' @name cfperformance-package
NULL
