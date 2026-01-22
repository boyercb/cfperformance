#' Simulated Cardiovascular Disease Data
#'
#' A simulated dataset for demonstrating counterfactual prediction model
#' performance estimation. The data mimics a scenario where patients receive
#' a cardiovascular treatment based on their risk factors, and we want to
#' evaluate how well a prediction model would perform under different 
#' treatment policies.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{age}{Standardized age (mean 0, SD 1)}
#'   \item{bp}{Standardized blood pressure (mean 0, SD 1)}
#'   \item{chol}{Standardized cholesterol level (mean 0, SD 1)}
#'   \item{treatment}{Binary treatment indicator (1 = treated, 0 = untreated).
#'     Treatment assignment is confounded by age and blood pressure.}
#'   \item{event}{Binary outcome indicating cardiovascular event (1 = event, 
#'     0 = no event). Risk depends on age, bp, chol, and is reduced by treatment.}
#'   \item{risk_score}{Predicted probability of event from a logistic regression
#'     model fit on observed data using age, bp, and chol as predictors.}
#' }
#'
#' @details
#' The data generating process is:
#' \itemize{
#'   \item Covariates are independent standard normal
#'   \item Treatment probability: \code{plogis(-0.3 + 0.4*age + 0.3*bp + 0.1*chol)}
#'   \item Outcome probability: \code{plogis(-2 + 0.6*age + 0.5*bp + 0.3*chol - 0.4*treatment)}
#' }
#'
#' This creates confounding because sicker patients (higher age, bp) are more
#' likely to receive treatment, but treatment also affects outcomes.
#'
#' @examples
#' data(cvd_sim)
#' head(cvd_sim)
#'
#' # Estimate counterfactual MSE under no treatment
#' result <- cf_mse(
#'   predictions = cvd_sim$risk_score,
#'   outcomes = cvd_sim$event,
#'   treatment = cvd_sim$treatment,
#'   covariates = cvd_sim[, c("age", "bp", "chol")],
#'   treatment_level = 0,
#'   estimator = "dr"
#' )
#' result
#'
#' @source Simulated data based on the framework in Boyer, Dahabreh & 
#'   Steingrimsson (2025). "Estimating and evaluating counterfactual prediction models."
#'   Statistics in Medicine, 44(23-24), e70287.
#'   \doi{10.1002/sim.70287}
#'
"cvd_sim"
