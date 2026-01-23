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


#' Simulated Transportability Data
#'
#' A simulated dataset for demonstrating transportability analysis of
#' prediction model performance. The data includes source (RCT) and target
#' populations, where treatment is randomized in the source but confounded
#' in the target.
#'
#' @format A data frame with 1500 rows and 7 variables:
#' \describe{
#'   \item{age}{Standardized age (mean 0, SD 1)}
#'   \item{biomarker}{Continuous biomarker value}
#'   \item{smoking}{Binary smoking status (1 = smoker, 0 = non-smoker)}
#'   \item{source}{Population indicator (1 = source/RCT, 0 = target).
#'     RCT patients tend to be younger and less likely to smoke.}
#'   \item{treatment}{Binary treatment indicator (1 = treated, 0 = untreated).
#'     Randomized (50/50) in source, confounded by age and biomarker in target.}
#'   \item{event}{Binary outcome indicating event (1 = event, 0 = no event).
#'     Risk depends on age, biomarker, smoking, and is reduced by treatment.}
#'   \item{risk_score}{Predicted probability of event from a model trained
#'     on the source population, using age, biomarker, and smoking as predictors.}
#' }
#'
#' @details
#' The data generating process creates realistic heterogeneity between
#' source (RCT) and target populations:
#'
#' \itemize{
#'   \item **Selection into RCT**: RCT patients are younger on average and
#'     less likely to be smokers, reflecting typical trial enrollment patterns.
#'   \item **Treatment assignment**:
#'     \itemize{
#'       \item Source: Randomized with \code{P(A=1) = 0.5}
#'       \item Target: Confounded with \code{P(A=1|X) = plogis(-0.3 + 0.015*age + 0.2*biomarker)}
#'     }
#'   \item **Outcome model**: \code{P(Y=1|X,A) = plogis(-2.5 + 0.03*age + 0.4*biomarker + 0.3*smoking - 0.5*A)}
#' }
#'
#' This creates a scenario where naive performance estimates from the RCT
#' will not generalize to the target population due to covariate shift.
#'
#' @examples
#' data(transport_sim)
#' head(transport_sim)
#'
#' # Population sizes
#' table(transport_sim$source)  # 0=target, 1=source
#'
#' # Estimate transportable MSE under no treatment
#' result <- tr_mse(
#'   predictions = transport_sim$risk_score,
#'   outcomes = transport_sim$event,
#'   treatment = transport_sim$treatment,
#'   source = transport_sim$source,
#'   covariates = transport_sim[, c("age", "biomarker", "smoking")],
#'   treatment_level = 0,
#'   analysis = "transport",
#'   estimator = "dr"
#' )
#' result
#'
#' @source Simulated data based on the framework in Voter et al. (2025).
#'   "Transportability of machine learning-based counterfactual prediction models."
#'   Diagnostic and Prognostic Research, 9(4).
#'   \doi{10.1186/s41512-025-00201-y}
#'
"transport_sim"
