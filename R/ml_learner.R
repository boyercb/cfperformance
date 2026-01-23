#' Specify a Machine Learning Learner for Nuisance Models
#'
#' Creates a learner specification that can be passed to `propensity_model`
#' or `outcome_model` arguments in [cf_mse()], [cf_auc()], [cf_calibration()],
#' and their transportability variants. When an `ml_learner` specification is
#' provided, cross-fitting is automatically used for valid inference.
#'
#' @param method Character string specifying the learner type:
#'   - `"ranger"`: Random forest via [ranger::ranger()]
#'   - `"xgboost"`: Gradient boosting via [xgboost::xgboost()]
#'   - `"grf"`: Generalized random forest via [grf::regression_forest()] or
#'     [grf::probability_forest()]
#'   - `"glmnet"`: Regularized regression via [glmnet::cv.glmnet()]
#'   - `"superlearner"`: Ensemble via [SuperLearner::SuperLearner()]
#'   - `"custom"`: User-supplied fit/predict functions
#' @param ... Additional arguments passed to the fitting function.
#' @param fit_fun For `method = "custom"`, a function with signature
#'   `function(formula, data, family, ...)` that returns a fitted model object.
#' @param predict_fun For `method = "custom"`, a function with signature
#'   `function(object, newdata, ...)` that returns predicted probabilities.
#'
#' @return An object of class `ml_learner` containing the learner specification.
#'
#' @details
#' ## Supported Learners
#'
#' **ranger**: Fast random forest implementation. Key arguments:
#' - `num.trees`: Number of trees (default: 500)
#' - `mtry`: Number of variables to sample at each split
#' - `min.node.size`: Minimum node size
#'
#' **xgboost**: Gradient boosting. Key arguments:
#' - `nrounds`: Number of boosting rounds (default: 100)
#' - `max_depth`: Maximum tree depth (default: 6)
#' - `eta`: Learning rate (default: 0.3)
#'
#' **grf**: Generalized random forests with built-in honesty. Key arguments:
#' - `num.trees`: Number of trees (default: 2000)
#' - `honesty`: Whether to use honest estimation (default: TRUE)
#'
#' **glmnet**: Elastic net regularization with cross-validation. Key arguments:
#' - `alpha`: Elastic net mixing parameter (0 = ridge, 1 = lasso, default: 1)
#' - `nfolds`: Number of CV folds for lambda selection (default: 10)
#'
#' **superlearner**: Ensemble of multiple learners. Key arguments:
#' - `SL.library`: Vector of learner names (default: `c("SL.glm", "SL.ranger")`)
#'
#' @seealso [cf_mse()], [cf_auc()], [cf_calibration()]
#'
#' @importFrom stats model.frame model.response
#' @importFrom utils modifyList
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Random forest for propensity score
#' cf_mse(
#'   Y = outcome, A = treatment, predictions = preds, data = df,
#'   propensity_formula = treatment ~ x1 + x2,
#'   propensity_model = ml_learner("ranger", num.trees = 500),
#'   cross_fit = TRUE
#' )
#'
#' # XGBoost with custom parameters
#' cf_mse(
#'   Y = outcome, A = treatment, predictions = preds, data = df,
#'   propensity_formula = treatment ~ x1 + x2,
#'   propensity_model = ml_learner("xgboost", nrounds = 200, max_depth = 4),
#'   cross_fit = TRUE
#' )
#'
#' # Custom learner
#' my_fit <- function(formula, data, family, ...) {
#'   glm(formula, data = data, family = binomial())
#' }
#' my_predict <- function(object, newdata, ...) {
#'   predict(object, newdata = newdata, type = "response")
#' }
#' cf_mse(
#'   ...,
#'   propensity_model = ml_learner("custom", fit_fun = my_fit,
#'                                  predict_fun = my_predict)
#' )
#' }
ml_learner <- function(method = c("ranger", "xgboost", "grf", "glmnet",
                                   "superlearner", "custom"),
                       ...,
                       fit_fun = NULL,
                       predict_fun = NULL) {

  method <- match.arg(method)

  # Validate custom method has required functions

if (method == "custom") {
    if (is.null(fit_fun) || !is.function(fit_fun)) {
      stop("`fit_fun` must be provided as a function for method = 'custom'")
    }
    if (is.null(predict_fun) || !is.function(predict_fun)) {
      stop("`predict_fun` must be provided as a function for method = 'custom'")
    }
  }

  structure(
    list(
      method = method,
      args = list(...),
      fit_fun = fit_fun,
      predict_fun = predict_fun
    ),
    class = "ml_learner"
  )
}


#' Print Method for ml_learner Objects
#' @param x An `ml_learner` object.
#' @param ... Ignored.
#' @export
print.ml_learner <- function(x, ...) {
  cat("ML Learner Specification\n")
  cat("------------------------\n")
  cat("Method:", x$method, "\n")
  if (length(x$args) > 0) {
    cat("Arguments:\n")
    for (nm in names(x$args)) {
      val <- x$args[[nm]]
      if (is.character(val)) {
        cat("  ", nm, ": ", paste(val, collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", nm, ": ", val, "\n", sep = "")
      }
    }
  }
  invisible(x)
}


# =============================================================================
# Internal: Fit ML Learner
# =============================================================================

#' Fit an ML Learner to Data
#'
#' @param learner An `ml_learner` object.
#' @param formula Model formula.
#' @param data Training data.
#' @param family Character: "binomial" for classification, "gaussian" for
#'   regression.
#'
#' @return A fitted model object with class `ml_fitted`.
#' @keywords internal
.fit_ml_learner <- function(learner, formula, data, family = "binomial") {

  fit <- switch(learner$method,
    ranger = .fit_ranger(formula, data, family, learner$args),
    xgboost = .fit_xgboost(formula, data, family, learner$args),
    grf = .fit_grf(formula, data, family, learner$args),
    glmnet = .fit_glmnet(formula, data, family, learner$args),
    superlearner = .fit_superlearner(formula, data, family, learner$args),
    custom = .fit_custom(formula, data, family, learner$args,
                         learner$fit_fun, learner$predict_fun)
  )

  # Wrap in ml_fitted class for consistent prediction dispatch
  structure(
    list(
      fit = fit,
      method = learner$method,
      formula = formula,
      predict_fun = learner$predict_fun  # For custom method
    ),
    class = "ml_fitted"
  )
}


#' Predict from Fitted ML Learner
#'
#' @param object An `ml_fitted` object.
#' @param newdata Data frame for prediction.
#' @param type Character: "response" for probabilities (default), "class" for
#'   class predictions.
#' @param ... Ignored.
#'
#' @return Numeric vector of predicted probabilities.
#' @keywords internal
.predict_ml_learner <- function(object, newdata, type = "response", ...) {

  if (!inherits(object, "ml_fitted")) {
    stop("Expected an 'ml_fitted' object")
  }

  pred <- switch(object$method,
    ranger = .predict_ranger(object$fit, newdata),
    xgboost = .predict_xgboost(object$fit, newdata),
    grf = .predict_grf(object$fit, newdata),
    glmnet = .predict_glmnet(object$fit, newdata, object$formula),
    superlearner = .predict_superlearner(object$fit, newdata, object$formula),
    custom = object$predict_fun(object$fit, newdata)
  )

  as.numeric(pred)
}


# =============================================================================
# ranger Implementation
# =============================================================================

.fit_ranger <- function(formula, data, family, args) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required. Install with: install.packages('ranger')",
         call. = FALSE)
  }

  # Set defaults for classification
  default_args <- list(
    formula = formula,
    data = data,
    probability = (family == "binomial"),
    num.trees = 500
  )

  # Merge user args (user args take precedence)
  call_args <- modifyList(default_args, args)

  do.call(ranger::ranger, call_args)
}

.predict_ranger <- function(fit, newdata) {
  pred <- predict(fit, data = newdata)

  if (fit$treetype == "Probability estimation") {
    # Return probability of positive class (second column)
    if (ncol(pred$predictions) >= 2) {
      return(pred$predictions[, 2])
    }
  }
  pred$predictions
}


# =============================================================================
# xgboost Implementation
# =============================================================================

.fit_xgboost <- function(formula, data, family, args) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' required. Install with: install.packages('xgboost')",
         call. = FALSE)
  }

  # Extract response and predictors from formula
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = data)[, -1, drop = FALSE]  # Remove intercept

  # Convert to DMatrix
  dtrain <- xgboost::xgb.DMatrix(data = x, label = y)

  # Set defaults based on family
  objective <- if (family == "binomial") "binary:logistic" else "reg:squarederror"

  default_args <- list(
    data = dtrain,
    nrounds = 100,
    objective = objective,
    verbose = 0
  )

  # Merge user args
  call_args <- modifyList(default_args, args)

  # Store formula info for prediction
  fit <- do.call(xgboost::xgboost, call_args)
  attr(fit, "formula") <- formula
  attr(fit, "feature_names") <- colnames(x)

  fit
}

.predict_xgboost <- function(fit, newdata) {
  formula <- attr(fit, "formula")

  # Create model matrix from newdata
  x <- model.matrix(formula, data = newdata)[, -1, drop = FALSE]

  # Make prediction
  predict(fit, newdata = x)
}


# =============================================================================
# grf Implementation
# =============================================================================

.fit_grf <- function(formula, data, family, args) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' required. Install with: install.packages('grf')",
         call. = FALSE)
  }

  # Extract response and predictors
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = data)[, -1, drop = FALSE]

  default_args <- list(
    X = x,
    Y = y,
    num.trees = 2000
  )

  # Merge user args
  call_args <- modifyList(default_args, args)

  # Use probability_forest for binary outcomes
  if (family == "binomial") {
    # Convert to factor for probability_forest
    call_args$Y <- as.factor(y)
    fit <- do.call(grf::probability_forest, call_args)
  } else {
    fit <- do.call(grf::regression_forest, call_args)
  }

  attr(fit, "formula") <- formula
  fit
}

.predict_grf <- function(fit, newdata) {
  # GRF stores the training data columns - use those names
  # For probability_forest and regression_forest, the X matrix column names are stored
  train_cols <- colnames(fit$X.orig)

  if (is.null(train_cols)) {
    # Fallback: use all columns from newdata
    x <- as.matrix(newdata)
  } else {
    # Use the same columns as training
    x <- as.matrix(newdata[, train_cols, drop = FALSE])
  }

  pred <- predict(fit, newdata = x)

  # probability_forest returns matrix with class probabilities
  if (is.matrix(pred$predictions) && ncol(pred$predictions) >= 2) {
    return(pred$predictions[, 2])
  }
  pred$predictions
}


# =============================================================================
# glmnet Implementation
# =============================================================================

.fit_glmnet <- function(formula, data, family, args) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' required. Install with: install.packages('glmnet')",
         call. = FALSE)
  }

  # Extract response and predictors
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = data)[, -1, drop = FALSE]

  default_args <- list(
    x = x,
    y = y,
    family = family,
    alpha = 1,  # Lasso by default
    nfolds = 10
  )

  # Merge user args
  call_args <- modifyList(default_args, args)

  # Use cv.glmnet for automatic lambda selection
  fit <- do.call(glmnet::cv.glmnet, call_args)
  attr(fit, "formula") <- formula

  fit
}

.predict_glmnet <- function(fit, newdata, formula) {
  if (is.null(formula)) formula <- attr(fit, "formula")
  x <- model.matrix(formula, data = newdata)[, -1, drop = FALSE]

  # Predict at lambda.min (optimal from CV)
  pred <- predict(fit, newx = x, s = "lambda.min", type = "response")
  as.vector(pred)
}


# =============================================================================
# SuperLearner Implementation
# =============================================================================

.fit_superlearner <- function(formula, data, family, args) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' required. Install with: install.packages('SuperLearner')",
         call. = FALSE)
  }

  # Extract response and predictors
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  x <- as.data.frame(model.matrix(formula, data = data)[, -1, drop = FALSE])

  # Default library
  default_library <- c("SL.glm", "SL.ranger")

  default_args <- list(
    Y = y,
    X = x,
    family = if (family == "binomial") binomial() else gaussian(),
    SL.library = default_library
  )

  # Merge user args
  call_args <- modifyList(default_args, args)

  fit <- do.call(SuperLearner::SuperLearner, call_args)
  attr(fit, "formula") <- formula
  attr(fit, "feature_names") <- colnames(x)

  fit
}

.predict_superlearner <- function(fit, newdata, formula) {
  if (is.null(formula)) formula <- attr(fit, "formula")
  x <- as.data.frame(model.matrix(formula, data = newdata)[, -1, drop = FALSE])

  pred <- predict(fit, newdata = x)
  as.vector(pred$pred)
}


# =============================================================================
# Custom Implementation
# =============================================================================

.fit_custom <- function(formula, data, family, args, fit_fun, predict_fun) {
  # Call user-supplied fit function
  do.call(fit_fun, c(list(formula = formula, data = data, family = family), args))
}


# =============================================================================
# Check if object is an ml_learner
# =============================================================================

#' Check if Object is an ML Learner Specification
#'
#' @param x Object to check.
#' @return Logical indicating if `x` is an `ml_learner` object.
#' @keywords internal
is_ml_learner <- function(x) {
  inherits(x, "ml_learner")
}


#' Check if Object is a Fitted ML Model
#'
#' @param x Object to check.
#' @return Logical indicating if `x` is an `ml_fitted` object.
#' @keywords internal
is_ml_fitted <- function(x) {
  inherits(x, "ml_fitted")
}


#' Unified Prediction from Nuisance Models
#'
#' Handles prediction from both standard glm/gam models and ml_fitted objects.
#'
#' @param model A fitted model (glm, gam, or ml_fitted).
#' @param newdata Data frame for prediction.
#' @param type Prediction type (default: "response").
#'
#' @return Numeric vector of predictions.
#' @keywords internal
.predict_nuisance <- function(model, newdata, type = "response") {
  if (is_ml_fitted(model)) {
    .predict_ml_learner(model, newdata, type = type)
  } else {
    # Standard model (glm, gam, etc.)
    predict(model, newdata = newdata, type = type)
  }
}
