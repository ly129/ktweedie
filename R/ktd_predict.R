#' Predict outcome using fitted kernel Tweedie model
#'
#' \code{ktd_predict()} predicts the outcome with fitted \code{ktweedie} or \code{sktweedie} model at the user supplied new data.
#'
#' @param model Fitted model from \code{\link{ktd_estimate}}
#' @param newdata New x matrix for the prediction. If not provided, it will be the x matrix used to fit \code{model}.
#' @param which.lam1 The index of the \code{lam1} in \code{model} used in the prediction. Default is 1.
#' @param type The type of prediction to be made - "\code{link}" for the linear predictor and "\code{response}" for the predicted outcome. Default is "\code{link}".
#' @details
#' \code{ktd_predict()} uses the fitted model from \code{\link{ktd_estimate}} to estimate the mean outcome for new data points.
#'
#' @returns A list named \code{prediction} containing the vector of predicted outcomes.
#' @seealso \code{\link{ktd_estimate}}, \code{\link{ktd_cv}}, \code{\link{ktd_cv2d}}
#'
#' @examples
#' # Fit a ktweedie model
#' fit <- ktd_estimate(x = dat$x, y = dat$y,
#'                     kern = rbfdot(sigma = 1e-6),
#'                     lam1 = 10^(-5:1))
#' # Generate newx at which predictions are to be made.
#' # The newdata should have the same dimension as the original trainig data.
#' newx <- matrix(rnorm(10 * ncol(dat$x)), nrow = 10)
#' pred <- ktd_predict(model = fit, newdata = newx,
#'                     which.lam1 = 3, type = "link")
#' @export
#'
#'
ktd_predict <- function(model, newdata, which.lam1 = 1, type = "link") {
  if (!missing(newdata)) {
    newdata <- as.matrix(newdata)
  } else {
    newdata <- model$data$x
  }
  sk <- model$sparsity
  fit <- model$estimates[[which.lam1]]
  hyperparameter <- list(lambda1 = model$data$lambda1[which.lam1])

  if (sk) {
    hyperparameter$lambda2 <- model$data$lambda2[which.lam1]
    sigma <- model$data$sigma[which.lam1]
    hyperparameter$sigma <- sigma
    kern <- rbfdot(sigma = sigma)
    wt <- fit$weight
    wx <- t(apply(model$data$x, MARGIN = 1, FUN = '*', wt))
    wnewdata <- t(apply(newdata, MARGIN = 1, FUN = '*', wt))
    K <- as.matrix(kernelMatrix(kernel = kern, x = wx, y = wnewdata))
  } else {
    kern <- model$data$kernel

    K <- as.matrix(kernelMatrix(kernel = kern, x = model$data$x, y = newdata))
  }

  # beta <- ifelse(is.null(fit$intercept), 0, fit$intercept)
  alpha <- fit$coefficient
  if (type == "link") {
    pred <- crossprod(K, alpha)
  } else {
    if (type == "response") {
      pred <- exp(crossprod(K, alpha))
    } else {
      pred <- NULL
    }
  }
  if (fit$convergence > 7) {
    return(NULL)
  } else {
    return(list(prediction = pred))
                # ,
                # fit = fit,
                # hyperparameter = hyperparameter,
                # kernel = kern,
                # kernel.matrix = K))
  }


}
