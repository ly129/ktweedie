#' Predict outcome using estimated kernel Tweedie model coefficients
#'
#' \code{TD_predict()} predicts the coefficients of the Tweedie compound poisson gamma model in the RKHS. Input is the estimated coefficients, covariate matrix of new data whose outcomes are to be predicted.
#'
#' @param x Covariate matrix.
#' @param fit Fitted model from \code{\link{TD_estimate}}
#'
#' @details
#' \code{TD_predict()} uses the fitted model from \code{\link{TD_estimate}} to estimate the mean claim cost for new data points.
#'
#' @examples
#' kern <- rbfdot(sigma = 0.01)
#' # X: N * p matrix
#' \dontrun{fit <- TD_estimate(x = X, z = Z, kern = kern,
#'                             LamReg = 10, phi = 1.8,
#'                             rho = 1.2, ftol = 1e-10,
#'                             partol = 1e-15, maxit = 1e6,
#'                             quiet = TRUE)}
#' # newdata: n * p matrix
#' \dontrun{pred <- TD_predict(x = newdata, model = fit)}
#'
#' @export
#'
#'
ktd_predict <- function(x.test, model, which.hyperparam = 1, type = "link", KernelMatrix) {
  if (!missing(x.test)) {
    x.test <- as.matrix(x.test)
  }
  sk <- model$sparsity
  fit <- model$estimates[[which.hyperparam]]
  hyperparameter <- list(lambda1 = model$data$lambda1[which.hyperparam])
  if (sk) {
    hyperparameter$lambda2 <- model$data$lambda2[which.hyperparam]
    sigma <- model$data$sigma[which.hyperparam]
    hyperparameter$sigma <- sigma
    kern <- rbfdot(sigma = sigma)
    wt <- fit$weight
    wx <- t(apply(model$data$x, MARGIN = 1, FUN = '*', wt))
    wx.test <- t(apply(x.test, MARGIN = 1, FUN = '*', wt))
    K <- as.matrix(kernelMatrix(kernel = kern, x = wx, y = wx.test))
  } else {
    kern <- model$data$kernel

    if (is.matrix(kern)) {
      if (missing(KernelMatrix)) {
        K <- kern
      } else {
        K <- KernelMatrix
      }
    } else {
      hyperparameter$sigma <- model$data$sigma[which.hyperparam]
      K <- as.matrix(kernelMatrix(kernel = kern, x = model$data$x, y = x.test))
    }
  }

  beta <- ifelse(is.null(fit$intercept), 0, fit$intercept)
  alpha <- fit$coefficient
  if (type == "link") {
    pred <- beta + crossprod(K, alpha)
  } else {
    if (type == "response") {
      pred <- exp(beta + crossprod(K, alpha))
    } else {
      pred <- NULL
    }
  }
  if (fit$convergence > 7) {
    return(NULL)
  } else {
    return(list(prediction = pred,
                fit = fit,
                hyperparameter = hyperparameter,
                kernel = kern,
                kernel.matrix = K))
  }


}
