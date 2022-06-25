#' Cross validation for tuning the regularization coefficient in the kernel Tweedie model
#'
#' \code{ktd_cv()} performs cross-validation to determine the regularization coefficient of the \code{ktweedie} model.
#'
#' @param x Covariate matrix.
#' @param y Outcome vector (e.g. insurance cost).
#' @param kern Choice of kernel. See \code{\link{dots}} for details on supported kernel functions.
#' @param lambda A vector of candidate regularization coefficients used in cross-validation.
#' @param nfolds Number of folds in cross-validation. Default is 5.
#' @param rho The power parameter of the Tweedie model. Default is 1.5 and can take any real value between 1 and 2.
#' @param loss Criterion used in cross-validation. "LL" for log likelihood, "RMSE" for root mean squared error, "MAD" for mean absolute difference. Default is "LL".
#' @param ... Optional arguments to be passed to \code{\link{ktd_estimate}()}.
#'
#' @details
#' \code{ktd_cv()} is a built-in wrapper for cross-validation for the choice of regularization coefficient.
#'
#' @seealso \code{\link{ktd_cv2d}}, \code{\link{ktd_estimate}}, \code{\link{ktd_predict}}
#' @examples
#' ( cv1d <- ktd_cv(x = dat$x, y = dat$y,
#'                  kern = rbfdot(sigma = 1e-8),
#'                  lambda = 10^(-8:-1),
#'                  nfolds = 5) )
#' @export
#'
#'
ktd_cv <- function(x, y, kern, lambda, nfolds= 5, rho = 1.5, loss = "LL", ...) {
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  N <- nrow(x)
  nz <- length(y)
  if (N != nz) stop("Sample sizes are different in x and y.")

  # sort lambda in decreasing order for warm-start
  lambda <- sort(lambda, decreasing = TRUE)
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds = 5 recommended")
  }

  loss.types <- c("LL", "RMSE", "MAD")
  if (! loss %in% loss.types) {
    stop("Invalid loss type for cross validation.")
  }

  foldid <- sample(rep(seq(nfolds), length = N))

  ###Now fit the nfold models and store them
  no.hyper <- length(lambda)
  predmat <- matrix(NA, nrow = N, ncol = no.hyper,
                    dimnames = list(NULL, lambda))
  for (i in seq(nfolds)) {
    which <- foldid == i
    x_train <- x[!which, , drop = FALSE]
    x_test <- x[which, , drop = FALSE]
    y_train <- y[!which]
    y_test <- y[which]

    model <- ktd_estimate(x = x_train, y = y_train, kern, lam1 = lambda, rho, ...)

    K <- as.matrix(kernelMatrix(kernel = kern, x = x_train, y = x_test))

    for (l in lambda) {
      fit <- model$estimates[[paste("lambda", l)]]
      alphaVec <- fit$coefficient
      K_alpha <- crossprod(K, alphaVec)
      # likelihood-based
      if (loss == "LL") {
        Kexp1 <- exp(-(rho-1) * K_alpha)
        Kexp2 <- exp( (2-rho) * K_alpha)
        beta.num <- drop(y_test * Kexp1)
        beta.den <- Kexp2

        if (is.null(fit$intercept)) {
          lkhd <- -(beta.num/(rho-1) + beta.den/(2-rho))
        } else {
          beta <- fit$intercept
          Bexp1 <- exp(-(rho-1) * beta)
          Bexp2 <- exp( (2-rho) * beta)
          lkhd <- -(beta.num * Bexp1/(rho-1) + beta.den * Bexp2/(2-rho))
        }
        predmat[which, paste(l)] <- lkhd
      }

      # rmse and mad based
      if (loss == "RMSE" | loss == "MAD") {
        beta <- ifelse(is.null(fit$intercept), 0, fit$intercept)
        pred <- exp(beta + K_alpha)
        predmat[which, paste(l)] <- pred
      }
    }
  }

  if (loss == "LL") {
    lossVec <- apply(predmat, MARGIN = 2, FUN = mean)
  }

  if (loss == "RMSE") {
    lossVec <- apply(predmat, MARGIN = 2, FUN = function(x) {sqrt(mean((x - y)^2))})
  }

  if (loss == "MAD") {
    lossVec <- apply(predmat, MARGIN = 2, FUN = function(x) {mean(abs(x - y))} )
  }

  names(lossVec) <- lambda
  if (loss == "LL") {
    lambda.min <- lambda[which.max(lossVec)]
  } else {
    lambda.min <- lambda[which.min(lossVec)]
  }
  results <- list(lossVec, Best_lambda = lambda.min)
  names(results)[1] <- loss
  return(results)
}
