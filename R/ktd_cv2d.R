#' Cross validation for jointly tuning the regularization coefficient and kernel parameter in the Kernel Tweedie Model
#'
#' \code{ktd_cv2d()} performs 2-dimensional random search from a user-specified range to determine the optimal pair of regularization coefficient and kernel parameter of the \code{ktweedie} model.
#'
#' @param x Covariate matrix.
#' @param y Outcome vector (e.g. insurance cost).
#' @param kernfunc Choice of kernel function. See \code{\link{dots}} for details on supported kernel functions.
#' @param lambda A vector of length two indicating the lower and upper bound from which candidate regularization coefficient values are sampled uniformly on the log scale.
#' @param sigma A vector of length two indicating the lower and upper bound from which candidate kernel parameter values are sampled uniformly on the log scale.
#' @param ncoefs The number of candidate \code{lambda} and \code{sigma} pairs to be evaluated.
#' @param nfolds Number of folds in cross-validation. Default is 5.
#' @param rho The power parameter of the Tweedie model. Default is 1.5 and can take any real value between 1 and 2.
#' @param loss Criterion used in cross-validation. "LL" for log likelihood, "RMSE" for root mean squared error, "MAD" for mean absolute difference. Default is "LL".
#' @param ... Optional arguments to be passed to \code{ktd_estimate()}.
#'
#' @details
#' \code{ktd_cv2d()} is a built-in wrapper for 2D random search for the regularization coefficient and kernel parameter. For kernel functions with greater than one parameters, \code{ktd_cv2d()} supports the tuning of the first one.
#'
#' @return A list of three items.
#' 1. LL or RMSE or MAD: a vector of validation error based on the user-specified \code{loss}, named by the corresponding \code{lambda} and \code{sigma} values;
#' 2. Best_lambda: the \code{lambda} value in the pair that generates the best loss;
#' 3. Best_sigma: the \code{sigma} value in the pair that generates the best loss.
#' @seealso \code{\link{ktd_cv}}, \code{\link{ktd_estimate}}, \code{\link{ktd_predict}}
#' @examples
#' ### Cross-validation
#' ( cv2d <- ktd_cv2d(x = dat$x, y = dat$y,
#'                    kernfunc = rbfdot,
#'                    lambda = c(1e-10, 1e0),
#'                    sigma = c(1e-10, 1e0),
#'                    ncoefs = 6) )
#' ### Followed by fitting
#' fit <- ktd_estimate(x = dat$x, y = dat$y,
#'                     kern = rbfdot(sigma = cv2d$Best_sigma),
#'                     lam1 = cv2d$Best_lambda)
#' @export
#'
#'
ktd_cv2d <- function(x, y, kernfunc, lambda, sigma, ncoefs, nfolds = 5, rho = 1.5, loss = "LL", ...) {
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  N <- nrow(x)
  nz <- length(y)
  if (N != nz) stop("Sample sizes are different in x and y.")

  loss.types <- c("LL", "RMSE", "MAD")
  if (! loss %in% loss.types) {
    stop("Invalid loss type for cross validation.")
  }

  foldid <- sample(rep(seq(nfolds), length = N))
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds = 5 recommended")
  }
  ###Now fit the nfold models and store them
  lam <- exp(runif(ncoefs, min = log(min(lambda)), max = log(max(lambda))))
  sig <- exp(runif(ncoefs, min = log(min(sigma)), max = log(max(sigma))))

  lam <- signif(lam, 6)
  sig <- signif(sig, 6)

  name2d <- paste("Lambda=", lam, ", Sigma=", sig, sep = "")

  predmat <- matrix(NA, nrow = N, ncol = ncoefs,
                    dimnames = list(NULL, name2d))

  for (i in seq(nfolds)) {
    which <- foldid == i
    x_train <- x[!which, , drop = FALSE]
    x_test <- x[which, , drop = FALSE]
    y_train <- y[!which]
    y_test <- y[which] #

    for (j in 1:ncoefs) {
      sigma.now <- sig[j]
      lambda.now <- lambda[j]
      kern <- sapply(sigma.now, FUN = kernfunc)[[1]]
      model <- ktd_estimate(x = x_train, y = y_train, kern, lam1 = lambda.now, rho, ...)

      K <- as.matrix(kernelMatrix(kernel = kern, x = x_train, y = x_test))
      fit <- model$estimates[[1]]
      alphaVec <- fit$coefficient
      K_alpha <- crossprod(K, alphaVec)

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
        predmat[which, j] <- lkhd
      }

      # rmse and mad based
      if (loss == "RMSE" | loss == "MAD") {
        beta <- ifelse(is.null(fit$intercept), 0, fit$intercept)
        pred <- exp(beta + K_alpha)
        predmat[which, j] <- pred
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

  names(lossVec) <- name2d
  if (loss == "LL") {
    max.index <- which.max(lossVec)
  } else {
    max.index <- which.min(lossVec)
  }
  lambda.max <- lam[max.index]
  sigma.max <- sig[max.index]
  results <- list(lossVec,
                  Best_lambda = lambda.max,
                  Best_sigma = sigma.max)
  names(results)[1] <- loss
  return(results)
}
