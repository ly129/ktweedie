#' A demo dataset
#'
#' A simulated dataset with covariate matrix \code{x} of size 200 x 50 and an outcome vector \code{y} of length 200.
#'
#'
#' @name dat
#' @keywords datasets
#' @usage data(dat)
#'
#' @format A list with 2 items:
#' \describe{
#'   \item{x}{Covariate matrix}
#'   \item{y}{Outcome vector}
#' }
#'
#' @details \code{x} is generated from standard normal distribution. \code{y} is generated from Tweedie distribution with mean equal to exp(sin(x) %*% (6, -4, 3, 2, -2,0,...,0)). Only the first five variables are associated with the outcome.
#' @docType data
NULL
