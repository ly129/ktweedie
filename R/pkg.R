#' ktweedie Package
#'
#' Kernel-based Tweedie compound Poisson gamma model using high-dimensional covariates for the analyses of zero-inflated response variables. The package features built-in estimation, prediction and cross-validation tools and supports choice of different kernel functions.
#'
#' The ktweedie package includes functions \code{ktd_estimate}, \code{ktd_cv}, \code{ktd_cv2d} and \code{ktd_predict}.
#'
#' @section Functions:
#' \code{\link{ktd_estimate}} Estimation of Tweedie model coefficients.
#'
#' \code{\link{ktd_cv}} Wrapper for cross-validation for regularization parameter tuning.
#'
#' \code{\link{ktd_cv2d}} Wrapper for cross-validation for simultaneous regularization and kernel parameter tuning.
#'
#' \code{\link{ktd_predict}} Prediction of Tweedie model outcomes.
#'
#' @docType package
#'
#' @author Yi Lian \email{yi.lian@mail.mcgill.ca}, Yi Yang, Boxiang Wang, Peng Shi, Robert W. Platt.
#'
#' @name ktweedie
NULL
