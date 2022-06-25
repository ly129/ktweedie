#' ktweedie Package
#'
#' Fast estimation and prediction of insurance cost using the Tweedie compound poisson gamma model. The package supports choice of different kernels and built-in cross-validation tools.
#'
#' The ktweedie package includes functions ktd_estimate, ktd_cv, ktd_cv2d and ktd_predict.
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
#' @author Yi Lian \email{yi.lian@mail.mcgill.ca}
#'
#' @name ktweedie
NULL
