#' Ornstein-Uhlenbeck Correlation Matrix and Derivatives
#'
#' @description
#' \code{computePsiMatrix} builds the \eqn{n \times n} correlation matrix using
#' the range-normalized OU kernel \eqn{\Psi_{jl} = \exp(-w|t_j-t_l|/R)}, where
#' \eqn{R = \max(X_t) - \min(X_t)}. This normalisation makes \eqn{w}
#' scale-invariant; note the estimated \eqn{\hat{w}} is on the normalized scale
#' and converts to the paper's scale via \eqn{w_{\text{paper}} = \hat{w} / R}.
#'
#' \code{computeCovMatrix} returns \eqn{\sigma^2 \Psi}.
#'
#' \code{devPsi} returns the first and second derivatives of \eqn{\Psi} with
#' respect to \eqn{w}, used by \code{dev_elbo} during the M-step.
#'
#' @param Xt,Xp Numeric vectors of time points (rows and columns of \eqn{\Psi}).
#'   In practice \code{Xp = Xt}.
#' @param w Positive scalar. Correlation decay parameter.
#' @param sigma Positive scalar. Standard deviation for covariance scaling
#'   (\code{computeCovMatrix} only).
#'
#' @return \code{computePsiMatrix}: an \eqn{n \times n} matrix with values in
#'   \eqn{(0,1]}. \code{computeCovMatrix}: an \eqn{n \times n} covariance matrix.
#'   \code{devPsi}: a list with matrices \code{dPsi} and \code{d2Psi}.
#' @keywords internal
#' @name psi_helpers
NULL

computePsiMatrix <-  function(Xt, Xp, w=1) {
  range_scale = max(Xt) - min(Xt)
  dist_mat <- abs(outer(Xt, Xp, "-"))
  Psi <- exp(-w*dist_mat/ range_scale)
  return (Psi)
}

computeCovMatrix <- function(Xt, Xp, sigma=1, w=1) {
  sigma^2*computePsiMatrix(Xt, Xp, w)
}


devPsi <- function(Xt, Xp, w = 1) {
  range_scale <- max(Xt) - min(Xt)
  dist_mat <- abs(outer(Xt, Xp, "-"))

  kernel_base <- exp(-w*dist_mat/range_scale)

  # First derivative ∂Ψ/∂w
  dPsi <- (-dist_mat/range_scale) * kernel_base

  # Second derivative ∂²Ψ/∂w²
  d2Psi <- (dist_mat/range_scale)^2 * kernel_base

  list(dPsi = dPsi, d2Psi = d2Psi)
}

