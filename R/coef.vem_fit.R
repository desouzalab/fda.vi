#' Extract Active Basis Coefficients from a VEM Fit
#'
#' @description
#' Returns a \eqn{K \times m} matrix of estimated basis coefficients.
#' Each column corresponds to one curve; each row to one basis function.
#' Coefficients are set to zero when the posterior inclusion probability
#' \eqn{p_{ki} \leq} \code{threshold} (inactive bases). When
#' \code{is.composite = TRUE}, the matrix has dimension \eqn{\max(K) \times m},
#' where \eqn{\max(K)} is the highest \eqn{K} selected by GCV across all
#' curves; coefficients for curves with smaller optimal \eqn{K} are
#' zero-padded (structural padding).
#'
#' @param object A \code{vem_fit} object from \code{\link{vem_fit}}.
#' @param threshold Numeric in \eqn{(0,1)}. Posterior inclusion probability
#'   below which a coefficient is set to zero. Default \code{0.5}.
#' @param ... Currently unused.
#'
#' @return A numeric matrix of dimension \eqn{\max(K) \times m}, with row
#'   names \code{B1, B2, \ldots} and column names \code{Curve_1, Curve_2, \ldots}.
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{predict.vem_fit}},
#'   \code{\link{summary.vem_fit}}
#' @export
#' @examples
#' \donttest{
#'   data(toy_curves)
#'   fit <- vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8)
#'
#'   # K x m matrix of active coefficients
#'   coefs <- coef(fit)
#'   dim(coefs)   # 8 x 3
#'
#'   # Compare estimated vs true coefficients for curve 1
#'   cbind(estimated = coefs[, 1], true = toy_curves$true_coef)
#'
#'   # Stricter threshold — only very confident inclusions
#'   coef(fit, threshold = 0.9)
#' }


coef.vem_fit <- function(object, threshold = 0.5, ...) {
  m <- length(object$data_orig)

  if (object$is_composite) {
    all_Ks <- vapply(object$model, function(x) as.integer(x$K), integer(1))
    max_K <- max(all_Ks)
  } else {
    max_K <- object$best_K
  }

  # rows = basis functions, cols = curves; zeros cover both inactive and structural padding
  mat <- matrix(0, nrow = max_K, ncol = m)
  colnames(mat) <- paste0("Curve_", 1:m)
  rownames(mat) <- paste0("B", 1:max_K)

  if (object$is_composite) {
    for (i in 1:m) {
      mod <- object$model[[i]]
      K_i <- mod$K
      active_vals <- mod$mu_beta * ifelse(mod$prob > threshold, 1, 0)
      mat[1:K_i, i] <- active_vals
    }
  } else {
    mu_vec   <- object$model$mu_beta
    prob_vec <- object$model$prob
    mat[] <- mu_vec * ifelse(prob_vec > threshold, 1, 0)
  }

  return(mat)
}
