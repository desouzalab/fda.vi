#' GCV Score for a VEM Smooth Fit
#'
#' @description
#' Computes the Generalized Cross-Validation (GCV) score for each curve from a
#' \code{vem_smooth} model object. GCV approximates leave-one-out prediction
#' error and is used by \code{\link{tune_vem_by_gcv}} to select the optimal
#' number of basis functions \eqn{K}.
#'
#' The smoother matrix \eqn{S_i} maps observed values to fitted values and is
#' constructed from the variational posteriors. Its trace provides the effective
#' degrees of freedom used in the GCV penalty.
#'
#' @param out A fitted object returned by \code{\link{vem_smooth}}.
#' @param threshold Numeric in \eqn{(0,1)}. Posterior inclusion probability
#'   threshold for treating a basis as active. Default \code{0.5}.
#'
#' @return A named numeric vector of length \eqn{m} of per-curve GCV scores.
#'   Lower scores indicate better fit relative to model complexity.
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{tune_vem_by_gcv}}, \code{\link{vem_smooth}}
#' @export

gcv_vem <- function(out, threshold = 0.5) {

  # validating inputs
  if (is.null(out$data) || is.null(out$B)) stop("`out` must include `data` and `B`.")
  if (is.null(out$Xt)) stop("`out` must include `Xt` (update vem_smooth() to return Xt).")
  if (is.null(out$prob)) stop("`out` must include `prob`.")
  if (is.null(out$m) || is.null(out$K)) stop("`out` must include `m` and `K` (update vem_smooth()).")
  if (is.null(out$w)) stop("`out` must include `w`.")

  y <- out$data
  B <- out$B
  Xt <- out$Xt
  m <- out$m
  K <- out$K

  #sanity check for prob vector dimensuions
  if (length(out$prob) != m * K) {
    stop(sprintf("Expected length(out$prob)= m*K = %d but got %d", m*K, length(out$prob)))
  }

  # reconstructing variational quantities
  # expectations
  E_inv_sigma2 <- expectedInvSigma2(out$delta1, out$delta2)
  E_inv_tau2   <- expectedInvTau2(out$lambda1, out$lambda2)

  # psi (correlation) matrices
  Psi <- computePsiMatrix(Xt, Xt, w = out$w)
  Psi_inv <- solve(Psi)

  seq_values <- lapply(seq(1, m*K, by = K), function(s) s:(s + K - 1))
  gcv_vals <- numeric(m)

  # gcv for each curve
  for (i in 1:m) {
    idx <- seq_values[[i]]

    p_i <- out$prob[idx] #inclusion probs
    Zi_hat <- ifelse(p_i > threshold, 1, 0) #threshold selection
    mu_i <- out$mu_beta[idx] #posterior mean coefficients

    # like design matrix
    EG_i <- p_i * t(B[[i]])

    # computing smoother matrix (maps y obs to y_hat values)
    S_i <- E_inv_sigma2 *
      B[[i]] %*% diag(Zi_hat) %*%
      solve(E_inv_sigma2 * (diag(E_inv_tau2, K) + EG_i %*% Psi_inv %*% t(EG_i))) %*%
      EG_i %*% Psi_inv

    # computing residual sum of squares (RSS)
    yhat <- as.numeric(B[[i]] %*% (Zi_hat * mu_i))
    resid <- y[[i]] - yhat
    RSS <- sum(resid^2)
    n_i <- length(y[[i]])
    df <- sum(diag(S_i))

    #computing GCV
    gcv_vals[i] <- (RSS / n_i) / (1 - df / n_i)^2
  }

  names(gcv_vals) <- paste0("curve_", seq_len(m))
  gcv_vals
}
