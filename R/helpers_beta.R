#' Expected Residual and Coefficient Squared Moments
#'
#' @description
#' Three expectations under \eqn{q(\beta_i)} and \eqn{q(Z_i)} appearing in the
#' updates for \eqn{q(\sigma^2)}, \eqn{q(\tau^2)}, and \eqn{q(Z_{ki})}.
#'
#' \code{expectedResidualSq} computes the expected weighted residual sum of
#' squares \eqn{E[(y_i - B_i\xi_i)^\top \Psi^{-1}(y_i - B_i\xi_i)]} via a
#' quadratic-form plus trace decomposition. Uses \code{prob[iter-1,]} because
#' \eqn{q(Z)} has not yet been updated at the point this is called (Step 2).
#'
#' \code{expectedSumBetaSq} computes
#' \eqn{\sum_k(\text{Var}_{\beta_{ki}} + \mu_{\beta_{ki}}^2)}, appearing in
#' the \eqn{q(\sigma^2)} and \eqn{q(\tau^2)} updates.
#'
#' \code{expectedBetaSq} computes the same residual quadratic form as
#' \code{expectedResidualSq} but with the \eqn{k}-th inclusion indicator fixed
#' at candidate value \eqn{z \in \{0,1\}}, used in the \eqn{q(Z_{ki})} update.
#'
#' @param B List of \eqn{n_i \times K} basis matrices.
#' @param i Integer. Curve index.
#' @param y List of observed curve vectors.
#' @param mu Matrix (\eqn{\text{maxIter} \times mK}) of posterior beta means.
#' @param Sigma Array (\eqn{K \times K \times m}) of posterior beta covariances.
#' @param prob Matrix (\eqn{\text{maxIter} \times mK}) of inclusion probabilities.
#' @param iter Integer. Current iteration index.
#' @param psi Correlation matrix \eqn{\Psi} from \code{computePsiMatrix}.
#' @param z Integer (0 or 1). Candidate inclusion value (\code{expectedBetaSq} only).
#' @param k Integer. Basis function index (\code{expectedBetaSq} only).
#' @param K Integer. Total number of basis functions (\code{expectedBetaSq} only).
#'
#' @return A numeric scalar.
#' @keywords internal
#' @name beta_expectations
NULL


expectedResidualSq <- function(B, i, y, mu, Sigma, prob, iter, psi) {

  #extracts beta-index vector for curve i
  idx_i <- grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))

  #pips
  prob_i <- prob[iter - 1, idx_i]
  mu_bi <- mu[iter, idx_i]
  Sigmai <- Sigma[,,i]

  #Var(zi) - diag mat of Bernoulli variances
  var_z <- diag(prob_i*(1-prob_i))

  Psi_inv <- solve(psi)

  Var_zi_betai <- var_z*Sigmai + var_z*(tcrossprod(mu_bi)) + Sigmai*(tcrossprod(prob_i))

  #quadratic form in residual mean
  resid <- y[[i]]-B[[i]]%*%(mu_bi * prob_i)
  quad_part <- t(resid)%*%Psi_inv%*%resid


  #trace term for residual variance
  trace_part <- sum(diag(Psi_inv %*% (B[[i]] %*% Var_zi_betai %*% t(B[[i]]))))

  as.numeric(quad_part+trace_part)
}


expectedSumBetaSq <- function(i, mu, Sigma, iter) {
  idx_i <- grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))
  mu_bi <- mu[iter, idx_i]
  var_bi <- diag(Sigma[,,i])
  sum(var_bi + mu_bi^2)
}


expectedBetaSq <- function(z, i, prob, mu, Sigma, B, y, k, K, iter, psi) {
  idx_i <- grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))
  mu_bi <- mu[iter, idx_i]
  Sigmai <- Sigma[,,i]
  Psi_inv <- solve(psi)

  # vector of probabilities w/ k-th element/position is replaced with z
  prob_i_notk <- prob
  prob_i_notk[k] <- z

  var_zi_notk <- prob_i_notk * (1 - prob_i_notk)
  var_zi_notk[k] <- 0
  Var_zi_notk <- diag(var_zi_notk)

  Var_zi_betai_notk <- Var_zi_notk * Sigmai +
    Var_zi_notk * (tcrossprod(mu_bi)) +
    Sigmai * (tcrossprod(prob_i_notk))

  # mean residual given prob_i_notk
  fitted <- t(prob_i_notk*t(B[[i]])) %*% mu_bi
  resid <- y[[i]] - fitted

  quad_part <- t(resid)%*%Psi_inv%*%resid
  trace_part <- sum(diag(Psi_inv%*% (B[[i]] %*% Var_zi_betai_notk %*% t(B[[i]]))))

  as.numeric(quad_part+trace_part)
}
