#' ELBO Component Functions
#'
#' @description
#' Internal functions computing the six additive terms of the Evidence Lower
#' Bound (ELBO), which serves as both the convergence criterion and the M-step
#' objective in the VEM algorithm (da Cruz et al., 2024, eq. 36--42).
#'
#' The per-curve terms (\code{expectedLogLikelihood}, \code{elboInclusionTerm},
#' \code{elboBetaTerm}, \code{elboThetaTerm}) are summed over \eqn{i=1,\ldots,m}.
#' The global terms (\code{elboSigmaPriorTerm}, \code{elboTauTerm}) contribute once.
#'
#' \code{elbo_corr} assembles all six terms into the total ELBO at fixed \eqn{w}.
#'
#' \code{elbo_omega} wraps \code{elbo_corr} as a function of \eqn{w} alone,
#' passed to \code{optim()} as the M-step objective.
#'
#' \code{dev_elbo} provides the analytic gradient \eqn{\partial\text{ELBO}/\partial w},
#' passed to \code{optim()} to enable L-BFGS-B optimisation.
#'
#' @param y List of observed curve vectors.
#' @param ni Integer vector of per-curve sample sizes.
#' @param B List of basis matrices.
#' @param i Integer. Curve index (per-curve terms only).
#' @param m Integer. Number of curves.
#' @param K Integer. Number of basis functions.
#' @param iter Integer. Current iteration index.
#' @param delta_1,delta_2 Prior hyperparameters for \eqn{\sigma^2}.
#' @param lambda_1,lambda_2 Prior hyperparameters for \eqn{\tau^2}.
#' @param delta1_q,delta2_values Shape and scale trajectory of \eqn{q(\sigma^2)}.
#' @param lambda1_q,lambda2_values Shape and scale trajectory of \eqn{q(\tau^2)}.
#' @param mu_beta_values Matrix of posterior beta mean trajectories.
#' @param Sigma_beta Array of posterior beta covariance matrices.
#' @param a1_values,a2_values Shape parameter trajectories for \eqn{q(\theta_{ki})}.
#' @param prob_values Matrix of inclusion probability trajectories.
#' @param mu_ki Scalar prior hyperparameter for \eqn{\theta_{ki}}.
#' @param psi Current correlation matrix \eqn{\Psi(w)}.
#' @param w Positive scalar. Decay parameter (\code{elbo_omega}, \code{dev_elbo} only).
#' @param Xt Numeric vector of evaluation points (\code{elbo_omega}, \code{dev_elbo} only).
#' @param logdet_psi Pre-computed \eqn{\log|\Psi|} or \code{NULL}
#'   (\code{expectedLogLikelihood} only).
#'
#' @return A numeric scalar.
#' @keywords internal
#' @name elbo_functions
NULL

# ELBO 1st Term: expected log-likelihood
expectedLogLikelihood <- function(y, ni, B, i, iter,
                                  delta1_q, delta2_values,
                                  mu_beta_values, Sigma_beta,
                                  prob_values, psi,
                                  logdet_psi=NULL) {

  #sigma2: variance scale parameter at iteration iter
  delta2_q <- delta2_values[iter]

  # log determinant computed once for numerical stability
  if (is.null(logdet_psi)) {
    logdet_psi <- as.numeric(determinant(psi, logarithm = TRUE)$modulus)
  }

  #expected squared residual term
  e_res <- expectedResidualSq(
    B   = B,
    i   = i,
    y   = y,
    mu  = mu_beta_values,
    Sigma = Sigma_beta,
    prob  = prob_values,
    iter = iter,
    psi  = psi
  )

    #ELBO contribution for curve i
  res <- (-ni[i] / 2) * (log(2 * pi) + expectedLogSigma2(delta1_q, delta2_q)) -
    0.5 * logdet_psi -
    0.5 * expectedInvSigma2(delta1_q, delta2_q) * e_res

  as.numeric(res)
}

#inclusion probability term (originally diff_z_i)
elboInclusionTerm <- function(i, iter, K, prob_values, mu_beta_values, a1_values, a2_values) {

  idx <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))

  p_qi <- prob_values[iter, idx]
  a1_qi <- a1_values[iter, idx]
  a2_qi <- a2_values[iter, idx]

  K_i <- length(idx)
  out <- numeric(K_i)

  # vectorized entropy, handling 0*log(0) = 0 by convention
  #term a = plogp
  term_a = ifelse(p_qi>0, p_qi*log(p_qi), 0)

  #term b = (1-p)log(1-p)
  term_b = ifelse(p_qi <1, (1-p_qi)*log(1-p_qi), 0)

  out <- p_qi  * expectedLogTheta(a1_qi, a2_qi) -
      term_a +
      (1 - p_qi) * expectedLogOneMinusTheta(a1_qi, a2_qi) -
      term_b

  sum(out)
}

#beta coefficient term
elboBetaTerm <- function(i, iter, K,
                         delta1_q, delta2_values,
                         lambda1_q, lambda2_values,
                         mu_beta_values, Sigma_beta) {

  #delta and lambda values at iteration iter
  delta2_q <- delta2_values[iter]
  lambda2_q <- lambda2_values[iter]

  #identify beta cols for curve i
  idx <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))

  mu_i_q <- mu_beta_values[iter, idx]
  Sigma_i_q <- Sigma_beta[,,i]
  sigma2_i_q <- diag(Sigma_i_q)

  #using modulus based log determinant
  #beta-term contribution to the ELBO
  res <- -(K/ 2) * (expectedLogSigma2(delta1_q, delta2_q) + expectedLogTau2(lambda1_q, lambda2_q)) -
    sum((sigma2_i_q + mu_i_q^2) / 2) *
    expectedInvSigma2(delta1_q, delta2_q) * expectedInvTau2(lambda1_q, lambda2_q) -
    (-0.5 * as.numeric(determinant(Sigma_i_q, logarithm=TRUE)$modulus) - 0.5 * K)

  as.numeric(res)
}

#Theta (Beta prior) term
elboThetaTerm <- function(i, iter, K, a1_values, a2_values, mu_beta_values, mu_ki) {
  idx <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))
  a1_qi <- a1_values[iter, idx]
  a2_qi <- a2_values[iter, idx]

  # vectorized; lgamma used for numerical stability with large parameters
  term_sum <- (mu_ki - a1_qi) * expectedLogTheta(a1_qi, a2_qi) +
      (1 - mu_ki - a2_qi) * expectedLogOneMinusTheta(a1_qi, a2_qi) +
      lgamma(a1_qi) + lgamma(a2_qi)

  sum(term_sum)
}

#variance term
elboSigmaPriorTerm <- function(iter, delta_1, delta_2, delta1_q, delta2_values) {

  delta2_q <- delta2_values[iter]

  res <- -delta1_q * log(delta2_q) +
    (delta1_q - delta_1) * expectedLogSigma2(delta1_q, delta2_q) +
    (delta2_q - delta_2) * expectedInvSigma2(delta1_q, delta2_q)

  as.numeric(res)
}

#tau variance term
elboTauTerm <- function(iter, lambda_1, lambda_2, lambda1_q, lambda2_values) {

  lambda2_q <- lambda2_values[iter]

  res <- -lambda1_q * log(lambda2_q) +
    (lambda1_q - lambda_1) * expectedLogTau2(lambda1_q, lambda2_q) +
    (lambda2_q - lambda_2) * expectedInvTau2(lambda1_q, lambda2_q)

  as.numeric(res)
}

#elb_corr() - computes the full ELBOP for the VEM algorithm under correlated errors
elbo_corr <-  function(
    y, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values,
    lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki, psi
) {

  # expected log likelihood term
  term_loglik <- sum(
    sapply(1:m, function(i)
      expectedLogLikelihood(
        y = y, ni = ni, B = B, i = i, iter = iter,
        delta1_q = delta1_q, delta2_values = delta2_values,
        mu_beta_values = mu_beta_values,
        Sigma_beta = Sigma_beta,
        prob_values = prob_values,
        psi = psi
      )
    )
  )

  # inclusion (Z) indicator term
  term_z <- sum(
    sapply(1:m, function(i)
      elboInclusionTerm(
        i = i, iter = iter, K = K,
        prob_values = prob_values,
        mu_beta_values = mu_beta_values,
        a1_values = a1_values, a2_values = a2_values
      )
    )
  )

  # beta coefficient term
  term_beta <- sum(
    sapply(1:m, function(i)
      elboBetaTerm(
        i = i, iter = iter, K = K,
        delta1_q = delta1_q, delta2_values = delta2_values,
        lambda1_q = lambda1_q, lambda2_values = lambda2_values,
        mu_beta_values = mu_beta_values,
        Sigma_beta = Sigma_beta
      )
    )
  )

  # theta prior term
  term_theta <- sum(
    sapply(1:m, function(i)
      elboThetaTerm(
        i = i, iter = iter, K = K,
        a1_values = a1_values, a2_values = a2_values,
        mu_beta_values = mu_beta_values,
        mu_ki = mu_ki
      )
    )
  )

  # sigma term
  term_sigma <- elboSigmaPriorTerm(
    iter = iter, delta_1 = delta_1, delta_2 = delta_2,
    delta1_q = delta1_q, delta2_values = delta2_values
  )

  # tau term
  term_tau <- elboTauTerm(
    iter = iter, lambda_1 = lambda_1, lambda_2 = lambda_2,
    lambda1_q = lambda1_q, lambda2_values = lambda2_values
  )

  # summing up to arrive at total ELBO
  total <- term_loglik + term_z + term_beta + term_theta + term_sigma + term_tau
  as.numeric(total)
}

# elbo_omega - ELBO as function of w only (used during M-step)
elbo_omega <- function(
    w, y, Xt, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values, lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki
) {

  # compute corr mat
  psi <- computePsiMatrix(Xt, Xt, w)

  # compute full ELBO at this w
  elbo_value <- elbo_corr(
    y = y,
    B = B,
    ni = ni,
    m = m, K = K,
    iter = iter,
    delta_1 = delta_1, delta_2 = delta_2,
    lambda_1 = lambda_1, lambda_2 = lambda_2,
    delta1_q = delta1_q,
    delta2_values = delta2_values,
    mu_beta_values = mu_beta_values,
    lambda1_q = lambda1_q,
    lambda2_values = lambda2_values,
    a1_values = a1_values, a2_values = a2_values,
    Sigma_beta = Sigma_beta,
    prob_values = prob_values,
    mu_ki = mu_ki,
    psi = psi
  )

  as.numeric(elbo_value)
}

# dev_elbo - derivative of ELBO wrt w
dev_elbo <- function(
    w, y, Xt, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values, lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki
) {

  # variational inverse sigma2 expectation
  delta2_q <- delta2_values[iter]
  E_inv_sigma2 <- expectedInvSigma2(delta1_q, delta2_q)

  #compute Psi(w) and inverse
  psi <- computePsiMatrix(Xt, Xt, w)
  Psi_inv <- solve(psi)

  #computing derivatives wrt w
  dPsi_list <- devPsi(Xt, Xt, w)
  dPsi  <- dPsi_list$dPsi
  #d2Psi <- dPsi_list$d2Psi

  #building A_i matrics
  n_pts = length(y[[1]])
  # accumulate A_sum in-place to avoid allocating m large matrices
  A_sum <- matrix(0, nrow=n_pts, ncol=n_pts)

  for (i in 1:m) {
    idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))

    p_i  <- prob_values[iter, idx_i]
    mu_i <- mu_beta_values[iter, idx_i]
    Sig_i <- Sigma_beta[,,i]

    #bernoulli variance
    var_z <- diag(p_i*(1-p_i))

    Var_zi_betai <-var_z * Sig_i +var_z * (mu_i %*% t(mu_i)) +Sig_i * (p_i %*% t(p_i))

    # mean residual
    fitted <- t(p_i * t(B[[i]])) %*% mu_i
    resid  <- y[[i]] - fitted

    A_sum <- A_sum + (resid %*% t(resid) + B[[i]] %*% Var_zi_betai %*% t(B[[i]]))
  }

  # derivative of ELBO wrt w
  # CHANGE using trace properties
  term1 <- -(m / 2) * sum(Psi_inv * dPsi)
  Q_dPsi_Q <- Psi_inv %*% dPsi %*% Psi_inv
  term2 <-  0.5 * E_inv_sigma2 * sum(Q_dPsi_Q * A_sum)

  as.numeric(term1+term2)
}


