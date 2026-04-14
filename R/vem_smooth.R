#' Variational EM Algorithm for Bayesian Basis Function Selection
#'
#' @description
#' Fits \eqn{m} functional curves simultaneously via Bayesian basis function
#' selection with an Ornstein-Uhlenbeck within-curve correlation structure.
#' This function is called internally by \code{\link{vem_fit}} and only runs
#' the VEM algorithm itself, without performing basis construction,
#' standardization, or GCV tuning. Most users should call
#' \code{\link{vem_fit}} instead, which handles those steps automatically.
#'
#' @details
#' The algorithm alternates between an E-step — sequential coordinate ascent variational inference (CAVI) updates for
#' \eqn{q(\beta_i)}, \eqn{q(\sigma^2)}, \eqn{q(\tau^2)}, \eqn{q(Z_{ki})},
#' and \eqn{q(\theta_{ki})} — and an M-step that maximizes the ELBO with
#' respect to the correlation decay parameter \eqn{w} via L-BFGS-B with an
#' analytic gradient. Convergence is declared when the absolute ELBO change
#' between iterations falls below \code{convergence_threshold}.
#'
#' For hyperparameter initialization, set \code{delta_1} and \code{delta_2}
#' such that \code{delta_2 / (delta_1 - 1)} is a rough estimate of the noise
#' variance, and initialize \code{w} consistent with the expected correlation
#' strength in the data.
#'
#' @param y List of length \eqn{m} of numeric vectors (observed curves,
#'   possibly standardized).
#' @param B List of length \eqn{m} of \eqn{n_i \times K} basis matrices,
#'   typically from \code{\link[fda]{getbasismatrix}}.
#' @param Xt Numeric vector of \eqn{n} evaluation points, common across curves.
#' @param m Integer. Number of curves. Defaults to \code{length(y)}.
#' @param K Integer. Number of basis functions.
#' @param mu_ki Numeric scalar in \eqn{(0,1)}. Beta prior hyperparameter for
#'   inclusion probabilities. Default \code{0.5}.
#' @param lambda_1,lambda_2 Positive scalars. Inverse-Gamma prior hyperparameters
#'   for \eqn{\tau^2}. Default \eqn{10^{-10}}.
#' @param delta_1,delta_2 Positive scalars. Inverse-Gamma prior hyperparameters
#'   for \eqn{\sigma^2}. Default \eqn{10^{-10}}.
#' @param maxIter Integer. Maximum VEM iterations. Default \code{1000}.
#' @param initial_values Named list with elements \code{p} (inclusion
#'   probabilities, length \eqn{mK}), \code{delta2}, \code{lambda2}, and
#'   \code{w}.
#' @param convergence_threshold Positive scalar. Absolute ELBO tolerance for
#'   convergence. Default \code{0.01}.
#' @param lower_opt Positive scalar. Lower bound for \eqn{w} in L-BFGS-B.
#'   Default \code{0.1}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{mu_beta}}{Posterior means \eqn{\mu_{\beta_{ki}}} (length \eqn{mK}).}
#'   \item{\code{Sigma_beta}}{Posterior covariance array (\eqn{K \times K \times m}).}
#'   \item{\code{prob}}{Posterior inclusion probabilities \eqn{p_{ki}} (length \eqn{mK}).
#'     Basis \eqn{k} is active for curve \eqn{i} when \eqn{p_{ki} > 0.5}.}
#'   \item{\code{delta1}, \code{delta2}}{Final \eqn{q(\sigma^2)} parameters.}
#'   \item{\code{lambda1}, \code{lambda2}}{Final \eqn{q(\tau^2)} parameters.}
#'   \item{\code{w}}{Estimated correlation decay parameter (range-normalized scale).}
#'   \item{\code{cor_mat}}{The \eqn{n \times n} Ornstein-Uhlenbeck correlation
#'     matrix \eqn{\Psi} evaluated at the final estimated decay parameter
#'     \eqn{\hat{w}}, as returned by \code{computePsiMatrix}.}
#'   \item{\code{elbo_values}}{ELBO trajectory across iterations.}
#'   \item{\code{converged}}{Logical. Whether the convergence criterion was met.}
#'   \item{\code{n_iterations}}{Number of iterations run.}
#' }
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{plot.vem_fit}},
#'   \code{\link{predict.vem_fit}}, \code{\link{coef.vem_fit}}
#' @export

vem_smooth <- function(
    y,
    B,
    Xt = Xt,
    m = length(y),
    K = K,
    mu_ki = 0.5,
    lambda_1 = 1e-10, lambda_2 = 1e-10,
    delta_1 = 1e-10, delta_2 = 1e-10,
    maxIter=1000,
    initial_values,
    convergence_threshold = 0.01,
    lower_opt = 0.1
) {

  converged <- FALSE
  iter <- 1
  ni <- unlist(lapply(y,length)) # stores sample sizes per curve

  #intialize w and covariance matrix (chosen by user) --> force positivity to prevent garbage values
  w_decay <- initial_values$w
  w_decay <- max(lower_opt, as.numeric(initial_values$w))

  Psi_mat <- computePsiMatrix(Xt, Xt, w=w_decay) #correlation kernel

  #storage for the matrices for variational parameters / iterative updates
  mu_beta_iter <- matrix(NA, nrow=maxIter, ncol = m*K) #posterior beta means (Norm dist)
  colnames(mu_beta_iter) <- as.vector(outer(paste0("beta_", 1:K, "_"), 1:m, paste0))

  Sigma_beta <- array(NA, dim=c(K, K, m)) #posterior beta covariances (inv gamma) (for each curve)
  prob_inclusion_iter <- matrix(NA, nrow=maxIter, ncol = m*K) # inclusion probabilities
  a1_iter <- a2_iter <- matrix(NA, nrow=maxIter, ncol=m*K) # theta parameters (Beta dist)
  tau2_var_iter <- sigma2_var_iter <- numeric(maxIter) # (rate) parameters for tau, sigma (inv gamma)
  # (w_trace removed; final correlation matrix returned as cor_mat)

  #Initialize (chosen by user)
  prob_inclusion_iter[1, ] <- initial_values$p
  tau2_var_iter[1] <- initial_values$lambda2
  sigma2_var_iter[1] <- initial_values$delta2

  # fixed variational shape parameters (ie. not getting updated in VEM --> are fixed)
  delta1_q <- (sum(ni)+ m*K + 2 * delta_1)/2
  lambda1_q <- (m*K+2*lambda_1)/2

  #initialize ELBO tracking
  elbo_prev <- 0
  ELBO_values <- -Inf

  ### MAIN VEM LOOP
  while (!converged && iter < maxIter) {
    iter <- iter +1

    # E-step
    ### Step 1 - Update q(Bi)

    # precompute Psi inverse once — common to all curves at this iteration
    Psi_inv <- solve(Psi_mat)

    for (i in 1:m) {
      #inclusion probabilities for curve i from prev iter
      idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_iter))
      prob_i   <- prob_inclusion_iter[iter - 1, idx_i]

      #Expected basis contributions under q(zi)

      #expected basis design
      Gi <- prob_i * t(B[[i]])
      GPsiGt <- Gi %*% Psi_inv %*% t(Gi)
      yPsiG <- t(y[[i]]) %*% Psi_inv %*% t(Gi)

      #Expected inverse variances: under variational posteriors
      E_tau_inv <- expectedInvTau2(lambda1_q, tau2_var_iter[iter-1])
      E_sigma_inv <- expectedInvSigma2(delta1_q, sigma2_var_iter[iter-1])

      #variational covar of Bi
      V_i <- solve(diag(E_tau_inv, K)+GPsiGt)

      #full covariance and mean of q(Bi)
      Sigma_bi <- (1/E_sigma_inv) * V_i
      mu_bi <- (E_sigma_inv * yPsiG) %*% Sigma_bi

      #storing updated variational parameters for curve i
      mu_beta_iter[iter, idx_i] <- mu_bi
      Sigma_beta[,,i] <- Sigma_bi
    }

    ### Step 2 - update q(sigma) - scale parameter for delta2_q (stored in sigma2_var_iter)
    # CHANGE using vapply (marginally faster and safer)
    sigma2_var_iter[iter] <- (
      sum(vapply(1:m, function(i){
        expectedResidualSq(B=B, i=i, y=y, mu=mu_beta_iter, Sigma=Sigma_beta, prob=prob_inclusion_iter, iter=iter, Psi_mat)
      }, numeric(1))) +
        sum(vapply(1:m, function(i){
          expectedSumBetaSq(i, mu_beta_iter, Sigma_beta, iter)
        }, numeric(1))) * expectedInvTau2(lambda1_q, tau2_var_iter[iter-1]) + 2 * delta_2
    ) / 2


    #Step 3 - update q(tau): scale parameter for lambda2_q (stored in tau2_var_iter)

    tau2_var_iter[iter] <- (
      sum(vapply(1:m, function(i){
        expectedSumBetaSq(i, mu_beta_iter, Sigma_beta, iter)
      }, numeric(1))) * expectedInvSigma2(delta1_q, sigma2_var_iter[iter]) + 2 * lambda_2
    ) / 2

    #Step 4 - update q(Z_Ki) and q(theta_Ki)
    for (i in 1:m) {
      idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_iter))

      # probabilities from previous iteration
      prob_prev <- prob_inclusion_iter[iter - 1, idx_i]

      prob_curr <-  prob_inclusion_iter[iter, idx_i]

      #fill NAs in curr row with previous iteration values
      na_idx <- is.na(prob_curr)
      if (any(na_idx)) prob_curr[na_idx] <- prob_prev[na_idx]

      #storage for updated values at this iter
      prob_new <- numeric(K)
      a1_new <- numeric(K)
      a2_new <- numeric(K)


      for (k in 1:K) {
        #variational parameters for theta_ki
        a1_q <- prob_prev[k] + mu_ki
        a2_q <- 2 - prob_prev[k] - mu_ki

        # CHANGE - (1-z) logic
        # unnormalized log posterior for Zki
        log_rho <- vapply(0:1, function(z) {
          (-ni[i] / 2) * expectedLogSigma2(delta1_q, sigma2_var_iter[iter]) -
            0.5 * expectedInvSigma2(delta1_q, sigma2_var_iter[iter]) *
            expectedBetaSq(z, i, prob_prev, mu_beta_iter, Sigma_beta, B, y, k, K, iter, Psi_mat) +
            z * expectedLogTheta(a1_q, a2_q) +
            (1-z) * expectedLogOneMinusTheta(a1_q, a2_q)
        }, numeric(1))

        # convert log_rho for z=0,1 into prob p(Zki)
        denom <- sum(exp(log_rho))
        #fallback if numerical issues
        if (denom == 0 || is.infinite(denom)) {
          prob_ki_1 <- as.numeric(which.max(log_rho) == 2) # 0 if max is index 1, 1 if index 2
        } else {
          prob_ki_1 <- exp(log_rho[2]) / denom
        }

        #coordinate ascent update
        prob_curr[k] <- prob_ki_1

        #storing output
        prob_new[k] <- prob_ki_1
        a1_new[k] <- a1_q
        a2_new[k] <- a2_q
      }

      #assigning updated values
      prob_inclusion_iter[iter, idx_i] <-prob_curr
      a1_iter[iter, idx_i] <- a1_new
      a2_iter[iter, idx_i] <- a2_new
    }

    w_c <- "Error"
    current_lower <- lower_opt

    while (is.character(w_c) && w_c == "Error") {
      w_result <- tryCatch(
        optim(
          par=w_decay,
          fn = elbo_omega,
          gr=dev_elbo,
          method='L-BFGS-B',
          lower=current_lower,
          upper=1e10,
          control=list(fnscale=-1), # maximize ELBO

          #inputs
          y=y, Xt=Xt, B=B, ni=ni, m=m, K=K, iter=iter,
          delta_1 = delta_1, delta_2 = delta_2,
          lambda_1 = lambda_1, lambda_2 = lambda_2,
          delta1_q = delta1_q, delta2_values = sigma2_var_iter,
          mu_beta_values = mu_beta_iter,
          lambda1_q = lambda1_q, lambda2_values = tau2_var_iter,
          a1_values = a1_iter, a2_values = a2_iter,
          Sigma_beta = Sigma_beta,
          prob_values = prob_inclusion_iter,
          mu_ki = mu_ki
        ),
        error = function(e) "Error"
      )

      if(is.list(w_result)) {
        w_c <- w_result$par #extract the parameter value

        if (!is.finite(w_c)) w_c <- "Error"
      } else {
        w_c <- "Error"
      }

      #if all fails, increase lower bound
      if (is.character(w_c) && w_c == "Error") {
        current_lower <- current_lower+0.5

        if (current_lower > 50) {
          warning("Optimization for w failed repeatedly, keeping previous w.")
          w_c <- w_decay
          break
        }
    }
  }

    w_decay <- as.numeric(w_c)



    Psi_mat <- computePsiMatrix(Xt, Xt, w_decay)

    #Step 6 Compute ELBO
    elbo_curr <-  elbo_corr(
      y, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2,
      delta1_q, sigma2_var_iter, mu_beta_iter, lambda1_q, tau2_var_iter,
      a1_iter, a2_iter, Sigma_beta, prob_inclusion_iter, mu_ki, Psi_mat
    )

    ELBO_values <- c(ELBO_values, elbo_curr)

    converged <- checkConvergence(elbo_curr, elbo_prev, convergence_threshold)
    elbo_prev <- elbo_curr
  }

  #OUTPUT
  result <- list(
    data = y, B = B, Xt=Xt, m=m, K=K,
    mu_beta = mu_beta_iter[iter, ], Sigma_beta=Sigma_beta,
    prob=prob_inclusion_iter[iter, ],
    lambda1=lambda1_q, lambda2 = tau2_var_iter[iter],
    delta1=delta1_q, delta2=sigma2_var_iter[iter],
    a1=a1_iter[iter, ], a2 = a2_iter[iter, ],
    elbo_values = ELBO_values, n_iterations = iter,
    converged =converged, w=w_decay, cor_mat=Psi_mat
  )

  #ALWAYS return best guess
  if (!converged) {
    warning("The algorithm did not converge within the maximum iterations...")
  }
  return(result)
}
