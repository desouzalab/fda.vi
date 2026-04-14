#' Plot a VEM Fit with Credible Band
#'
#' @description
#' Plots observed data, the posterior mean fitted curve, and an optional 95%
#' credible band for a single curve from a \code{vem_fit} object. The credible
#' band provides uncertainty quantification by sampling from the
#' variational posteriors: \eqn{\beta_i \sim \text{MVN}(\boldsymbol{\mu}_{\boldsymbol{\beta}_i}, \boldsymbol{\Sigma}_{\boldsymbol{\beta}_i})} and
#' \eqn{Z_{ki} \sim \text{Bernoulli}(p_{ki})}.
#' Predictions are automatically back-transformed if the model was fitted with
#' \code{center = TRUE} or \code{scale = TRUE}.
#'
#' @param x A \code{vem_fit} object from \code{\link{vem_fit}}.
#' @param curve_idx Integer. Index of the curve to plot. Default \code{1}.
#' @param type Character. Credible band style: \code{"polygon"} (shaded region)
#'   or \code{"lines"} (dashed lines). Default \code{"polygon"}.
#' @param show_CI Logical. If \code{TRUE}, compute and display the credible band.
#'   Default \code{TRUE}.
#' @param n_samples Integer. Number of posterior draws used to construct the
#'   credible band. Default \code{200}.
#' @param alpha_shade Numeric in \eqn{(0,1)}. Opacity of the shaded credible
#'   band (\code{type = "polygon"} only). Default \code{0.25}.
#' @param ylim Optional numeric vector of length 2. If \code{NULL}, axis limits
#'   are set to cover the data and credible band.
#' @param xlab Character. Label for the horizontal axis. Default \code{"t"}
#'   (the evaluation point). Supply a domain-specific label where appropriate
#'   (e.g., \code{"Time"}, \code{"Age"}, \code{"Frequency"}).
#' @param ylab Character. Label for the vertical axis. Default \code{"Value"}.
#' @param show_basis Logical. If \code{TRUE}, adds a subplot below the main
#'   plot showing each basis function coloured by inclusion status (blue =
#'   active, grey = inactive). Default \code{FALSE}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a plot.
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{predict.vem_fit}}
#' @export
#' @examples
#' \donttest{
#'   data(toy_curves)
#'   fit <- vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8)
#'
#'   # Default: shaded credible band for curve 1
#'   plot(fit)
#'
#'   # Dashed credible band for curve 2
#'   plot(fit, curve_idx = 2, type = "lines")
#'
#'   # With basis selection subplot
#'   plot(fit, curve_idx = 1, show_basis = TRUE)
#'
#'   # Suppress credible band
#'   plot(fit, show_CI = FALSE, main = "Mean fit only")
#' }

plot.vem_fit <- function(x, curve_idx=1, type=c("polygon", "lines"),
                         show_CI = TRUE, n_samples=200, alpha_shade=0.25,
                         ylim=NULL, xlab="t", ylab="Value",
                         show_basis = FALSE, ...) {

  type <- match.arg(type)

  # Data validation
  if(is.null(x$data_orig)) stop("Original data not found in fit object.")
  y_raw <- x$data_orig[[curve_idx]]
  Xt <- x$Xt

  # Extract parameters based on model structure
  if (isTRUE(x$is_composite)) {
    mod <- x$model[[curve_idx]]
    K <- mod$K
    mu_i <- mod$mu_beta
    prob_i <- mod$prob
    Sigma_i <- mod$Sigma_beta
    B_i <- mod$B
  } else {
    mod <- x$model
    K <- if(length(x$selected_K) > 1) x$selected_K[curve_idx] else x$best_K

    idx_start <- (curve_idx - 1) * K + 1
    idx_end   <- curve_idx * K

    mu_i <- mod$mu_beta[idx_start:idx_end]
    prob_i <- mod$prob[idx_start:idx_end]
    Sigma_i <- mod$Sigma_beta[,,curve_idx]
    B_i <- mod$B[[curve_idx]]
  }

  # Mean fit (Slab selection)
  mu_active <- mu_i * (prob_i > 0.5)
  y_hat_std <- as.vector(B_i %*% mu_active)

  # Destandardize logic (Safe for center = FALSE)
  if (!is.null(x$scaling_params)) {
    s_sd <- x$scaling_params$sds[[curve_idx]]
    # Check if 'center' was TRUE in the original call
    # If x$call$center isn't available, we check if means exist
    s_mean <- if(!is.null(x$scaling_params$means)) x$scaling_params$means[[curve_idx]] else 0
    y_hat <- (y_hat_std * s_sd) + s_mean
  } else {
    y_hat <- y_hat_std
  }

  # Simulation for Credible Bands
  LL <- NULL; UL <- NULL
  if (show_CI) {
    betas_sample <- MASS::mvrnorm(n_samples, mu = mu_i, Sigma = Sigma_i)
    # Spike-and-Slab Simulation
    z_sample <- matrix(rbinom(K * n_samples, 1, prob_i), nrow = n_samples, ncol = K, byrow = TRUE)
    sim_curves <- (betas_sample * z_sample) %*% t(B_i)

    if (!is.null(x$scaling_params)) {
      sim_curves <- (sim_curves * s_sd) + s_mean
    }

    LL <- apply(sim_curves, 2, quantile, probs=0.025)
    UL <- apply(sim_curves, 2, quantile, probs=0.975)
  }

  # Set axis limits
  if (is.null(ylim)) {
    ylim <- range(c(y_raw, y_hat, LL, UL), na.rm = TRUE)
  }

  # Layout for basis plot
  if (show_basis) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    layout(matrix(c(1, 2), 2, 1), heights = c(2, 1))
    par(mar = c(3, 4, 2, 1))
  }

  # Main Plot
  dot_args <- list(...)
  if (!"main" %in% names(dot_args)) dot_args$main <- paste0("VEM Fit: ", names(x$data_orig)[curve_idx])

  do.call(plot, c(list(x = Xt, y = y_raw, type="n", ylim = ylim, xlab = if(show_basis) "" else xlab, ylab = ylab), dot_args))

  # Shaded Credible Band
  if (show_CI) {
    if (type == "polygon") {
      polygon(c(Xt, rev(Xt)), c(LL, rev(UL)), col = scales::alpha("grey70", alpha_shade), border = NA)
    } else {
      lines(Xt, LL, col = "grey40", lty = 2); lines(Xt, UL, col = "grey40", lty = 2)
    }
  }

  points(Xt, y_raw, pch = 16, col = scales::alpha("dodgerblue", 0.6), cex = 0.6)
  lines(Xt, y_hat, lwd = 2.5, col = "firebrick")

  legend("topright", bty = "n", cex = 0.8,
         legend = c("Observed", "Mean Fit", "95% Credible Band"),
         col = c("dodgerblue", "firebrick", "grey70"),
         lty = c(NA, 1, 1), pch = c(16, NA, 15))

  # Basis Selection Subplot
  if (show_basis) {
    par(mar = c(4, 4, 1, 1))
    plot(Xt, rep(0, length(Xt)), type="n", ylim=range(B_i), xlab=xlab, ylab="Basis")
    for(k in 1:K) {
      col_b <- if(prob_i[k] > 0.5) "dodgerblue" else "grey80"
      lty_b <- if(prob_i[k] > 0.5) 1 else 3
      lines(Xt, B_i[, k], col=col_b, lty=lty_b, lwd = if(prob_i[k]>0.5) 1.5 else 1)
    }
  }
}
