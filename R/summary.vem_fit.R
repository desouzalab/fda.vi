#' Summary Method for VEM Fits
#'
#' @description
#' Provides a displayed summary of the results from \code{\link{vem_fit}} and
#' invisibly returns a list of summary statistics, including the basis type,
#' number of curves, selected \eqn{K}, active basis counts per curve,
#' estimated model parameters, and GCV tuning results if applicable.
#'
#' Reported variational posterior parameters for \eqn{\sigma^2} and
#' \eqn{\tau^2} are the shape and scale of their respective
#' Inverse-Gamma variational distributions:
#' \eqn{q(\sigma^2) = \text{IG}(\delta_1^*, \delta_2^*)} and
#' \eqn{q(\tau^2) = \text{IG}(\lambda_1^*, \lambda_2^*)}.
#' For composite fits (\code{selection_metric = "per_curve"}), parameters
#' from the first curve are shown as representative values.
#'
#' @param object A \code{vem_fit} object from \code{\link{vem_fit}}.
#' @param ... Currently unused.
#'
#' @return Invisibly returns a list with element \code{active_bases}: an
#'   integer vector of active basis counts per curve.
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{coef.vem_fit}}
#' @export
#' @examples
#' \donttest{
#'   data(toy_curves)
#'   fit <- vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8)
#'
#'   summary(fit)
#'
#'   # Active basis counts are returned invisibly
#'   s <- summary(fit)
#'   s$active_bases
#' }

summary.vem_fit <- function(object, ...) {
  m <- length(object$data_orig)
  b_type <- object$basis_type

  # active counts and representative parameter extraction
  if (object$is_composite) {
    active_counts <- integer(m)
    K_vals <- integer(m)

    for (i in 1:m) {
      mod <- object$model[[i]]
      active_counts[i] <- sum(mod$prob > 0.5)
      K_vals[i] <- mod$K
    }
    K_display <- paste0("Varies (", min(K_vals), "-", max(K_vals), ")")

    # representative parameters: fit at most frequently selected K (smallest on tie)
    modal_K  <- as.numeric(names(which.max(table(object$best_K))))
    rep_fit  <- object$tuning$fits[[which(eval(object$call$K) == modal_K)]]
    delta1_val  <- rep_fit$delta1
    delta2_val  <- rep_fit$delta2
    lambda1_val <- rep_fit$lambda1
    lambda2_val <- rep_fit$lambda2
    w           <- rep_fit$w

  } else {
    K <- object$best_K
    K_display <- as.character(K)

    # count active bases per curve
    prob_vec <- object$model$prob
    prob_mat <- matrix(prob_vec, nrow = m, ncol = K, byrow = TRUE)
    active_counts <- rowSums(prob_mat > 0.5)

    # representative parameters: fit at best_K when tuning, else model directly
    if (!is.null(object$tuning)) {
      rep_fit     <- object$tuning$fits[[which(eval(object$call$K) == object$best_K)]]
      delta1_val  <- rep_fit$delta1
      delta2_val  <- rep_fit$delta2
      lambda1_val <- rep_fit$lambda1
      lambda2_val <- rep_fit$lambda2
      w           <- rep_fit$w
    } else if (!is.null(object$model$delta2)) {
      delta1_val  <- object$model$delta1
      delta2_val  <- object$model$delta2
      lambda1_val <- object$model$lambda1
      lambda2_val <- object$model$lambda2
      w           <- object$model$w
    } else {
      delta1_val <- delta2_val <- lambda1_val <- lambda2_val <- w <- NA
    }
  }

  # summary table
  cat("------------------------------------------------\n")
  cat("  VEM Smooth Fit Summary\n")
  cat("------------------------------------------------\n")
  cat("Basis Type           ", b_type, "\n")
  cat("Curves (m):          ", m, "\n")
  cat("Basis Functions (K): ", K_display, "\n")

  if(m <= 10) {
    cat("Active Bases (per curve):", paste(active_counts, collapse=", "), "\n")
  } else {
    cat("Active Bases (summary):  Min:", min(active_counts),
        " Median:", median(active_counts),
        " Max:", max(active_counts), "\n")
  }

  if (!is.na(w)) {
    cat("\nModel Parameters (Representative):\n")
    cat("  Point estimate for decay parameter (w):", format(w, scientific=FALSE, digits = 4), "\n")
    cat("\n  Posterior q(sigma^2) ~ IG(delta1, delta2):\n")
    cat("    delta1 (shape): ", format(delta1_val,  scientific=FALSE, digits = 6), "\n")
    cat("    delta2 (scale): ", format(delta2_val,  scientific=FALSE, digits = 6), "\n")
    cat("\n  Posterior q(tau^2) ~ IG(lambda1, lambda2):\n")
    cat("    lambda1 (shape):", format(lambda1_val, scientific=FALSE, digits = 6), "\n")
    cat("    lambda2 (scale):", format(lambda2_val, scientific=FALSE, digits = 6), "\n")
  }

  # gcv results
  if (!is.null(object$tuning)) {
    cat("\nGCV Tuning Results (Mean Score):\n")
    gcv_mat <- object$tuning$gcv_matrix
    if (!is.null(gcv_mat)) {
      mean_scores <- colMeans(gcv_mat, na.rm = TRUE)
      print(round(mean_scores, 4))
    }
  }

  cat("------------------------------------------------\n")
  invisible(list(active_bases = active_counts))
}
