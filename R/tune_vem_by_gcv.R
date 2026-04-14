#' Tune Basis Complexity via GCV
#'
#' @description
#' Fits \code{\link{vem_smooth}} across a grid of candidate basis sizes
#' \code{K_grid} and selects the best \eqn{K} using GCV scores from
#' \code{\link{gcv_vem}}. Called internally by \code{\link{vem_fit}} when
#' a vector of \eqn{K} values is supplied; not typically called directly.
#'
#' Two selection modes are supported: \code{"mean"} selects the single \eqn{K}
#' minimizing the mean GCV across all curves; \code{"per_curve"} selects the
#' \eqn{K} that minimizes the GCV criterion for each individual curve, producing a composite fit.
#'
#' @param y List of curves.
#' @param Xt Numeric vector of time points.
#' @param K_grid Integer vector of candidate \eqn{K} values.
#' @param build_B Function with signature \code{function(K, Xt, y)} that returns
#'   a list of \eqn{n \times K} basis matrices.
#' @param initial_values_fn Function with signature \code{function(K, m)} that
#'   returns an \code{initial_values} list for \code{\link{vem_smooth}}.
#' @param threshold Posterior inclusion probability (PIP) threshold passed to \code{\link{gcv_vem}}. Default \code{0.5}.
#' @param mode Character. \code{"mean"} for a single global \eqn{K};
#'   \code{"per_curve"} for curve-specific \eqn{K}. Default \code{"mean"}.
#' @param ... Additional arguments passed to \code{\link{vem_smooth}}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{fits}}{Named list of fitted \code{vem_smooth} objects, one per
#'     candidate \eqn{K}.}
#'   \item{\code{gcv_matrix}}{Numeric matrix (\eqn{m \times} \code{length(K_grid)})
#'     of per-curve GCV scores.}
#'   \item{\code{best_K_mean}}{Integer. Best \eqn{K} by mean GCV.}
#'   \item{\code{best_K_per_curve}}{Integer vector of length \eqn{m}. Best
#'     \eqn{K} for each curve.}
#' }
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{gcv_vem}}
#' @export

tune_vem_by_gcv <- function(
    y,
    Xt,
    K_grid,
    build_B,
    initial_values_fn,
    threshold = 0.5,
    mode = c("mean", "per_curve"),
    ...
) {

  mode <- match.arg(mode)
  m <- length(y)

  # fits = models and gcv scores storage
  fits <- vector("list", length(K_grid))
  names(fits) <- paste0("K_", K_grid)

  # GCV matrix - rows = curves, cols = k vals
  gcv_mat <- matrix(NA_real_, nrow=m, ncol=length(K_grid))
  colnames(gcv_mat) <- paste0("K_", K_grid)

  # grid search through K
  for (j in seq_along(K_grid)) {
    K_curr <- K_grid[j]

    # basis functions
    B_list <- build_B(K = K_curr, Xt = Xt, y = y)
    K_effective <- ncol(B_list[[1]])

    init_vals <- initial_values_fn(K = K_effective, m = m)

    # fit model
    fit <- vem_smooth(y = y, B = B_list, Xt = Xt, m = m, K = K_effective,
                      initial_values = init_vals, ...)

    # store fit and k values
    fit$K <- K_effective
    fits[[j]] <- fit

    # calc gcv and store
    gcv_vec <- gcv_vem(fit, threshold = threshold)
    gcv_mat[, j] <- gcv_vec
  }

  # selection logic

  # average gcv across all curves
  mean_scores <- colMeans(gcv_mat, na.rm = TRUE)
  best_idx_mean <- which.min(mean_scores)

  # per curve gcv
  best_idx_per_curve <- apply(gcv_mat, 1, which.min)
  best_K_per_curve <- K_grid[best_idx_per_curve]

  # final model returned
  list(
    fits = fits,
    gcv_matrix = gcv_mat,
    best_K_mean = K_grid[best_idx_mean],
    best_K_per_curve = best_K_per_curve
  )
}
