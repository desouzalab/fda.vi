#' Fit a VEM Smooth Model
#'
#' @description
#' Fits one or more functional curves using Bayesian basis function selection
#' via the Variational EM algorithm, with an Ornstein-Uhlenbeck within-curve
#' correlation structure. Internally calls \code{\link{vem_smooth}} to run the
#' VEM algorithm.
#'
#' If a single value is provided for K, the model is fitted using that fixed
#' number of basis functions. If a vector of candidate values is supplied, the
#' function \code{\link{tune_vem_by_gcv}} is called to automatically select the
#' optimal K based on the GCV criterion. The resulting fitted object provides
#' methods for plot, predict, coef, and summary via the corresponding S3
#' methods \code{\link{plot.vem_fit}}, \code{\link{predict.vem_fit}},
#' \code{\link{coef.vem_fit}}, and \code{\link{summary.vem_fit}}.
#'
#' @param y Named list of numeric vectors, one per curve.
#' @param Xt Numeric vector of time points, common across all curves.
#' @param K Integer or integer vector of candidate basis sizes. If a single
#'   value, fits directly at that \code{K}. If a vector, selects best \code{K}
#'   via GCV. If \code{NULL}, defaults to \code{c(10, 15, 20, 30)} for B-splines
#'   and Fourier.
#' @param basis_type Character. One of \code{"cubic_bspline"} (default), or
#'   \code{"fourier"}.
#' @param selection_metric Character. \code{"mean"} selects a single global
#'   \code{K} minimizing mean GCV across curves. \code{"per_curve"} selects
#'   the best \code{K} independently for each curve, returning a composite fit.
#'   Only relevant when \code{K} is a vector. Default \code{"mean"}.
#' @param threshold Numeric in \code{(0,1)}. Posterior inclusion probability (PIP) threshold for active basis
#'   functions in GCV calculation. Default \code{0.5}.
#' @param center Logical. If \code{TRUE}, subtract each curve's mean before
#'   fitting. If TRUE, the function automatically centralizes the curves before
#'   model fitting. Default \code{FALSE}.
#' @param scale Logical. If \code{TRUE}, the function automatically standardizes
#'   the curves before model fitting, by dividing each curve by its standard
#'   deviation. Predictions are automatically back-transformed.
#'   Default \code{FALSE}.
#' @param period Numeric. Period for Fourier bases. Defaults to the domain
#'   range \code{diff(range(Xt))} if \code{NULL}.
#' @param initial_values_fn Function with signature \code{function(K, m)}
#'   returning an \code{initial_values} list for \code{\link{vem_smooth}}.
#'   If \code{NULL}, an empirical initialization based on a dense regression
#'   spline fit is used.
#' @param lambda_1,lambda_2 Positive scalars. Inverse-Gamma prior hyperparameters
#'   for \eqn{\tau^2}. Defaults: \code{lambda_1 = 0.001}, \code{lambda_2 = 0.001}.
#' @param delta_1,delta_2 Positive scalars. Inverse-Gamma prior hyperparameters
#'   for \eqn{\sigma^2}. Defaults: \code{delta_1 = 10}, \code{delta_2 = 0.09}.
#' @param ... Additional arguments passed to \code{\link{vem_smooth}}, such as
#'   \code{maxIter}, \code{convergence_threshold}, and \code{mu_ki}.
#'
#' @return An object of class \code{"vem_fit"} containing:
#' \describe{
#'   \item{\code{model}}{The fitted \code{vem_smooth} object (global fit), or a
#'     named list of per-curve model objects (composite fit).}
#'   \item{\code{selected_K}}{Integer vector of length \code{m}. The \code{K}
#'     used for each curve.}
#'   \item{\code{best_K}}{The single selected \code{K} (global fit), or a vector
#'     (composite fit).}
#'   \item{\code{tuning}}{Output of \code{\link{tune_vem_by_gcv}}, including
#'     the full GCV matrix and all candidate fits. \code{NULL} if a single
#'     \code{K} was supplied.}
#'   \item{\code{scaling_params}}{List with \code{means} and \code{sds} used
#'     for standardization. Used by \code{\link{predict.vem_fit}} and
#'     \code{\link{plot.vem_fit}} to back-transform predictions.}
#'   \item{\code{data_orig}}{The input curves in their original scales.}
#'   \item{\code{basis_type}, \code{is_composite}, \code{Xt}, \code{call}}{
#'     Metadata stored for use by S3 methods.}
#' }
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_smooth}}, \code{\link{plot.vem_fit}},
#'   \code{\link{predict.vem_fit}}, \code{\link{coef.vem_fit}},
#'   \code{\link{summary.vem_fit}}
#' @export
#' @examples
#' \donttest{
#'   data(toy_curves)
#'
#'   fit <- vem_fit(
#'     y    = toy_curves$y,
#'     Xt   = toy_curves$Xt,
#'     K    = 8
#'   )
#'
#'   summary(fit)
#'   plot(fit, curve_idx = 1)
#'   coef(fit)
#'   predict(fit)
#'
#'   # GCV tuning over a grid of K values
#'   fit_gcv <- vem_fit(
#'     y    = toy_curves$y,
#'     Xt   = toy_curves$Xt,
#'     K    = c(6, 8, 10)
#'   )
#'   fit_gcv$best_K
#' }
vem_fit <- function(y, Xt, K = NULL,
                    basis_type = c("cubic_bspline", "fourier"),
                    selection_metric = c("mean", "per_curve"),
                    threshold = 0.5,
                    center = FALSE, scale = FALSE,
                    period = NULL,
                    initial_values_fn = NULL,
                    lambda_1 = NULL, lambda_2 = NULL,
                    delta_1 = NULL, delta_2 = NULL,
                    ...) {

  # validating that y is a curve
  if (!is.list(y) || length(y) == 0) stop("Input 'y' must be a list of numeric vectors.")
  basis_type <- match.arg(basis_type)
  selection_metric <- match.arg(selection_metric)
  m <- length(y)

  # default hyperparams
  if (is.null(lambda_1)) lambda_1 <- 0.001
  if (is.null(lambda_2)) lambda_2 <- 0.001
  if (is.null(delta_1))  delta_1  <- 10
  if (is.null(delta_2))  delta_2  <- 0.09

  # standardizing daya
  curve_means <- lapply(y, mean, na.rm = TRUE)
  curve_sds <- lapply(y, sd, na.rm = TRUE)

  if (!center) curve_means <- lapply(seq_len(m), function(x) 0)
  if (!scale)  curve_sds   <- lapply(seq_len(m), function(x) 1)

  y_proc <- lapply(seq_len(m), function(i) {
    denom <- if (scale && curve_sds[[i]] != 0) curve_sds[[i]] else 1
    (y[[i]] - curve_means[[i]]) / denom
  })

  scaling_params <- list(means = curve_means, sds = curve_sds)

  # K grid set up
  if (is.null(K) || length(K) == 0) {
    K_grid <- c(10, 15, 20, 30)
  } else {
    K_grid <- sort(unique(K))
  }
  if (basis_type == "fourier") {
    K_adjusted <- ifelse(K_grid %% 2 != 0, K_grid + 1, K_grid)
    if (!identical(K_grid, K_adjusted)) K_grid <- sort(unique(K_adjusted))
  }
  message("\nRunning VEM Smooth | ")
  message(sprintf("Basis Type: %s | Selection: %s ---\n", basis_type, selection_metric))

  build_B_default <- function(K_req, Xt, y) {
    range_val <- range(Xt)
    if (basis_type == "cubic_bspline") {
      basis <- fda::create.bspline.basis(range_val, nbasis = K_req, norder = 4)
      B_mat <- fda::getbasismatrix(Xt, basis, nderiv = 0)
    } else {
      p_val <- if(is.null(period)) diff(range_val) else period
      basis <- fda::create.fourier.basis(range_val, nbasis = K_req + 1, period = p_val, dropind = 1)
      B_mat <- fda::getbasismatrix(Xt, basis, nderiv = 0)
    }
    lapply(y, function(.x) B_mat)
  }

  # initialization logic
  init_fn_empirical <- function(K_val, m_val) {
    # fitting a linear model to estimate noise variance roughly for starting value
    K_dense <- 50
    if (K_dense >= length(Xt)) K_dense <- length(Xt) - (if(length(Xt) %% 2 == 0) 2 else 1)

    B_est <- build_B_default(K_dense, Xt, y_proc)[[1]]

    sigma2_ests <- vapply(1:m_val, function(i) {
      fit_lm <- tryCatch(lm(y_proc[[i]] ~ B_est - 1), error = function(e) NULL)
      if(is.null(fit_lm)) return(var(y_proc[[i]]))
      summary(fit_lm)$sigma^2
    }, numeric(1))

    mse_avg <- mean(sigma2_ests, na.rm = TRUE)
    if(is.na(mse_avg) || mse_avg <= 0) mse_avg <- 1.0 # Fallback

    # delta2 scale
    delta2_val <- (((length(Xt) + K_val + 2 * delta_1) / 2) - 1) * mse_avg

    list(p = rep(1, m_val * K_val), delta2 = delta2_val, lambda2 = 100, w = 10)
  }

  fn_to_use <- if (!is.null(initial_values_fn)) initial_values_fn else init_fn_empirical

  tuning_info <- NULL
  final_model <- NULL
  selected_K <- NULL
  is_composite <- FALSE

  if (length(K_grid) == 1) {
    # SINGLE K
    k_val <- K_grid
    B_list <- build_B_default(k_val, Xt, y_proc)
    K_effective <- ncol(B_list[[1]])

    final_model <- vem_smooth(
      y = y_proc, B = B_list, Xt = Xt, m = m, K = K_effective,
      initial_values = fn_to_use(K_effective, m),
      lambda_1 = lambda_1, lambda_2 = lambda_2,
      delta_1 = delta_1, delta_2 = delta_2, ...
    )
    final_model$K <- K_effective
    selected_K <- rep(K_effective, m)

  } else {
    # MULTIPLE K
    tuning_info <- tune_vem_by_gcv(
      y = y_proc, Xt = Xt, K_grid = K_grid,
      build_B = build_B_default, initial_values_fn = fn_to_use,
      threshold = threshold,
      mode = selection_metric,
      lambda_1 = lambda_1, lambda_2 = lambda_2,
      delta_1 = delta_1, delta_2 = delta_2, ...
    )

    if (selection_metric == "mean") {
      # GLOBAL MEAN
      k_opt <- tuning_info$best_K_mean
      idx_best <- match(k_opt, K_grid)
      final_model <- tuning_info$fits[[idx_best]]
      selected_K <- rep(k_opt, m)

    } else {
      # PER CURVE
      is_composite <- TRUE
      selected_K <- tuning_info$best_K_per_curve

      # composite model list
      final_model <- vector("list", m)
      names(final_model) <- paste0("Curve_", 1:m)

      for(i in 1:m) {
        k_opt <- selected_K[i]
        idx_fit <- match(k_opt, K_grid)

        # retrieve model fit at optimal K
        full_fit <- tuning_info$fits[[idx_fit]]

        idx_i <- ((i - 1) * k_opt + 1):(i * k_opt)

        final_model[[i]] <- list(
          mu_beta = full_fit$mu_beta[idx_i],
          Sigma_beta = full_fit$Sigma_beta[,,i], # Sigma is array [K,K,m]
          prob = full_fit$prob[idx_i],
          B = full_fit$B[[i]],
          K = k_opt,
          w = full_fit$w,
          delta2 = full_fit$delta2,
          lambda2 = full_fit$lambda2
        )
      }
    }
  }

  # return object
  structure(
    list(
      model = final_model,
      selected_K = selected_K,
      best_K = if(length(unique(selected_K))==1) selected_K[1] else selected_K,
      tuning = tuning_info,
      data_orig = y,
      scaling_params = scaling_params,
      basis_type = basis_type,
      is_composite = is_composite,
      Xt = Xt,
      call = match.call()
    ),
    class = "vem_fit"
  )
}
