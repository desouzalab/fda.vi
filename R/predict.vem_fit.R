#' Predict Method for VEM Fits
#'
#' @description
#' Returns posterior mean curve estimates from a \code{vem_fit} object.
#' Active basis functions are selected by thresholding the posterior inclusion
#' probabilities at 0.5. If \code{newdata} is supplied, a new basis matrix is
#' constructed at those time points; otherwise the original fitted time points
#' are used. Predictions are automatically back-transformed if the model was
#' fitted with \code{center = TRUE} or \code{scale = TRUE}.
#'
#' @param object A \code{vem_fit} object from \code{\link{vem_fit}}.
#' @param newdata Optional numeric vector of new time points at which to
#'   evaluate the fitted curves. Must lie within the original domain
#'   \code{range(Xt)}. If \code{NULL}, predictions are returned at the
#'   original \code{Xt}.
#' @param ... Currently unused.
#'
#' @return A list of length \eqn{m}. Each element is a numeric vector of
#'   predicted values on the original (back-transformed) scale.
#'
#' @references
#' da Cruz, A. C., de Souza, C. P. E., & Sousa, P. H. T. O. (2024).
#' Fast Bayesian basis selection for functional data representation with
#' correlated errors. \emph{arXiv:2405.20758}.
#' \url{https://arxiv.org/abs/2405.20758}
#'
#' @seealso \code{\link{vem_fit}}, \code{\link{plot.vem_fit}},
#'   \code{\link{coef.vem_fit}}
#' @export
#' @examples
#' \donttest{
#'   data(toy_curves)
#'   fit <- vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8)
#'
#'   # Predictions at original time points
#'   preds <- predict(fit)
#'   length(preds)       # 3 — one vector per curve
#'
#'   # Predictions at a denser grid
#'   Xt_new <- seq(0, 1, length.out = 200)
#'   preds_dense <- predict(fit, newdata = Xt_new)
#'
#'   # Plot observed vs predicted for curve 1
#'   plot(toy_curves$Xt, toy_curves$y[[1]],
#'        pch = 16, cex = 0.6, col = "grey50",
#'        xlab = "t", ylab = "y")
#'   lines(Xt_new, preds_dense[[1]], col = "firebrick", lwd = 2)
#' }

predict.vem_fit <- function(object, newdata = NULL, ...) {
  m <- length(object$data_orig)
  preds <- vector("list", m)


  for (i in 1:m) {
    # model params
    if (object$is_composite) {
      # for per curve, params are in list [[i]]
      mod <- object$model[[i]]
      K_i <- mod$K
      mu_i <- mod$mu_beta
      prob_i <- mod$prob

      # basis function setup
      if(is.null(newdata)) {
        B_i <- mod$B
      } else {
        range_val <- range(object$Xt)
        if (object$basis_type == "cubic_bspline") {
          basis <- fda::create.bspline.basis(range_val, nbasis = K_i, norder = 4)
        } else {
          basis <- fda::create.fourier.basis(range_val, nbasis = K_i + 1, dropind = 1)
        }
        B_i <- fda::getbasismatrix(newdata, basis, nderiv = 0)
      }

    } else {
      # GLOBAL logic
      K <- object$best_K
      idx_start <- (i - 1) * K + 1
      idx_end   <- i * K

      mu_i <- object$model$mu_beta[idx_start:idx_end]
      prob_i <- object$model$prob[idx_start:idx_end]

      if(is.null(newdata)) {
        B_i <- object$model$B[[i]]
      } else {
        range_val <- range(object$Xt)
        if (object$basis_type == "cubic_bspline") {
          basis <- fda::create.bspline.basis(range_val, nbasis = K, norder = 4)
        } else {
          basis <- fda::create.fourier.basis(range_val, nbasis = K + 1, dropind = 1)
        }
        B_i <- fda::getbasismatrix(newdata, basis, nderiv = 0)
      }
    }

    #computing active basis
    z_mode <- ifelse(prob_i > 0.5, 1, 0)
    mu_active <- mu_i * z_mode
    y_pred <- as.vector(B_i %*% mu_active)

    #destandardzing the data
    if (!is.null(object$scaling_params)) {
      s_mean <- object$scaling_params$means[[i]]
      s_sd   <- object$scaling_params$sds[[i]]
      y_pred <- (y_pred * s_sd) + s_mean
    }
    preds[[i]] <- y_pred
  }
  return(preds)
}
