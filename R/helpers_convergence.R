#' ELBO Convergence Check helper method
#'
#' @description
#' Returns \code{TRUE} when the absolute ELBO change between consecutive
#' iterations falls below \code{convergence_threshold}. Returns \code{FALSE}
#' if either ELBO value is \code{NULL} or \code{NA}.
#'
#' @param elbo_c Numeric scalar. ELBO at the current iteration.
#' @param elbo_prev Numeric scalar. ELBO at the previous iteration.
#' @param convergence_threshold Positive scalar. Tolerance for convergence.
#'
#' @return Logical scalar.
#' @keywords internal


checkConvergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if (is.null(elbo_prev) || is.na(elbo_prev)) return(FALSE)
  if (is.null(elbo_c)    || is.na(elbo_c))    return(FALSE)

  abs(elbo_c - elbo_prev) <= convergence_threshold
}
