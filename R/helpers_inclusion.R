#' Theta-related Expectation Helpers
#'
#' @description
#' Computes \eqn{E[\log\theta]} and \eqn{E[\log(1-\theta)]} for a
#' \eqn{\text{Beta}(a_1, a_2)} distribution. These terms appear in the
#' variational updates for the inclusion indicators \eqn{Z_{ki}}.
#'
#' @param a1_ki shape1 parameter of the Beta distribution.
#' @param a2_ki shape2 parameter of the Beta distribution.
#'
#' @return Numeric expected values (scalar).
#' @keywords internal

expectedLogTheta <- function(a1_ki, a2_ki) {
  digamma(a1_ki) - digamma(a1_ki+a2_ki)
}

expectedLogOneMinusTheta <- function(a1_ki, a2_ki) {
  digamma(a2_ki) - digamma(a1_ki+a2_ki)
}
