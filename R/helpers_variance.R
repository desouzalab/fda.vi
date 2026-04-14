#' Expected Values Under Inverse-Gamma Variational Variance Parameters (sigma and tau)
#'
#' @description
#' Scalar expectations used in the VEM update equations for the noise variance
#' \eqn{\sigma^2 \sim \text{IG}(\delta_1^*, \delta_2^*)} and coefficient scale
#' \eqn{\tau^2 \sim \text{IG}(\lambda_1^*, \lambda_2^*)}.
#'
#' @param delta1,delta2 Shape and scale of \eqn{q(\sigma^2)}.
#' @param lambda1,lambda2 Shape and scale of \eqn{q(\tau^2)}.
#'
#' @return A numeric scalar.
#' @keywords internal
#' @name variance_expectations
NULL

expectedInvTau2 <- function(lambda1, lambda2) lambda1/lambda2
expectedInvSigma2 <- function(delta1, delta2) delta1/delta2
expectedLogSigma2 <- function(delta1, delta2) log(delta2) - digamma(delta1)
expectedLogTau2 <- function(lambda1, lambda2) log(lambda2) - digamma(lambda1)
