# data-raw/generate_toy_curves.R
# Run this script once to generate the toy_curves package dataset.
# Output is saved to data/toy_curves.rda via usethis::use_data().
#
# To regenerate:
#   source("data-raw/generate_toy_curves.R")

library(fda)
library(MASS)

set.seed(1234)

# --- Setup ---
m     <- 3      # number of curves
n     <- 50     # observations per curve
K     <- 8      # basis functions used to generate data
Xt    <- seq(0, 1, length.out = n)
sigma <- 0.1    # noise standard deviation
w     <- 6      # OU correlation decay (normalised scale)

# True coefficients: positions 2, 5 are zero (inactive bases)
true_coef <- c(1.5, 0, -1, 0.8, 0, -0.5, 1.2, -0.9)

# --- Basis matrix ---
basis <- fda::create.bspline.basis(range(Xt), nbasis = K, norder = 4)
B_mat <- fda::getbasismatrix(Xt, basis, nderiv = 0)

# --- OU covariance ---
range_scale <- max(Xt) - min(Xt)
dist_mat    <- abs(outer(Xt, Xt, "-"))
Cov_i       <- (sigma^2) * exp(-w * dist_mat / range_scale)

# --- Generate curves ---
set.seed(1234)
errors <- replicate(m, MASS::mvrnorm(1, mu = rep(0, n), Sigma = Cov_i))
y_list <- lapply(seq_len(m), function(i) {
  as.numeric(B_mat %*% true_coef) + errors[, i]
})
names(y_list) <- paste0("curve_", seq_len(m))

# --- Bundle ---
toy_curves <- list(
  y         = y_list,
  Xt        = Xt,
  true_coef = true_coef,
  K         = K,
  m         = m,
  sigma     = sigma,
  w         = w
)

usethis::use_data(toy_curves, overwrite = TRUE)
