# tests for vem_fit and S3 methods

data(toy_curves)
set.seed(1234)
fit_base <- vem_fit(
  y      = toy_curves$y,
  Xt     = toy_curves$Xt,
  K      = 8,
  center = FALSE,
  scale  = FALSE
)

# basic output structure

test_that("vem_fit returns correct structure on toy_curves", {
  expect_s3_class(fit_base, "vem_fit")
  expect_equal(fit_base$best_K,       8L)
  expect_equal(fit_base$is_composite, FALSE)
  expect_equal(fit_base$basis_type,   "cubic_bspline")
  expect_true(fit_base$model$converged)
  expect_true(is.finite(tail(fit_base$model$elbo_values, 1)))
  expect_length(fit_base$data_orig, 3L)
})

# parameter estimation

test_that("vem_fit recovers inactive bases and plausible parameters", {
  K <- fit_base$best_K
  m <- length(toy_curves$y)
  pip_mat <- matrix(fit_base$model$prob, nrow = K, ncol = m)

  # on average across curves, PIPs for inactive bases (2 and 5) should be low
  expect_true(mean(pip_mat[2, ]) < 0.5)
  expect_true(mean(pip_mat[5, ]) < 0.5)

  # PIPs in [0, 1]
  expect_true(all(fit_base$model$prob >= 0 & fit_base$model$prob <= 1))

  # decay parameter w should be positive
  expect_true(fit_base$model$w > 0)

  # posterior IG parameters should be strictly positive
  expect_true(fit_base$model$delta1 > 0)
  expect_true(fit_base$model$delta2 > 0)
  expect_true(fit_base$model$lambda1 > 0)
  expect_true(fit_base$model$lambda2 > 0)
})

# coef()

test_that("coef returns correctly shaped matrix", {
  coefs <- coef(fit_base)
  expect_true(is.matrix(coefs))
  expect_equal(dim(coefs), c(8L, 3L))
  expect_equal(rownames(coefs), paste0("B", 1:8))
  expect_equal(colnames(coefs), paste0("Curve_", 1:3))
})

# predict()

test_that("predict returns list of correct length and dimension", {
  preds <- predict(fit_base)
  expect_type(preds, "list")
  expect_length(preds, 3L)
  expect_length(preds[[1]], length(toy_curves$Xt))

  # predictions at a new dense grid
  Xt_new <- seq(0, 1, length.out = 200)
  preds_new <- predict(fit_base, newdata = Xt_new)
  expect_length(preds_new, 3L)
  expect_length(preds_new[[1]], 200L)
  expect_true(all(sapply(preds_new, is.finite)))
})

# summary()

test_that("summary runs without error and returns active_bases", {
  s <- summary(fit_base)
  expect_named(s, "active_bases")
  expect_length(s$active_bases, 3L)
  expect_true(all(s$active_bases >= 0 & s$active_bases <= 8))
})

# plot()

test_that("plot runs without error for each curve", {
  for (i in 1:3) {
    expect_silent(plot(fit_base, curve_idx = i))
  }
  expect_silent(plot(fit_base, curve_idx = 1, show_CI   = FALSE))
  expect_silent(plot(fit_base, curve_idx = 1, show_basis = TRUE))
})

# GCV automatic K selection

test_that("vem_fit selects K via GCV when given a vector", {
  set.seed(42)
  fit_gcv <- vem_fit(
    y    = toy_curves$y,
    Xt   = toy_curves$Xt,
    K    = c(6, 8, 10)
  )
  expect_s3_class(fit_gcv, "vem_fit")
  expect_true(fit_gcv$best_K %in% c(6L, 8L, 10L))
  expect_false(is.null(fit_gcv$tuning))
  expect_equal(dim(fit_gcv$tuning$gcv_matrix), c(3L, 3L))
})

# per-curve composite K selection

test_that("vem_fit per_curve produces a composite fit", {
  set.seed(42)
  fit_pc <- vem_fit(
    y                = toy_curves$y,
    Xt               = toy_curves$Xt,
    K                = c(6, 8, 10),
    selection_metric = "per_curve"
  )
  expect_true(fit_pc$is_composite)
  expect_length(fit_pc$selected_K, 3L)
  expect_true(all(fit_pc$selected_K %in% c(6L, 8L, 10L)))

  # coef matrix has max(K) rows
  coefs <- coef(fit_pc)
  expect_equal(nrow(coefs), max(fit_pc$selected_K))
  expect_equal(ncol(coefs), 3L)

  # summary should not error on composite fit
  expect_no_error(summary(fit_pc))
})

# Fourier basis

test_that("vem_fit works with Fourier basis", {
  set.seed(42)
  fit_f <- vem_fit(
    y          = toy_curves$y,
    Xt         = toy_curves$Xt,
    K          = 8,
    basis_type = "fourier"
  )
  expect_s3_class(fit_f, "vem_fit")
  expect_equal(fit_f$basis_type, "fourier")
  expect_true(fit_f$model$converged)

  preds <- predict(fit_f)
  expect_length(preds, 3L)
  expect_true(all(sapply(preds, function(p) all(is.finite(p)))))
})

# center and scale back-transformation

test_that("predictions are back-transformed when center/scale are used", {
  set.seed(42)
  fit_scaled <- vem_fit(
    y      = toy_curves$y,
    Xt     = toy_curves$Xt,
    K      = 8,
    center = TRUE,
    scale  = TRUE
  )
  preds_scaled <- predict(fit_scaled)

  # predictions should be on the original scale — compare range to raw data
  raw_range  <- range(unlist(toy_curves$y))
  pred_range <- range(unlist(preds_scaled))
  expect_true(pred_range[1] >= raw_range[1] - 1)
  expect_true(pred_range[2] <= raw_range[2] + 1)
  expect_true(all(sapply(preds_scaled, function(p) all(is.finite(p)))))
})

# input validation

test_that("vem_fit errors on invalid input", {
  expect_error(
    vem_fit(y = matrix(rnorm(50), 10, 5), Xt = toy_curves$Xt, K = 8),
    "must be a list"
  )
  expect_error(
    vem_fit(y = list(), Xt = toy_curves$Xt, K = 8),
    "must be a list"
  )
  expect_error(
    vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8,
            basis_type = "invalid_basis")
  )
  expect_error(
    vem_fit(y = toy_curves$y, Xt = toy_curves$Xt, K = 8,
            selection_metric = "invalid_metric")
  )
})
