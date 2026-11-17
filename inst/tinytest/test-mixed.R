if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
n <- 2000
p <- 6
set.seed(1)
covariates <- cbind(rbinom(n, 1, .5), rnorm(n))
eta <- as.vector(covariates %*% c(.5, -.3)) + rnorm(n)
latent <- sapply(seq_len(p), function(j) .7 * eta + rnorm(n, 0, sqrt(.51)))
numeric_items <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
responses <- latent
responses[, 1] <- 10 * latent[, 1] + 3
responses[, 2] <- 0.1 * latent[, 2] - 1
for (j in 3:p) responses[, j] <- findInterval(latent[, j], c(-.8, 0, .8))
loadings_structure <- matrix(NA, p, 1)
data_mixed <- lavaan.pl:::lavaanpl_data(responses, NUMERIC = numeric_items)
model_mixed <- lavaan.pl:::lavaanpl_model(loadings_structure, NUMERIC = numeric_items)
fit_mixed <- lavaan.pl:::fit_lavaanpl(data_mixed, model_mixed, CONTROL = set_control("QN"))
data_both <- lavaan.pl:::lavaanpl_data(responses, NUMERIC = numeric_items, X = covariates)
model_both <- lavaan.pl:::lavaanpl_model(
  loadings_structure,
  GAMMA = matrix(NA, 1, 2),
  NUMERIC = numeric_items
)
fit_both <- lavaan.pl:::fit_lavaanpl(data_both, model_both, CONTROL = set_control("QN"))
frame <- as.data.frame(responses)
names(frame) <- paste0("y", seq_len(p))
for (j in 3:p) frame[[j]] <- ordered(frame[[j]])
frame$x1 <- covariates[, 1]
frame$x2 <- covariates[, 2]
syntax <- "f =~ y1 + y2 + y3 + y4 + y5 + y6"

#### Numeric items have no thresholds and free intercepts and residual standard deviations ####
local({
  expect_equal(data_mixed$cat[1:2], c(1, 1), check.attributes = FALSE)
  expect_equal(fit_mixed$dims$nthr, 12)
  expect_equal(fit_mixed$dims$nint, 2)
  expect_equal(fit_mixed$dims$nscale, 2)
  expect_equal(
    lavaan.pl:::parameter_names(data_mixed, model_mixed)[c(19, 20, 21, 22)],
    c("nu[1]", "nu[2]", "logrsd[1]", "logrsd[2]")
  )
  expect_equal(length(data_mixed$freq), 6 * 16)
  par <- lavaan.pl:::vector_to_matrices(fit_mixed$model, fit_mixed$parameter_vector)[[1]]
  expect_equal(diag(par$theta)[3:p], rep(0, 4))
  expect_true(all(diag(par$theta)[1:2] > 0))
  expect_equal(par$delta, rep(1, p))
})

#### The parameters are recovered, on the scale of each item ####
local({
  par <- lavaan.pl:::vector_to_matrices(fit_mixed$model, fit_mixed$parameter_vector)[[1]]
  total <- par$lambda^2 + c(diag(par$theta)[1:2], 1, 1, 1, 1)
  expect_equal(par$nu[1:2], colMeans(responses[, 1:2]), tolerance = .02)
  expect_equal(sqrt(total[1:2]), apply(responses[, 1:2], 2, sd), tolerance = .02)
  par_both <- lavaan.pl:::vector_to_matrices(fit_both$model, fit_both$parameter_vector)[[1]]
  se <- sqrt(lavaan.pl:::compute_var(fit_both)$asymptotic_variance / n)
  gamma <- as.vector(par_both$gamma)
  expect_true(all(abs(gamma / tail(se, 2)) > 3))
  expect_equal(sign(gamma), sign(c(.5, -.3)))
})

#### Linear rescaling of a numeric item only rescales its parameters and shifts the likelihood ####
local({
  rescaled <- responses
  rescaled[, 1] <- 100 * rescaled[, 1] + 7
  fit_rescaled <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(rescaled, NUMERIC = numeric_items), model_mixed, CONTROL = set_control("QN"))
  par <- lavaan.pl:::vector_to_matrices(fit_mixed$model, fit_mixed$parameter_vector)[[1]]
  par_rescaled <- lavaan.pl:::vector_to_matrices(fit_rescaled$model, fit_rescaled$parameter_vector)[[1]]
  expect_equal(par_rescaled$lambda[1], 100 * par$lambda[1], tolerance = 1e-3)
  expect_equal(par_rescaled$lambda[-1], par$lambda[-1], tolerance = 1e-3)
  expect_equal(par_rescaled$nu[1], 100 * par$nu[1] + 7, tolerance = 1e-3)
  pl <- function(x) lavaan.pl:::cpp_pl(x$parameter_vector, x$cpp_inputs$model, x$cpp_inputs$data)$pl
  expect_equal(pl(fit_rescaled), pl(fit_mixed) - n * (p - 1) * log(100), tolerance = 1e-6)
})

#### Frequency tables and individual observations agree with numeric items ####
local({
  discrete <- cbind(rbinom(n, 1, .5), rbinom(n, 1, .3))
  data_discrete <- lavaan.pl:::lavaanpl_data(responses, NUMERIC = numeric_items, X = discrete)
  expect_true(data_discrete$tables)
  model_discrete <- lavaan.pl:::lavaanpl_model(
    loadings_structure,
    GAMMA = matrix(NA, 1, 2),
    NUMERIC = numeric_items
  )
  theta <- suppressWarnings(lavaan.pl:::fit_lavaanpl(
    data_discrete,
    model_discrete,
    CONTROL = set_control("QN", maxeval = 1)
  ))$init + 0.02
  with_tables <- lavaan.pl:::cpp_inputs(data_discrete, lavaan.pl:::lavaanpl_groups(list(model_discrete)), TRUE)
  with_units <- lavaan.pl:::cpp_inputs(data_discrete, lavaan.pl:::lavaanpl_groups(list(model_discrete)), FALSE)
  pl_tables <- lavaan.pl:::cpp_pl(theta, with_tables$model, with_tables$data, 1L)
  pl_units <- lavaan.pl:::cpp_pl(theta, with_units$model, with_units$data, 1L)
  hj_tables <- lavaan.pl:::cpp_sandwich(theta, with_tables$model, with_tables$data, TRUE)
  hj_units <- lavaan.pl:::cpp_sandwich(theta, with_units$model, with_units$data, TRUE)
  expect_equal(pl_units$pl, pl_tables$pl, tolerance = 1e-10)
  expect_equal(pl_units$gradient, pl_tables$gradient, tolerance = 1e-8)
  expect_equal(hj_units$H, hj_tables$H, tolerance = 1e-8)
  expect_equal(hj_units$J, hj_tables$J, tolerance = 1e-8)
  expect_equal(hj_units$pairvar, hj_tables$pairvar, tolerance = 1e-8)
})

#### Mixed items with covariates: the objective is lavaan's, at the optimum of its gradient ####
local({
  fit <- pmlsem(paste(syntax, "f ~ x1 + x2", sep = "\n"), frame, std.lv = TRUE, ordered = FALSE, control = set_control("QN"))
  engine_fit <- fit@external$lavaanpl_internal
  lavaan_logl <- function(x) {
    table <- lavaan::parTable(x)
    table$start <- table$est
    lavaan::lavaan(table, frame, estimator = "PML", do.fit = FALSE)@optim$logl
  }
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_equal(fit@optim$logl, lavaan_logl(fit), tolerance = 1e-7)
  gradient <- lavaan.pl:::cpp_pl(
    engine_fit$parameter_vector,
    engine_fit$cpp_inputs$model,
    engine_fit$cpp_inputs$data,
    1L
  )$gradient
  expect_true(max(abs(gradient)) / n < 1e-4)
  fit_direct <- pmlsem(
    paste(syntax, "f ~ x1\ny1 ~ x2\ny3 ~ x2", sep = "\n"),
    frame,
    std.lv = TRUE,
    ordered = FALSE
  )
  engine_direct <- fit_direct@external$lavaanpl_internal
  expect_equal(engine_direct$dims$nkap, 2)
  expect_equal(fit_direct@optim$logl, lavaan_logl(fit_direct), tolerance = 1e-7)
  expect_true(is.finite(coef(fit_direct)["y1~x2"]))
})

#### Tests with numeric items ####
local({
  test <- lavaan.pl:::compute_plrt(fit_both, lavaan.pl:::compute_var(fit_both))
  expect_equal(test$model$df, 12 + 2 * p + 2 * 2 + p * (p - 1) / 2 - fit_both$dims$npar)
  expect_equal(test$baseline$df, p * (p - 1) / 2)
  expect_true(is.finite(test$model$stat.scaled))
  expect_true(is.finite(test$baseline$stat.scaled))
  fit_nested <- lavaan.pl:::fit_lavaanpl(
    data_both,
    lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(c(NA, 0), 1, 2), NUMERIC = numeric_items),
    CONTROL = set_control("QN")
  )
  difference <- lavaan.pl:::compute_plrt(fit_nested, H1 = fit_both)$model
  expect_equal(difference$df, 1)
  expect_true(difference$pvalue.scaled < .05)
})

#### SA agrees with PML on mixed items, with and without covariates ####
local({
  for (fit in list(fit_mixed, fit_both)) {
    fit_sa <- lavaan.pl:::fit_lavaanpl(
      fit$data,
      fit$model,
      CONTROL = set_control("SA", pairs_per_iteration = 4)
    )
    se <- sqrt(lavaan.pl:::compute_var(fit)$asymptotic_variance / n)
    expect_true(fit_sa$stoFit$converged)
    expect_true(max(abs(
      lavaan.pl:::transform_parameters(fit_sa$model, fit_sa$parameter_vector) - lavaan.pl:::transform_parameters(fit$model, fit$parameter_vector)
    ) / se) < 1)
    expect_true(all(lavaan.pl:::compute_var(fit_sa)$optimisation_noise >= 0))
  }
})

#### Threads give the same numbers on every path ####
local({
  serial <- lavaan.pl:::cpp_sandwich(fit_both$parameter_vector, fit_both$cpp_inputs$model, fit_both$cpp_inputs$data, TRUE)
  serial_pl <- lavaan.pl:::cpp_pl(fit_both$parameter_vector, fit_both$cpp_inputs$model, fit_both$cpp_inputs$data, 1L)
  serial_test <- lavaan.pl:::compute_plrt(fit_both, lavaan.pl:::compute_var(fit_both))
  threaded <- fit_both
  threaded$cores <- 4
  RcppParallel::setThreadOptions(numThreads = 4)
  parallel <- lavaan.pl:::cpp_sandwich(fit_both$parameter_vector, fit_both$cpp_inputs$model, fit_both$cpp_inputs$data, TRUE)
  parallel_pl <- lavaan.pl:::cpp_pl(fit_both$parameter_vector, fit_both$cpp_inputs$model, fit_both$cpp_inputs$data, 1L)
  RcppParallel::setThreadOptions(numThreads = 1)
  expect_equal(parallel_pl$pl, serial_pl$pl, tolerance = 1e-12)
  expect_equal(parallel_pl$gradient, serial_pl$gradient, tolerance = 1e-8)
  expect_equal(parallel$H, serial$H, tolerance = 1e-10)
  expect_equal(parallel$J, serial$J, tolerance = 1e-10)
  expect_equal(parallel$pairvar, serial$pairvar, tolerance = 1e-10)
  parallel_test <- lavaan.pl:::compute_plrt(threaded, lavaan.pl:::compute_var(threaded))
  expect_equal(parallel_test$model$stat.scaled, serial_test$model$stat.scaled, tolerance = 1e-8)
})

#### Inputs are validated ####
local({
  expect_error(
    lavaan.pl:::lavaanpl_data(responses, NUMERIC = numeric_items, MISSING = "available.cases"),
    "available.cases"
  )
  expect_error(lavaan.pl:::lavaanpl_data(responses, NUMERIC = c(TRUE, FALSE)))
  expect_error(
    lavaan.pl:::fit_lavaanpl(data_mixed, lavaan.pl:::lavaanpl_model(loadings_structure), CONTROL = set_control("QN")),
    "same numeric items"
  )
  expect_error(pmlsem(paste(syntax, "f ~ x1", sep = "\n"), frame, ordered = FALSE, missing = "available.cases"), "available.cases")
})

#### Observed mediation with numeric variables equals ML when just identified ####
local({
  set.seed(11)
  x <- rnorm(500)
  m <- 0.5 * x + rnorm(500)
  data <- data.frame(x, m, y = 0.7 * m + rnorm(500))
  model <- "y ~ c*x\nm ~ a*x\ny ~ b*m\nab := a*b"
  expect_equal(coef(pmlsem(model, data, ordered = FALSE)), coef(lavaan::sem(model, data, conditional.x = TRUE)), tolerance = 1e-5)
})
