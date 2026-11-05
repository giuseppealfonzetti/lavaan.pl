if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
p <- 12
n <- 800
loadings_structure <- build_constr_mat(p, 2, "simple")
latent_correlations <- matrix(.3, 2, 2)
diag(latent_correlations) <- 1
items <- sim_data(
  n,
  gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
  c(-1, 0, 1),
  latent_correlations,
  SEED = 2
)
data <- set_latpl_data(items)
model <- set_latpl_model(loadings_structure)
fit <- fit_latpl(data, model, CONTROL = set_sa_control())
frame <- as.data.frame(lapply(as.data.frame(items), ordered))
names(frame) <- paste0("y", seq_len(p))
syntax <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6\nf2 =~ y7 + y8 + y9 + y10 + y11 + y12"
fit_lavaan <- lavaan.pl::cfa(syntax, frame, std.lv = TRUE, optim.method = "SA")
fit_exact <- fit_latpl(data, model, CONTROL = set_qn_control())

#### sa_diagnostics() reports the run, per parameter and by block ####
local({
  d <- sa_diagnostics(fit)
  expect_inherits(d, "latpl_sa_diagnostics")
  expect_identical(d$converged, fit$stoFit$converged)
  expect_equal(nrow(d$parameters), fit$dims$npar)
  expect_equal(d$parameters$name, theta_names(data, model)[fit$lincon$free_idx])
  expect_equal(d$parameters$se, sqrt(compute_var(fit)$asymptotic_variance / n)[fit$lincon$free_idx])
  expect_equal(d$parameters$ratio, d$parameters$noise / d$parameters$se)
  expect_equal(max(d$parameters$distance), fit$stoFit$check$distance)
  expect_equal(sum(d$blocks$parameters), fit$dims$npar)
  expect_equal(d$settings$epoch_length, fit$stoFit$epoch_length)
  expect_equal(d$phases$iterations[3], fit$control$AVERAGE_EPOCHS * fit$stoFit$epoch_length)
  expect_equal(d$epochs_needed, ceiling(fit$control$AVERAGE_EPOCHS * (max(d$parameters$ratio) / .1)^2))
  expect_stdout(print(d), "converged")
  expect_stdout(print(d), "Optimisation noise")
})

#### Diagnostics take lavaan objects, and reject fits of the exact estimate ####
local({
  d <- sa_diagnostics(fit_lavaan)
  expect_equal(nrow(d$parameters), fit_lavaan@external$lavaan.pl$dims$npar)
  expect_error(sa_diagnostics(fit_exact), "not fitted by stochastic approximation")
  expect_error(plot_sa(fit_exact), "not fitted by stochastic approximation")
  expect_error(sa_stability(fit_exact), "not fitted by stochastic approximation")
  expect_error(sa_diagnostics(1), "must come from")
})

#### A saddle point still gives a report ####
local({
  start <- fit_exact$init
  start[lavaan.pl:::theta_block(model, sum(data$cat - 1)) == "lambda"] <- 0
  saddle <- suppressWarnings(fit_latpl(data, model, CONTROL = set_sa_control(), INIT = set_init(THETA = start)))
  d <- sa_diagnostics(saddle)
  expect_false(d$converged)
  expect_true(all(is.na(d$parameters$ratio)))
  expect_stdout(print(d), "NOT converged")
})

#### plot_sa() returns the data of each view ####
local({
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  trajectory <- plot_sa(fit, "trajectory")
  expect_equal(nrow(trajectory), length(fit$theta) * length(fit$stoFit$path$iter))
  expect_equal(unique(trajectory$iter), fit$stoFit$path$iter)
  expect_equal(unique(plot_sa(fit, "trajectory", BLOCKS = "lambda")$block), "lambda")
  objective <- plot_sa(fit, "objective")
  expect_true(all(is.finite(objective$objective)))
  noise <- plot_sa(fit_lavaan, "noise")
  expect_equal(nrow(noise), fit_lavaan@external$lavaan.pl$dims$npar)
  expect_error(plot_sa(fit, "other"))
})

#### sa_stability() compares the seeds with the predicted noise ####
local({
  stability <- sa_stability(fit_lavaan, SEEDS = 1:6)
  expect_equal(ncol(stability$estimates), 6)
  expect_equal(stability$converged, 6)
  expect_true(all(stability$parameters$observed > 0))
  ratio <- stability$parameters$predicted / stability$parameters$observed
  expect_true(stats::median(ratio) > .5 && stats::median(ratio) < 2)
  expect_stdout(print(stability), "6 seeds")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(plot(stability), stability$parameters)
  expect_error(sa_stability(fit, SEEDS = 1))
})
