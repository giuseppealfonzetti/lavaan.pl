source("fixtures/simulate.R")
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
data <- lavaan.pl:::lavaanpl_data(items)
model <- lavaan.pl:::lavaanpl_model(loadings_structure)
fit <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control())
frame <- as.data.frame(lapply(as.data.frame(items), ordered))
names(frame) <- paste0("y", seq_len(p))
syntax <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6\nf2 =~ y7 + y8 + y9 + y10 + y11 + y12"
fit_lavaan <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "SA")
fit_exact <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_qn_control())
fit_saqn <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE)
fit_qn <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "QN")

#### sa_diagnostics() reports the run, per parameter and by block ####
local({
  d <- sa_diagnostics(fit)
  expect_inherits(d, "lavaanpl_sa_diagnostics")
  expect_identical(d$converged, fit$stoFit$converged)
  expect_equal(nrow(d$parameters), fit$dims$npar)
  expect_equal(d$parameters$name, lavaan.pl:::parameter_names(data, model)[fit$lincon$free_idx])
  expect_equal(d$parameters$se, sqrt(lavaan.pl:::compute_var(fit)$asymptotic_variance / n)[fit$lincon$free_idx])
  expect_equal(d$parameters$ratio, d$parameters$noise / d$parameters$se)
  expect_equal(max(d$parameters$distance), fit$stoFit$check$distance)
  expect_equal(sum(d$blocks$parameters), fit$dims$npar)
  expect_equal(d$settings$epoch_length, fit$stoFit$epoch_length)
  expect_equal(d$phases$iterations[3], fit$control$average_epochs * fit$stoFit$epoch_length)
  expect_equal(d$epochs_needed, ceiling(fit$control$average_epochs * (max(d$parameters$ratio) / .1)^2))
  expect_stdout(print(d), "converged")
  expect_stdout(print(d), "Optimisation noise")
})

#### Diagnostics take lavaan objects, and reject fits of the exact estimate ####
local({
  d <- sa_diagnostics(fit_lavaan)
  expect_equal(nrow(d$parameters), fit_lavaan@external$lavaanpl_internal$dims$npar)
  expect_error(sa_diagnostics(fit_exact), "not fitted by stochastic approximation")
  expect_error(plot(fit_qn, "noise"), "not fitted by stochastic approximation")
  expect_error(sa_stability(fit_exact), "not fitted by stochastic approximation")
  expect_error(sa_diagnostics(1), "must come from")
})

#### A saddle point still gives a report ####
local({
  start <- fit_exact$init
  start[lavaan.pl:::parameter_blocks(model, sum(data$cat - 1)) == "lambda"] <- 0
  saddle <- suppressWarnings(lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control(), START = start))
  d <- sa_diagnostics(saddle)
  expect_false(d$converged)
  expect_true(all(is.na(d$parameters$ratio)))
  expect_stdout(print(d), "NOT converged")
})

#### plot() returns the data of each view, with every phase ####
local({
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  internal <- fit_lavaan@external$lavaanpl_internal
  trace <- plot(fit_lavaan)
  expect_equal(nrow(trace), length(internal$parameter_vector) * length(internal$stoFit$path$iter))
  expect_equal(unique(trace$epoch), internal$stoFit$path$iter / internal$stoFit$epoch_length)
  expect_equal(levels(droplevels(trace$phase)), c("Adam", "Preconditioning", "Averaged"))
  expect_equal(unique(plot(fit_lavaan, "trace", blocks = "lambda")$block), "lambda")
  expect_true(all(is.finite(plot(fit_lavaan, "objective")$objective)))
  expect_equal(nrow(plot(fit_lavaan, "noise")), internal$dims$npar)

  saqn <- fit_saqn@external$lavaanpl_internal
  path <- saqn$numFit$path
  expect_equal(tail(path$parameter_vector, 1)[[1]], saqn$parameter_vector)
  expect_true(all(diff(path$nll) < 0))
  trace <- plot(fit_saqn)
  expect_equal(levels(droplevels(trace$phase)), c("Adam", "Preconditioning", "QN"))
  expect_true(all(diff(unique(trace$epoch)) > 0))
  expect_true(all(is.finite(plot(fit_saqn, "objective")$objective)))
  expect_error(plot(fit_saqn, "noise"), "not fitted by stochastic approximation")

  qn <- fit_qn@external$lavaanpl_internal
  expect_equal(qn$numFit$path$parameter_vector[[1]], qn$init)
  expect_equal(tail(qn$numFit$path$parameter_vector, 1)[[1]], qn$parameter_vector)
  trace <- plot(fit_qn)
  expect_equal(levels(droplevels(trace$phase)), "QN")
  expect_equal(min(trace$epoch), 0)
  expect_error(plot(fit_qn, "other"))
})

#### sa_stability() compares the seeds with the predicted noise ####
local({
  stability <- sa_stability(fit_lavaan, seeds = 1:6)
  expect_equal(ncol(stability$estimates), 6)
  expect_equal(stability$converged, 6)
  expect_true(all(stability$parameters$observed > 0))
  ratio <- stability$parameters$predicted / stability$parameters$observed
  expect_true(stats::median(ratio) > .5 && stats::median(ratio) < 2)
  expect_stdout(print(stability), "6 seeds")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(plot(stability), stability$parameters)
  expect_error(sa_stability(fit, seeds = 1))
})
