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
frame <- as.data.frame(lapply(as.data.frame(items), ordered))
names(frame) <- paste0("y", seq_len(p))
syntax <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6\nf2 =~ y7 + y8 + y9 + y10 + y11 + y12"
fit_lavaan <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "SA")
fit_saqn <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE)
fit_qn <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "QN")

#### diagnostics() reports the phases and checks of each method ####
local({
  for (x in list(fit_lavaan, fit_saqn, fit_qn)) {
    internal <- x@external$lavaanpl_internal
    d <- diagnostics(x)
    expect_inherits(d, "lavaanpl_diagnostics")
    expect_identical(d$method, internal$method)
    expect_identical(d$converged, lavaan::lavInspect(x, "converged"))
    expect_equal(nrow(d$parameters), internal$dims$npar)
    expect_true(d$checks$ok[d$checks$check == "Observed information positive definite"])
    expect_stdout(print(d), "converged")
  }
  sa <- diagnostics(fit_lavaan)
  expect_equal(sa$phases$phase, c("Adam", "Preconditioning", "Averaged"))
  expect_equal(sa$parameters$ratio, sa$parameters$noise / sa$parameters$se)
  expect_equal(max(sa$parameters$distance), fit_lavaan@external$lavaanpl_internal$stoFit$check$distance)
  expect_true("Skipped steps in the averaged phase" %in% sa$checks$check)

  saqn <- diagnostics(fit_saqn)
  expect_equal(saqn$phases$phase, c("Adam", "Preconditioning", "QN"))
  expect_null(saqn$parameters$noise)

  qn <- diagnostics(fit_qn)
  expect_equal(qn$phases$phase, "QN")
  expect_equal(qn$phases$epochs, fit_qn@external$lavaanpl_internal$numFit$info[["neval"]])
  expect_true(all(qn$checks$ok %in% c(TRUE, NA)))
  expect_true(max(qn$parameters$distance) < 1e-3)
  expect_error(diagnostics(1), "must come from")
})

#### A saddle point and an exhausted quasi-Newton are flagged ####
local({
  start <- lavaan::coef(fit_lavaan)
  start[grepl("=~", names(start))] <- 0
  saddle <- suppressWarnings(lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "SA", start = start))
  d <- diagnostics(saddle)
  expect_false(d$converged)
  expect_true(all(is.na(d$parameters$ratio)))
  expect_false(d$checks$ok[d$checks$check == "Observed information positive definite"])
  expect_stdout(print(d), "NOT converged")

  short <- suppressWarnings(lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, optim.method = "QN", se = "none", test = "none", control = set_qn_control(maxeval = 5)))
  d <- diagnostics(short)
  expect_false(d$converged)
  expect_false(d$checks$ok[d$checks$check == "QN evaluations"])
  expect_stdout(print(d), "! QN evaluations")
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
  objective <- plot(fit_lavaan, "objective")
  expect_true(all(is.finite(objective$objective)))
  expect_equal(nrow(objective), length(internal$stoFit$path$iter))
  expect_true("Averaged" %in% objective$phase)
  expect_equal(tail(objective$objective, 1), -lavaan.pl:::cpp_pl(internal$parameter_vector, internal$cpp_inputs$model, internal$cpp_inputs$data)$pl / internal$data$n)
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
  expect_error(plot(fit_qn, "noise"), "not fitted by stochastic approximation")
  expect_error(plot(fit_qn, "other"))
})
