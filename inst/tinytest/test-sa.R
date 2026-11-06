source("fixtures/simulate.R")
if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
p <- 40
q <- 4
n <- 1000
loadings_structure <- build_constr_mat(p, q, "simple")
loadings_markers <- replace(loadings_structure, cbind(seq(1, p, p / q), 1:q), 1)
latent_correlations <- matrix(.3, q, q)
diag(latent_correlations) <- 1
items <- sim_data(
  n,
  gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
  c(-1.2, 0, 1.2),
  latent_correlations,
  SEED = 2
)
data <- lavaan.pl:::lavaanpl_data(items)
model <- lavaan.pl:::lavaanpl_model(loadings_structure)
fit_pml <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_qn_control())
fit_sa <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control())
epoch_length <- fit_sa$stoFit$epoch_length

#### Helpers ####
distance_in_se <- function(fit_sa, fit_pml) {
  se <- sqrt(lavaan.pl:::compute_var(fit_pml)$asymptotic_variance / fit_pml$data$n)
  difference <- abs(
    lavaan.pl:::transform_parameters(fit_sa$model, fit_sa$parameter_vector) -
      lavaan.pl:::transform_parameters(fit_pml$model, fit_pml$parameter_vector)
  )
  max((difference / se)[se > 0])
}

#### SA estimates are within one standard error of PML ####
local({
  fit_marker_pml <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA, q, q)),
    CONTROL = set_qn_control()
  )
  fit_marker_sa <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA, q, q)),
    CONTROL = set_sa_control()
  )
  for (fits in list(
    list(fit_sa, fit_pml),
    list(fit_marker_sa, fit_marker_pml)
  )) {
    expect_true(fits[[1]]$stoFit$converged)
    expect_identical(fits[[1]]$stoFit$gain, "hessian")
    expect_true(distance_in_se(fits[[1]], fits[[2]]) < 1)
  }
})

#### The diagonal gain skips the Newton phase ####
local({
  fit_diagonal <- lavaan.pl:::fit_lavaanpl(
    data,
    model,
    CONTROL = set_sa_control(gain = "diagonal")
  )
  expect_identical(fit_diagonal$stoFit$gain, "diagonal")
  expect_identical(fit_diagonal$stoFit$precondition_epochs, 0L)
  expect_true(distance_in_se(fit_diagonal, fit_pml) < 1)
})

#### SA is reproducible with a fixed seed and changes with another seed ####
local({
  expect_identical(
    lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control())$parameter_vector,
    fit_sa$parameter_vector
  )
  expect_false(identical(
    lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control(seed = 7))$parameter_vector,
    fit_sa$parameter_vector
  ))
})

#### The averaged phase lasts average_epochs epochs and averages from avg_from of it ####
local({
  expect_equal(epoch_length, ceiling(fit_sa$data$pairs / 16))
  expect_equal(
    fit_sa$stoFit$averaged,
    fit_sa$control$average_epochs * epoch_length - floor(0.25 * fit_sa$control$average_epochs * epoch_length)
  )
  expect_equal(
    fit_sa$stoFit$iterations,
    (fit_sa$stoFit$adam_epochs + fit_sa$stoFit$precondition_epochs + fit_sa$control$average_epochs) *
      epoch_length
  )
  expect_equal(tail(fit_sa$stoFit$path$phase, 1), 3)
})

#### The optimisation noise is the pair-selection covariance over the averaged iterations ####
local({
  vars <- lavaan.pl:::compute_var(fit_sa)
  pairs <- fit_sa$data$pairs
  pairs_per_iteration <- fit_sa$control$pairs_per_iteration
  hj <- lavaan.pl:::cpp_sandwich(
    fit_sa$parameter_vector,
    fit_sa$cpp_inputs$model,
    fit_sa$cpp_inputs$data,
    TRUE
  )
  pair_selection <- (pairs - pairs_per_iteration) /
    (pairs_per_iteration * (pairs - 1)) *
    (pairs * hj$pairvar / n^2 - tcrossprod(hj$gradient / n))
  inverse_hessian <- solve(vars$H)
  expect_equal(
    vars$noise,
    inverse_hessian %*%
      pair_selection %*%
      inverse_hessian /
      fit_sa$stoFit$averaged
  )
  expect_true(all(
    eigen(vars$noise, symmetric = TRUE, only.values = TRUE)$values > -1e-12
  ))
  expect_equal(
    vars$optimisation_noise,
    diag(vars$trJacob %*% vars$noise %*% t(vars$trJacob))
  )
  expect_null(lavaan.pl:::compute_var(fit_pml)$noise)
})

#### The optimisation noise matches the spread of SA estimates across seeds in two groups ####
local({
  data_two_groups <- lavaan.pl:::lavaanpl_data(items, GROUP = rep(1:2, each = n / 2))
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(model, model))
  optimum <- lavaan.pl:::fit_lavaanpl(data_two_groups, model_two_groups, CONTROL = set_qn_control())$parameter_vector
  fits <- lapply(1:10, function(seed) {
    lavaan.pl:::fit_lavaanpl(
      data_two_groups,
      model_two_groups,
      CONTROL = set_sa_control(seed = seed)
    )
  })
  deviations <- t(vapply(
    fits,
    function(fit) fit$parameter_vector - optimum,
    numeric(length(optimum))
  ))
  trace_ratio <- sum(colMeans(deviations^2)) /
    sum(diag(lavaan.pl:::compute_var(fits[[1]])$noise))
  expect_true(trace_ratio > 0.6 && trace_ratio < 1.5)
})

#### Linear constraints hold exactly ####
local({
  names_free <- lavaan.pl:::parameter_names(data, model)
  constraints <- rbind(
    replace(
      numeric(length(names_free)),
      match(c("lambda[2,1]", "lambda[3,1]"), names_free),
      c(1, -1)
    ),
    replace(numeric(length(names_free)), match("tau[4,2]", names_free), 1)
  )
  model_constrained <- lavaan.pl:::lavaanpl_model(
    loadings_structure,
    CONSTRAINTS = list(A = constraints, b = c(0, 0.25))
  )
  fit_constrained <- lavaan.pl:::fit_lavaanpl(
    data,
    model_constrained,
    CONTROL = set_sa_control()
  )
  expect_equal(as.vector(constraints %*% fit_constrained$parameter_vector), c(0, 0.25))
  expect_true(
    distance_in_se(fit_constrained, lavaan.pl:::fit_lavaanpl(data, model_constrained, CONTROL = set_qn_control())) < 1
  )
})

#### Direct covariances, structural models and the theta parameterization are within one standard error ####
local({
  regressions <- matrix(0, q, q)
  regressions[2, 1] <- regressions[3, 2] <- NA
  residual_covariances <- matrix(0, p, p)
  residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
  for (other_model in list(
    lavaan.pl:::lavaanpl_model(loadings_structure, LATCOV = "direct"),
    lavaan.pl:::lavaanpl_model(
      loadings_structure,
      PSI = diag(q),
      BETA = regressions,
      THETA = residual_covariances
    ),
    lavaan.pl:::lavaanpl_model(loadings_structure, PARAMETERIZATION = "theta")
  )) {
    expect_true(
      distance_in_se(
        lavaan.pl:::fit_lavaanpl(data, other_model, CONTROL = set_sa_control()),
        lavaan.pl:::fit_lavaanpl(data, other_model, CONTROL = set_qn_control())
      ) <
        1
    )
  }
})

#### Two groups are within one standard error ####
local({
  data_two_groups <- lavaan.pl:::lavaanpl_data(items, GROUP = rep(1:2, each = n / 2))
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(
    lavaan.pl:::lavaanpl_model(loadings_structure),
    lavaan.pl:::lavaanpl_model(loadings_structure, PARAMETERIZATION = "theta")
  ))
  fit_two_groups <- lavaan.pl:::fit_lavaanpl(
    data_two_groups,
    model_two_groups,
    CONTROL = set_sa_control()
  )
  expect_true(fit_two_groups$stoFit$converged)
  expect_true(
    distance_in_se(
      fit_two_groups,
      lavaan.pl:::fit_lavaanpl(data_two_groups, model_two_groups, CONTROL = set_qn_control())
    ) <
      1
  )
})

#### Two groups are within one standard error ####
local({
  data_two_groups <- lavaan.pl:::lavaanpl_data(items, GROUP = rep(1:2, each = n / 2))
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(
    lavaan.pl:::lavaanpl_model(loadings_structure),
    lavaan.pl:::lavaanpl_model(loadings_structure)
  ))
  fit_two_groups <- lavaan.pl:::fit_lavaanpl(
    data_two_groups,
    model_two_groups,
    CONTROL = set_sa_control()
  )
  expect_true(fit_two_groups$stoFit$converged)
  expect_true(
    distance_in_se(
      fit_two_groups,
      lavaan.pl:::fit_lavaanpl(data_two_groups, model_two_groups, CONTROL = set_qn_control())
    ) <
      1
  )
})

#### Starting values outside the parameter space are rejected ####
local({
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      data,
      model,
      CONTROL = set_sa_control(),
      START = replace(fit_pml$parameter_vector, 121, 5)
    ),
    "outside the parameter space"
  )
})

#### pmlcfa() with SA adds the optimisation noise to the covariance matrix ####
local({
  colnames(items) <- paste0("y", 1:p)
  items_ordered <- as.data.frame(lapply(as.data.frame(items), ordered))
  syntax <- paste(
    sapply(1:q, function(h) {
      paste0(
        "f",
        h,
        " =~ ",
        paste0("y", which(is.na(loadings_structure[, h])), collapse = " + ")
      )
    }),
    collapse = "\n"
  )
  fit <- pmlcfa(syntax, items_ordered, std.lv = TRUE, control = set_sa_control())
  template <- lavaan::cfa(
    syntax,
    items_ordered,
    std.lv = TRUE,
    estimator = "PML",
    se = "robust.huber.white",
    information = "observed",
    test = "none",
    do.fit = FALSE
  )
  without_noise <- lavaan.pl:::create_lav_from_lavaanpl(
    fit@external$lavaanpl_internal,
    modifyList(fit@external$vars, list(noise = NULL)),
    template
  )
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_true(all(is.finite(lavaan::parameterEstimates(fit)$se)))
  expect_true(all(diag(lavaan::vcov(fit)) >= diag(lavaan::vcov(without_noise))))
  expect_true(any(diag(lavaan::vcov(fit)) > diag(lavaan::vcov(without_noise))))
  expect_equal(
    attr(fit@vcov$vcov, "E.inv"),
    attr(without_noise@vcov$vcov, "E.inv")
  )
})

#### SA handles available cases ####
local({
  items_missing <- unname(items)
  set.seed(6)
  for (h in 1:q) {
    indicators <- which(is.na(loadings_structure[, h]))
    for (j in indicators[-1]) {
      items_missing[
        runif(n) < c(.1, .25, .25, .5)[items[, indicators[1]] + 1],
        j
      ] <- NA
    }
  }
  data_missing <- lavaan.pl:::lavaanpl_data(items_missing, MISSING = "available.cases")
  expect_true(
    distance_in_se(
      lavaan.pl:::fit_lavaanpl(data_missing, model, CONTROL = set_sa_control()),
      lavaan.pl:::fit_lavaanpl(data_missing, model, CONTROL = set_qn_control())
    ) <
      1
  )
})

#### SA-QN is the default of lavaan.pl:::fit_lavaanpl(), pmlcfa() and pmlsem() ####
local({
  expect_identical(lavaan.pl:::fit_lavaanpl(data, model)$method, "SA-QN")
  mimic <- as.data.frame(lapply(as.data.frame(items[, 1:6]), ordered))
  names(mimic) <- paste0("y", 1:6)
  for (fit in list(
    lavaan.pl::pmlcfa("f =~ y1 + y2 + y3 + y4 + y5 + y6", mimic),
    lavaan.pl::pmlsem("f =~ y1 + y2 + y3 + y4 + y5 + y6", mimic)
  )) {
    expect_identical(fit@external$lavaanpl_internal$method, "SA-QN")
    expect_identical(fit@Options$optim.method, "SA-QN")
    expect_true(lavaan::lavInspect(fit, "converged"))
  }
})

#### The end-of-run check stores the variance, predicts the distance and invalidates on a new estimate ####
local({
  vars <- lavaan.pl:::compute_var(fit_sa)
  expect_true(fit_sa$stoFit$check$done)
  expect_true(fit_sa$stoFit$check$pd)
  expect_identical(vars, fit_sa$vars)
  expect_equal(max(vars$distance), fit_sa$stoFit$check$distance)
  moved_vars <- lavaan.pl:::compute_var(replace(fit_sa, "parameter_vector", list(fit_sa$parameter_vector + 0.02)))
  gradient <- lavaan.pl:::cpp_pl(fit_sa$parameter_vector + 0.02, fit_sa$cpp_inputs$model, fit_sa$cpp_inputs$data, 1L)$gradient
  expect_equal(moved_vars$step, as.vector(solve(moved_vars$H, gradient)) / n, tolerance = 1e-8)
  se_pml <- sqrt(lavaan.pl:::compute_var(fit_pml)$asymptotic_variance / n)
  true_distance <- abs(lavaan.pl:::transform_parameters(fit_sa$model, fit_sa$parameter_vector) - lavaan.pl:::transform_parameters(fit_pml$model, fit_pml$parameter_vector)) / se_pml
  expect_true(abs(max(vars$distance) - max(true_distance[se_pml > 0])) < .05)
  moved <- replace(fit_sa, "parameter_vector", list(fit_sa$parameter_vector + 1e-3))
  expect_false(identical(lavaan.pl:::compute_var(moved)$H, vars$H))
  expect_true(max(lavaan.pl:::compute_var(fit_pml)$distance) < 1e-2)
})

#### The noise vanishes when every iteration uses all the pairs ####
local({
  small <- lavaan.pl:::lavaanpl_data(items[, 1:5])
  fit_all <- lavaan.pl:::fit_lavaanpl(small, lavaan.pl:::lavaanpl_model(build_constr_mat(5, 1, "simple")), CONTROL = set_sa_control())
  expect_equal(fit_all$stoFit$epoch_length, 1)
  expect_equal(lavaan.pl:::compute_var(fit_all)$optimisation_noise, rep(0, length(fit_all$parameter_vector)))
})

#### The maximum epochs of each phase are respected, and the default Adam budget is in iterations ####
local({
  small <- lavaan.pl:::lavaanpl_data(items[, 1:6])
  model_small <- lavaan.pl:::lavaanpl_model(build_constr_mat(6, 1, "simple"))
  start <- suppressWarnings(lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_qn_control(maxeval = 1)))$init
  start[lavaan.pl:::parameter_blocks(model_small, sum(small$cat - 1)) == "lambda"] <- .01
  fit_default <- lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_sa_control(tol = 1e-8), START = start)
  expect_equal(fit_default$stoFit$pairs_per_iteration, 15)
  expect_equal(fit_default$stoFit$epoch_length, 1)
  expect_equal(fit_default$stoFit$adam_epochs, 150)
  fit_capped <- suppressWarnings(lavaan.pl:::fit_lavaanpl(
    small,
    model_small,
    CONTROL = set_sa_control(adam_epochs = 2, precondition_epochs = 1, max_precondition_epochs = 1, average_epochs = 2, tol = 1e-8),
    START = start
  ))
  expect_equal(fit_capped$stoFit$adam_epochs, 2)
  expect_equal(fit_capped$stoFit$precondition_epochs, 1)
  expect_equal(sum(fit_capped$stoFit$path$phase == 3), 2)
  fit_large <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control(tol = 1e-8))
  expect_equal(fit_large$stoFit$adam_epochs, 5)
  fit_none <- suppressWarnings(lavaan.pl:::fit_lavaanpl(
    data,
    model,
    CONTROL = set_sa_control(precondition_epochs = 0)
  ))
  expect_equal(fit_none$stoFit$precondition_epochs, 0)
  expect_identical(fit_none$stoFit$gain, "hessian")
  expect_true(fit_none$stoFit$converged)
  flat <- lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_qn_control())$parameter_vector
  flat[lavaan.pl:::parameter_blocks(model_small, sum(small$cat - 1)) == "lambda"] <- 0
  still <- function(cap) suppressWarnings(lavaan.pl:::fit_lavaanpl(
    small,
    model_small,
    CONTROL = set_sa_control(
      alpha = 1e-4, adam_epochs = 1, precondition_epochs = 0, max_precondition_epochs = cap, average_epochs = 1
    ),
    START = flat
  ))$stoFit$precondition_epochs
  expect_equal(still(0), 0)
  expect_true(still(20) > 0)
  expect_true(still(20) <= 20)
})

#### A saddle point is reported, and pmlcfa() returns without standard errors ####
local({
  small <- lavaan.pl:::lavaanpl_data(items[, 1:6])
  model_small <- lavaan.pl:::lavaanpl_model(build_constr_mat(6, 1, "simple"))
  start <- suppressWarnings(lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_qn_control(maxeval = 1)))$init
  start[lavaan.pl:::parameter_blocks(model_small, sum(small$cat - 1)) == "lambda"] <- 0
  expect_warning(
    fit_saddle <- lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_sa_control(), START = start),
    "not positive definite"
  )
  expect_false(fit_saddle$stoFit$converged)
  expect_false(fit_saddle$stoFit$check$pd)
  expect_null(fit_saddle$vars)
  expect_error(lavaan.pl:::compute_var(fit_saddle), "not positive definite")
  frame <- as.data.frame(lapply(as.data.frame(items[, 1:6]), ordered))
  names(frame) <- paste0("y", 1:6)
  syntax <- "f =~ y1 + y2 + y3 + y4 + y5 + y6"
  template <- lavaan::cfa(syntax, frame, estimator = "PML", std.lv = TRUE, do.fit = FALSE)
  zero_loadings <- replace(lavaan::coef(template), grepl("=~", names(lavaan::coef(template))), 0)
  expect_warning(
    fit_lavaan <- lavaan.pl::pmlcfa(syntax, frame, std.lv = TRUE, start = zero_loadings, test = "mean.var.adjusted", optim.method = "SA"),
    "not positive definite"
  )
  expect_false(lavaan::lavInspect(fit_lavaan, "converged"))
  expect_true(all(is.na(lavaan::parameterEstimates(fit_lavaan)$se) | lavaan::parameterEstimates(fit_lavaan)$se == 0))
})

#### Default starts have the sign of the data and respect the marker ####
local({
  reversed <- items[, 1:8]
  reversed[, c(2, 5)] <- max(reversed) - reversed[, c(2, 5)]
  data_reversed <- lavaan.pl:::lavaanpl_data(reversed)
  structure_one <- build_constr_mat(8, 1, "simple")
  block <- lavaan.pl:::parameter_blocks(lavaan.pl:::lavaanpl_model(structure_one), sum(data_reversed$cat - 1))
  start_std <- lavaan.pl:::fit_lavaanpl(data_reversed, lavaan.pl:::lavaanpl_model(structure_one), CONTROL = set_qn_control())$init
  expect_equal(sign(start_std[block == "lambda"]), c(1, -1, 1, 1, -1, 1, 1, 1))
  model_marker <- lavaan.pl:::lavaanpl_model(replace(structure_one, 1, 1), PSI = matrix(NA))
  block_marker <- lavaan.pl:::parameter_blocks(model_marker, sum(data_reversed$cat - 1))
  start_marker <- lavaan.pl:::fit_lavaanpl(data_reversed, model_marker, CONTROL = set_qn_control())$init
  expect_equal(sign(start_marker[block_marker == "lambda"]), c(-1, 1, 1, -1, 1, 1, 1))
  expect_true(start_marker[block_marker == "logsd"] > log(.25))
  fit_exact <- lavaan.pl:::fit_lavaanpl(data_reversed, lavaan.pl:::lavaanpl_model(structure_one), CONTROL = set_qn_control())
  fit_stochastic <- lavaan.pl:::fit_lavaanpl(data_reversed, lavaan.pl:::lavaanpl_model(structure_one), CONTROL = set_sa_control())
  expect_true(fit_stochastic$stoFit$converged)
  expect_equal(sign(fit_stochastic$parameter_vector[block == "lambda"]), sign(fit_exact$parameter_vector[block == "lambda"]))
})

#### A distance explained by the noise warns and converges, one it does not explain fails ####
local({
  structure_many <- build_constr_mat(p, q, "simple")
  misfit <- items
  set.seed(14)
  misfit[, 1:10] <- pmin(pmax(misfit[, 1:10] + sample(-1:1, 10 * n, TRUE, c(.15, .7, .15)), 0), 3)
  misfit[, 21:30] <- pmin(pmax(misfit[, 21:30] + misfit[, 1:10] - 1, 0), 3)
  data_misfit <- lavaan.pl:::lavaanpl_data(misfit)
  fits <- lapply(1:6, function(seed) {
    suppressWarnings(lavaan.pl:::fit_lavaanpl(
      data_misfit,
      lavaan.pl:::lavaanpl_model(structure_many),
      CONTROL = set_sa_control(seed = seed, average_epochs = 2)
    ))
  })
  for (fit in fits) {
    expect_identical(
      fit$stoFit$converged,
      fit$stoFit$skipped_averaging == 0 && isTRUE(fit$stoFit$check$pd) &&
        fit$stoFit$check$distance <= max(1, 4 * fit$stoFit$check$noise)
    )
  }
  noisy <- fits[[which.max(vapply(fits, function(x) if (x$stoFit$converged) x$stoFit$check$distance else 0, numeric(1)))]]
  expect_true(noisy$stoFit$check$distance > 1)
  expect_true(noisy$stoFit$converged)
  expect_warning(
    lavaan.pl:::fit_lavaanpl(data_misfit, lavaan.pl:::lavaanpl_model(structure_many), CONTROL = noisy$control),
    "account for"
  )
  small <- lavaan.pl:::lavaanpl_data(items[, 1:8])
  model_small <- lavaan.pl:::lavaanpl_model(build_constr_mat(8, 1, "simple"))
  far <- lavaan.pl:::fit_lavaanpl(small, model_small, CONTROL = set_qn_control())$parameter_vector
  thresholds <- lavaan.pl:::parameter_blocks(model_small, sum(small$cat - 1)) == "tau"
  far[thresholds] <- far[thresholds] + .3
  expect_warning(
    fit_far <- lavaan.pl:::fit_lavaanpl(
      small,
      model_small,
      CONTROL = set_sa_control(alpha = 1e-4, adam_epochs = 1, precondition_epochs = 0, average_epochs = 1),
      START = far
    ),
    "more than the optimisation noise explains"
  )
  expect_false(fit_far$stoFit$converged)
})

#### SA-QN starts the exact optimiser from the burn-in and its curvature ####
local({
  fit_saqn <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_saqn_control())
  expect_identical(fit_saqn$method, "SA-QN")
  expect_equal(fit_saqn$parameter_vector, fit_pml$parameter_vector, tolerance = 1e-4)
  expect_equal(lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_saqn_control(seed = 7))$parameter_vector, fit_saqn$parameter_vector, tolerance = 1e-4)
  expect_null(lavaan.pl:::compute_var(fit_saqn)$noise)
  expect_null(fit_saqn$stoFit$invh)
  burn <- lavaan.pl:::run_sa(fit_saqn$init, data, fit_saqn$cpp_inputs, fit_saqn$lincon, set_saqn_control(), FALSE)
  expect_equal(dim(burn$invh), rep(length(fit_saqn$parameter_vector), 2))
  fit_flat <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_qn_control(), START = burn$parameter_vector)
  expect_true(fit_saqn$numFit$info[["neval"]] < fit_flat$numFit$info[["neval"]])
  expect_equal(length(lavaan.pl:::run_sa(fit_saqn$init, data, fit_saqn$cpp_inputs, fit_saqn$lincon, set_sa_control(average_epochs = 1), FALSE)$invh), 0)
})

#### Several SA-QN starts keep the highest pairwise log-likelihood of the runs with those seeds ####
local({
  logl <- function(fit) lavaan.pl:::cpp_pl(fit$parameter_vector, fit$cpp_inputs$model, fit$cpp_inputs$data)$pl
  runs <- lapply(5:7, function(seed) lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_saqn_control(seed = seed)))
  fit_three <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_saqn_control(seed = 5, nstarts = 3))
  expect_equal(fit_three$stoFit$start_logl, vapply(runs, logl, numeric(1)))
  expect_equal(logl(fit_three), max(vapply(runs, logl, numeric(1))))
  expect_equal(fit_three$parameter_vector, runs[[which.max(vapply(runs, logl, numeric(1)))]]$parameter_vector)
  expect_identical(length(runs[[1]]$stoFit$start_logl), 1L)
})

#### SA needs the averaged phase, whose noise enters the variance ####
local({
  expect_error(lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_sa_control(average_epochs = 0)), "average_epochs")
})
