if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
p <- 40
q <- 4
n <- 1000
loadings_structure <- build_constr_mat(p, q, "simple")
latent_correlations <- matrix(.3, q, q)
diag(latent_correlations) <- 1
items <- sim_data(
  n,
  gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
  c(-1.2, 0, 1.2),
  latent_correlations,
  SEED = 2
)
data <- set_latpl_data(items)
model <- set_latpl_model(loadings_structure)
fit_pml <- fit_latpl(data, model, CONTROL = set_ucminf_control())
fit_sa <- fit_latpl(data, model, CONTROL = set_sa_control())
epoch_length <- fit_sa$stoFit$epoch_length

#### Helpers ####
distance_in_se <- function(fit_sa, fit_pml) {
  se <- sqrt(compute_var(fit_pml)$asymptotic_variance / fit_pml$dims$n)
  difference <- abs(
    get_par(fit_sa, OPTION = "transformed") -
      get_par(fit_pml, OPTION = "transformed")
  )
  max((difference / se)[se > 0])
}

#### SA estimates are within one standard error of PML ####
local({
  fit_marker_pml <- fit_latpl(
    data,
    set_latpl_model(loadings_structure, STDLV = FALSE),
    CONTROL = set_ucminf_control()
  )
  fit_marker_sa <- fit_latpl(
    data,
    set_latpl_model(loadings_structure, STDLV = FALSE),
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
  fit_diagonal <- fit_latpl(
    data,
    model,
    CONTROL = set_sa_control(GAIN = "diagonal")
  )
  expect_identical(fit_diagonal$stoFit$gain, "diagonal")
  expect_identical(fit_diagonal$stoFit$newton_epochs, 0L)
  expect_true(distance_in_se(fit_diagonal, fit_pml) < 1)
})

#### SA is reproducible with a fixed seed and changes with another seed ####
local({
  expect_identical(
    fit_latpl(data, model, CONTROL = set_sa_control())$theta,
    fit_sa$theta
  )
  expect_false(identical(
    fit_latpl(data, model, CONTROL = set_sa_control(SEED = 7))$theta,
    fit_sa$theta
  ))
})

#### The averaged phase lasts EPOCHS epochs and averages from AVG_FROM of it ####
local({
  expect_equal(epoch_length, ceiling(fit_sa$dims$pairs / 16))
  expect_equal(
    fit_sa$stoFit$averaged,
    fit_sa$control$EPOCHS * epoch_length - floor(0.25 * fit_sa$control$EPOCHS * epoch_length)
  )
  expect_equal(
    fit_sa$stoFit$iterations,
    (fit_sa$stoFit$adam_epochs + fit_sa$stoFit$newton_epochs + fit_sa$control$EPOCHS) *
      epoch_length
  )
  expect_equal(tail(fit_sa$stoFit$path$phase, 1), 3)
})

#### The optimisation noise is the pair-selection covariance over the averaged iterations ####
local({
  vars <- compute_var(fit_sa)
  pairs <- fit_sa$dims$pairs
  pairs_per_iteration <- fit_sa$control$PAIRS_PER_ITERATION
  hj <- lavaan.pl:::cpp_HJ(
    fit_sa$theta,
    fit_sa$engine$model,
    fit_sa$engine$data,
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
  expect_null(compute_var(fit_pml)$noise)
})

#### The optimisation noise matches the spread of SA estimates across seeds in two groups ####
local({
  data_two_groups <- set_latpl_data(items, GROUP = rep(1:2, each = n / 2))
  model_two_groups <- set_latpl_groups(list(model, model))
  optimum <- fit_latpl(data_two_groups, model_two_groups, CONTROL = set_ucminf_control())$theta
  fits <- lapply(1:10, function(seed) {
    fit_latpl(
      data_two_groups,
      model_two_groups,
      CONTROL = set_sa_control(SEED = seed)
    )
  })
  deviations <- t(vapply(
    fits,
    function(fit) fit$theta - optimum,
    numeric(length(optimum))
  ))
  trace_ratio <- sum(colMeans(deviations^2)) /
    sum(diag(compute_var(fits[[1]])$noise))
  expect_true(trace_ratio > 0.6 && trace_ratio < 1.5)
})

#### Linear constraints hold exactly ####
local({
  names_free <- theta_names(data, model)
  constraints <- rbind(
    replace(
      numeric(length(names_free)),
      match(c("lambda[2,1]", "lambda[3,1]"), names_free),
      c(1, -1)
    ),
    replace(numeric(length(names_free)), match("tau[4,2]", names_free), 1)
  )
  model_constrained <- set_latpl_model(
    loadings_structure,
    CONSTRAINTS = list(A = constraints, b = c(0, 0.25))
  )
  fit_constrained <- fit_latpl(
    data,
    model_constrained,
    CONTROL = set_sa_control()
  )
  expect_equal(as.vector(constraints %*% fit_constrained$theta), c(0, 0.25))
  expect_true(
    distance_in_se(fit_constrained, fit_latpl(data, model_constrained, CONTROL = set_ucminf_control())) < 1
  )
})

#### Direct covariances, structural models and the theta parameterization are within one standard error ####
local({
  regressions <- matrix(0, q, q)
  regressions[2, 1] <- regressions[3, 2] <- NA
  residual_covariances <- matrix(0, p, p)
  residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
  for (other_model in list(
    set_latpl_model(loadings_structure, LATCOV = "direct"),
    set_latpl_model(
      loadings_structure,
      CORRFLAG = FALSE,
      REGMAT = regressions,
      RESCOV = residual_covariances
    ),
    set_latpl_model(loadings_structure, PARAMETERIZATION = "theta")
  )) {
    expect_true(
      distance_in_se(
        fit_latpl(data, other_model, CONTROL = set_sa_control()),
        fit_latpl(data, other_model, CONTROL = set_ucminf_control())
      ) <
        1
    )
  }
})

#### Two groups are within one standard error ####
local({
  data_two_groups <- set_latpl_data(items, GROUP = rep(1:2, each = n / 2))
  model_two_groups <- set_latpl_groups(list(
    set_latpl_model(loadings_structure),
    set_latpl_model(loadings_structure, PARAMETERIZATION = "theta")
  ))
  fit_two_groups <- fit_latpl(
    data_two_groups,
    model_two_groups,
    CONTROL = set_sa_control()
  )
  expect_true(fit_two_groups$stoFit$converged)
  expect_true(
    distance_in_se(
      fit_two_groups,
      fit_latpl(data_two_groups, model_two_groups, CONTROL = set_ucminf_control())
    ) <
      1
  )
})

#### Two groups with validation data are within one standard error ####
local({
  validation <- sim_data(
    400,
    gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
    c(-1.2, 0, 1.2),
    latent_correlations,
    SEED = 9
  )
  data_two_groups <- set_latpl_data(
    items,
    GROUP = rep(1:2, each = n / 2),
    VALDATA = validation,
    VALGROUP = rep(1:2, each = 200)
  )
  model_two_groups <- set_latpl_groups(list(
    set_latpl_model(loadings_structure),
    set_latpl_model(loadings_structure)
  ))
  fit_two_groups <- fit_latpl(
    data_two_groups,
    model_two_groups,
    CONTROL = set_sa_control()
  )
  expect_true(fit_two_groups$stoFit$converged)
  expect_true(
    distance_in_se(
      fit_two_groups,
      fit_latpl(data_two_groups, model_two_groups, CONTROL = set_ucminf_control())
    ) <
      1
  )
})

#### Starting values outside the parameter space are rejected ####
local({
  expect_error(
    fit_latpl(
      data,
      model,
      CONTROL = set_sa_control(),
      INIT = set_init(THETA = replace(fit_pml$theta, 121, 5))
    ),
    "outside the parameter space"
  )
})

#### The Adam phase initialises ucminf ####
local({
  expect_equal(
    fit_latpl(data, model, INIT = set_init(METHOD = "SA"), CONTROL = set_ucminf_control())$theta,
    fit_pml$theta,
    tolerance = 1e-4
  )
})

#### From several starts, SA continues from the best and counts the epochs of all ####
local({
  negative_pl <- function(theta) {
    -lavaan.pl:::cpp_pl(theta, fit_sa$engine$model, fit_sa$engine$data)$pl
  }
  starts <- list(fit_sa$init, fit_pml$theta)
  run_sa <- function(starts, control, adam_only) {
    lavaan.pl:::run_sa(
      starts,
      data,
      fit_sa$engine,
      fit_sa$lincon,
      control,
      adam_only,
      FALSE
    )
  }
  one_start <- run_sa(starts[1], set_sa_init_control(), TRUE)
  two_starts <- run_sa(starts, set_sa_init_control(), TRUE)
  full_run <- run_sa(starts, set_sa_control(), FALSE)
  expect_true(negative_pl(two_starts$theta) <= negative_pl(one_start$theta))
  expect_equal(sum(two_starts$path$phase == 0), 2)
  expect_equal(
    full_run$iterations,
    (full_run$adam_epochs + full_run$newton_epochs + set_sa_control()$EPOCHS) * epoch_length
  )
})

#### cfa() with SA adds the optimisation noise to the covariance matrix ####
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
  fit <- cfa(syntax, items_ordered, std.lv = TRUE, control = set_sa_control())
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
  without_noise <- as_lavaan(
    fit@external$lavaan.pl,
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
  data_missing <- set_latpl_data(items_missing, MISSING = "available.cases")
  expect_true(
    distance_in_se(
      fit_latpl(data_missing, model, CONTROL = set_sa_control()),
      fit_latpl(data_missing, model, CONTROL = set_ucminf_control())
    ) <
      1
  )
})

#### SA is the default of fit_latpl(), cfa() and sem() ####
local({
  expect_identical(fit_latpl(data, model)$method, "SA")
  mimic <- as.data.frame(lapply(as.data.frame(items[, 1:6]), ordered))
  names(mimic) <- paste0("y", 1:6)
  for (fit in list(
    lavaan.pl::cfa("f =~ y1 + y2 + y3 + y4 + y5 + y6", mimic),
    lavaan.pl::sem("f =~ y1 + y2 + y3 + y4 + y5 + y6", mimic)
  )) {
    expect_identical(fit@external$lavaan.pl$method, "SA")
    expect_identical(fit@Options$optim.method, "SA")
    expect_true(lavaan::lavInspect(fit, "converged"))
  }
})

#### The end-of-run check stores the variance, predicts the distance and invalidates on a new estimate ####
local({
  vars <- compute_var(fit_sa)
  expect_true(fit_sa$stoFit$check$done)
  expect_true(fit_sa$stoFit$check$pd)
  expect_identical(vars, fit_sa$vars)
  expect_equal(max(vars$distance), fit_sa$stoFit$check$distance)
  moved_vars <- compute_var(replace(fit_sa, "theta", list(fit_sa$theta + 0.02)))
  gradient <- numDeriv::grad(
    function(x) lavaan.pl:::cpp_pl(x, fit_sa$engine$model, fit_sa$engine$data)$pl,
    fit_sa$theta + 0.02
  )
  expect_equal(moved_vars$step, as.vector(solve(moved_vars$H, gradient)) / n, tolerance = 1e-5)
  se_pml <- sqrt(compute_var(fit_pml)$asymptotic_variance / n)
  true_distance <- abs(get_par(fit_sa, OPTION = "transformed") - get_par(fit_pml, OPTION = "transformed")) / se_pml
  expect_true(abs(max(vars$distance) - max(true_distance[se_pml > 0])) < .05)
  moved <- replace(fit_sa, "theta", list(fit_sa$theta + 1e-3))
  expect_false(identical(compute_var(moved)$H, vars$H))
  expect_true(max(compute_var(fit_pml)$distance) < 1e-2)
})

#### The noise vanishes when every iteration uses all the pairs ####
local({
  small <- set_latpl_data(items[, 1:5])
  fit_all <- fit_latpl(small, set_latpl_model(build_constr_mat(5, 1, "simple")))
  expect_equal(fit_all$stoFit$epoch_length, 1)
  expect_equal(compute_var(fit_all)$optimisation_noise, rep(0, length(fit_all$theta)))
})

#### The Adam phase has a budget in iterations when epochs are short ####
local({
  small <- set_latpl_data(items[, 1:6])
  model_small <- set_latpl_model(build_constr_mat(6, 1, "simple"))
  start <- suppressWarnings(fit_latpl(small, model_small, CONTROL = set_ucminf_control(MAXEVAL = 1)))$init
  start[lavaan.pl:::theta_block(model_small, sum(small$cat - 1)) == "lambda"] <- .01
  fit_short <- fit_latpl(
    small,
    model_small,
    CONTROL = set_sa_control(ADAM_EPOCHS = 2, TOL = 1e-8, BURN_EPOCHS = 0),
    INIT = set_init(THETA = start)
  )
  expect_true(fit_short$stoFit$adam_epochs > 2)
  expect_true(fit_short$stoFit$adam_epochs <= 150)
})

#### REFINE takes the average to the exact estimate and keeps it ####
local({
  fit_refined <- fit_latpl(data, model, CONTROL = set_sa_control(REFINE = TRUE))
  expect_equal(fit_refined$theta, fit_pml$theta, tolerance = 1e-3)
  expect_equal(fit_refined$stoFit$theta_average, fit_sa$theta)
  expect_true(fit_refined$stoFit$check$refine_steps > 0)
  expect_true(max(compute_var(fit_refined)$distance) < .01)
  expect_equal(compute_var(fit_refined)$optimisation_noise, rep(0, length(fit_sa$theta)))
  expect_true(fit_refined$stoFit$converged)
})

#### A saddle point is reported, and cfa() returns without standard errors ####
local({
  small <- set_latpl_data(items[, 1:6])
  model_small <- set_latpl_model(build_constr_mat(6, 1, "simple"))
  start <- suppressWarnings(fit_latpl(small, model_small, CONTROL = set_ucminf_control(MAXEVAL = 1)))$init
  start[lavaan.pl:::theta_block(model_small, sum(small$cat - 1)) == "lambda"] <- 0
  expect_warning(
    fit_saddle <- fit_latpl(small, model_small, INIT = set_init(THETA = start)),
    "not positive definite"
  )
  expect_false(fit_saddle$stoFit$converged)
  expect_false(fit_saddle$stoFit$check$pd)
  expect_null(fit_saddle$vars)
  expect_error(compute_var(fit_saddle), "not positive definite")
  frame <- as.data.frame(lapply(as.data.frame(items[, 1:6]), ordered))
  names(frame) <- paste0("y", 1:6)
  syntax <- "f =~ y1 + y2 + y3 + y4 + y5 + y6"
  template <- lavaan::cfa(syntax, frame, estimator = "PML", std.lv = TRUE, do.fit = FALSE)
  zero_loadings <- replace(lavaan::coef(template), grepl("=~", names(lavaan::coef(template))), 0)
  expect_warning(
    fit_lavaan <- lavaan.pl::cfa(syntax, frame, std.lv = TRUE, start = zero_loadings, test = "mean.var.adjusted"),
    "not positive definite"
  )
  expect_false(lavaan::lavInspect(fit_lavaan, "converged"))
  expect_true(all(is.na(lavaan::parameterEstimates(fit_lavaan)$se) | lavaan::parameterEstimates(fit_lavaan)$se == 0))
})

#### Default starts have the sign of the data and respect the marker ####
local({
  reversed <- items[, 1:8]
  reversed[, c(2, 5)] <- max(reversed) - reversed[, c(2, 5)]
  data_reversed <- set_latpl_data(reversed)
  structure_one <- build_constr_mat(8, 1, "simple")
  block <- lavaan.pl:::theta_block(set_latpl_model(structure_one), sum(data_reversed$cat - 1))
  start_std <- fit_latpl(data_reversed, set_latpl_model(structure_one), CONTROL = set_ucminf_control())$init
  expect_equal(sign(start_std[block == "lambda"]), c(1, -1, 1, 1, -1, 1, 1, 1))
  model_marker <- set_latpl_model(structure_one, STDLV = FALSE)
  block_marker <- lavaan.pl:::theta_block(model_marker, sum(data_reversed$cat - 1))
  start_marker <- fit_latpl(data_reversed, model_marker, CONTROL = set_ucminf_control())$init
  expect_equal(sign(start_marker[block_marker == "lambda"]), c(-1, 1, 1, -1, 1, 1, 1))
  expect_true(start_marker[block_marker == "logsd"] > log(.25))
  fit_exact <- fit_latpl(data_reversed, set_latpl_model(structure_one), CONTROL = set_ucminf_control())
  fit_stochastic <- fit_latpl(data_reversed, set_latpl_model(structure_one))
  expect_true(fit_stochastic$stoFit$converged)
  expect_equal(sign(fit_stochastic$theta[block == "lambda"]), sign(fit_exact$theta[block == "lambda"]))
})

#### A distance explained by the noise warns and converges, one it does not explain fails ####
local({
  structure_many <- build_constr_mat(p, q, "simple")
  misfit <- items
  set.seed(14)
  misfit[, 1:10] <- pmin(pmax(misfit[, 1:10] + sample(-1:1, 10 * n, TRUE, c(.15, .7, .15)), 0), 3)
  misfit[, 21:30] <- pmin(pmax(misfit[, 21:30] + misfit[, 1:10] - 1, 0), 3)
  data_misfit <- set_latpl_data(misfit)
  fits <- lapply(1:6, function(seed) {
    suppressWarnings(fit_latpl(
      data_misfit,
      set_latpl_model(structure_many),
      CONTROL = set_sa_control(SEED = seed, EPOCHS = 2)
    ))
  })
  for (fit in fits) {
    expect_identical(
      fit$stoFit$converged,
      fit$stoFit$plateau && fit$stoFit$skipped_averaging == 0 && isTRUE(fit$stoFit$check$pd) &&
        fit$stoFit$check$distance <= max(1, 4 * fit$stoFit$check$noise)
    )
  }
  noisy <- fits[[which.max(vapply(fits, function(x) if (x$stoFit$converged) x$stoFit$check$distance else 0, numeric(1)))]]
  expect_true(noisy$stoFit$check$distance > 1)
  expect_true(noisy$stoFit$converged)
  expect_warning(
    fit_latpl(data_misfit, set_latpl_model(structure_many), CONTROL = noisy$control),
    "account for"
  )
  small <- set_latpl_data(items[, 1:8])
  model_small <- set_latpl_model(build_constr_mat(8, 1, "simple"))
  far <- fit_latpl(small, model_small, CONTROL = set_ucminf_control())$theta
  thresholds <- lavaan.pl:::theta_block(model_small, sum(small$cat - 1)) == "tau"
  far[thresholds] <- far[thresholds] + .3
  expect_warning(
    fit_far <- fit_latpl(
      small,
      model_small,
      CONTROL = set_sa_control(ALPHA = 1e-4, ADAM_EPOCHS = 1, BURN_EPOCHS = 0, EPOCHS = 0),
      INIT = set_init(THETA = far)
    ),
    "more than the optimisation noise explains"
  )
  expect_false(fit_far$stoFit$converged)
})
