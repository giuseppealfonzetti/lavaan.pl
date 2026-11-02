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
fit_pml <- fit_latpl(data, model)
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
    set_latpl_model(loadings_structure, STDLV = FALSE)
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
    20 * epoch_length - floor(0.25 * 20 * epoch_length)
  )
  expect_equal(
    fit_sa$stoFit$iterations,
    (fit_sa$stoFit$adam_epochs + fit_sa$stoFit$newton_epochs + 20) *
      epoch_length
  )
  expect_equal(tail(fit_sa$stoFit$path$phase, 1), 3)
})

#### The optimisation noise uses the number of averaged iterations ####
local({
  vars <- compute_var(fit_sa)
  pairs <- fit_sa$dims$pairs
  pairs_per_iteration <- fit_sa$control$PAIRS_PER_ITERATION
  hessian_weight <- pairs *
    (pairs - pairs_per_iteration) /
    (pairs_per_iteration * (pairs - 1))
  variability_weight <- (pairs - pairs_per_iteration) /
    (pairs_per_iteration * (pairs - 1))
  inverse_hessian <- solve(vars$H)
  noise <- inverse_hessian %*%
    (hessian_weight * vars$H - variability_weight * vars$J) %*%
    inverse_hessian /
    (n * fit_sa$stoFit$averaged)
  eigen_noise <- eigen((noise + t(noise)) / 2, symmetric = TRUE)
  expect_equal(
    vars$noise,
    eigen_noise$vectors %*%
      (pmax(eigen_noise$values, 0) * t(eigen_noise$vectors))
  )
  expect_equal(
    vars$optimisation_noise,
    diag(vars$trJacob %*% vars$noise %*% t(vars$trJacob))
  )
  expect_null(compute_var(fit_pml)$noise)
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
    distance_in_se(fit_constrained, fit_latpl(data, model_constrained)) < 1
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
        fit_latpl(data, other_model)
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
      fit_latpl(data_two_groups, model_two_groups)
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
    fit_latpl(data, model, INIT = set_init(METHOD = "SA"))$theta,
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
    (full_run$adam_epochs + full_run$newton_epochs + 20) * epoch_length
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
      fit_latpl(data_missing, model)
    ) <
      1
  )
})
