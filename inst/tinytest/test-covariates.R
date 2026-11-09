if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
simulate_mimic <- function(n, p, covariates, gamma, seed) {
  set.seed(seed)
  eta <- as.vector(covariates %*% gamma) + rnorm(n)
  items <- sapply(seq_len(p), function(j) {
    findInterval(.7 * eta + rnorm(n, 0, sqrt(.51)), c(-.8, 0, .8))
  })
  storage.mode(items) <- "double"
  items
}
n <- 1500
p <- 6
set.seed(1)
binary <- cbind(rbinom(n, 1, .5), rbinom(n, 1, .3))
items <- simulate_mimic(n, p, binary, c(.5, -.3), 2)
loadings_structure <- matrix(NA, p, 1)
loadings_markers <- replace(loadings_structure, 1, 1)
data <- lavaan.pl:::lavaanpl_data(items, X = binary)
model <- lavaan.pl:::lavaanpl_model(
  loadings_markers,
  PSI = matrix(NA),
  GAMMA = matrix(NA, 1, 2)
)
fit <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_control("QN"))

#### Covariates define the strata of the frequency tables ####
local({
  expect_equal(nrow(data$strata[[1]]), 4)
  expect_equal(length(data$freq), 4 * length(lavaan.pl:::cpp_pairs_freq(items, data$cat)))
  expect_equal(sum(data$freq), n * p * (p - 1) / 2)
  expect_equal(
    data$strata[[1]][data$stratum + 1, , drop = FALSE],
    binary,
    check.attributes = FALSE
  )
  first <- data$stratum == 0
  expect_equal(
    data$freq[seq_along(lavaan.pl:::cpp_pairs_freq(items, data$cat))],
    lavaan.pl:::cpp_pairs_freq(items[first, ], data$cat)
  )
  expect_equal(ncol(fit$data$X), 2)
  expect_equal(fit$dims$ngam, 2)
  expect_equal(
    tail(lavaan.pl:::parameter_names(data, model), 2),
    c("gamma[1,1]", "gamma[1,2]")
  )
})

#### A model with discrete covariates is a constrained multigroup model ####
local({
  key <- paste(binary[, 1], binary[, 2])
  patterns <- c("0 0", "1 0", "0 1", "1 1")
  ordering <- order(match(key, patterns))
  data_groups <- lavaan.pl:::lavaanpl_data(items[ordering, ], GROUP = key[ordering])
  groups <- lavaan.pl:::lavaanpl_groups(list(
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), ALPHA = 0),
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), ALPHA = NA),
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), ALPHA = NA),
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), ALPHA = NA)
  ))
  names_groups <- lavaan.pl:::parameter_names(data_groups, groups)
  shared <- names_groups[!grepl("\\.g|^alpha", names_groups)]
  constraint <- function(names, values) {
    replace(numeric(length(names_groups)), match(names, names_groups), values)
  }
  equalities <- do.call(rbind, lapply(2:4, function(g) {
    t(sapply(shared, function(x) {
      constraint(c(x, paste0(x, ".g", g)), c(1, -1))
    }))
  }))
  additive <- constraint(
    c("alpha[1].g4", "alpha[1].g2", "alpha[1].g3"),
    c(1, -1, -1)
  )
  groups$CONSTRAINTS <- list(
    A = unname(rbind(equalities, additive)),
    b = numeric(nrow(equalities) + 1)
  )
  fit_groups <- lavaan.pl:::fit_lavaanpl(data_groups, groups, CONTROL = set_control("QN"))
  par_groups <- lavaan.pl:::vector_to_matrices(groups, fit_groups$parameter_vector)
  par <- lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]
  expect_equal(par$tau, par_groups[[1]]$tau, tolerance = 1e-4)
  expect_equal(par$lambda, par_groups[[1]]$lambda, tolerance = 1e-4)
  expect_equal(
    par$psi,
    par_groups[[1]]$psi,
    tolerance = 1e-4
  )
  expect_equal(
    as.vector(par$gamma),
    c(par_groups[[2]]$alpha, par_groups[[3]]$alpha),
    tolerance = 1e-4
  )
  expect_equal(
    lavaan.pl:::cpp_pl(fit$parameter_vector, fit$cpp_inputs$model, fit$cpp_inputs$data)$pl,
    lavaan.pl:::cpp_pl(fit_groups$parameter_vector, fit_groups$cpp_inputs$model, fit_groups$cpp_inputs$data)$pl,
    tolerance = 1e-8
  )
})

#### Covariates in the data with no effect in the model leave the fit unchanged ####
local({
  plain <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items), lavaan.pl:::lavaanpl_model(loadings_structure), CONTROL = set_control("QN"))
  ignored <- lavaan.pl:::fit_lavaanpl(data, lavaan.pl:::lavaanpl_model(loadings_structure), CONTROL = set_control("QN"))
  expect_equal(ignored$parameter_vector, plain$parameter_vector, tolerance = 1e-5)
})

#### The effects of the covariates are recovered ####
local({
  gamma <- as.vector(lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]$gamma)
  se <- sqrt(lavaan.pl:::compute_var(fit)$asymptotic_variance / n)
  expect_true(all(abs(tail(lavaan.pl:::transform_parameters(fit$model, fit$parameter_vector), 2) - gamma) < 1e-12))
  expect_true(all(gamma / tail(se, 2) * sign(c(.5, -.3)) > 3))
})

#### Direct effects shift a single item ####
local({
  shifted <- items
  shifted[, 1] <- pmin(pmax(shifted[, 1] + binary[, 2], 0), 3)
  kappa <- matrix(0, p, 2)
  kappa[1, 2] <- NA
  fit_direct <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(shifted, X = binary),
    lavaan.pl:::lavaanpl_model(
      loadings_markers,
      PSI = matrix(NA),
      GAMMA = matrix(NA, 1, 2),
      KAPPA = kappa
    ),
    CONTROL = set_control("QN")
  )
  par <- lavaan.pl:::vector_to_matrices(fit_direct$model, fit_direct$parameter_vector)[[1]]
  expect_true(par$kappa[1, 2] > .5)
  expect_equal(sum(par$kappa != 0), 1)
})

#### Observations with missing covariates are removed ####
local({
  with_missing <- binary
  with_missing[1:20, 1] <- NA
  expect_message(
    data_missing <- lavaan.pl:::lavaanpl_data(items, X = with_missing),
    "20 observations"
  )
  expect_equal(data_missing$n, n - 20)
  expect_equal(
    data_missing$freq,
    lavaan.pl:::lavaanpl_data(items[-(1:20), ], X = binary[-(1:20), ])$freq
  )
})

#### Frequency weights equal replicated rows ####
local({
  times <- rep(1:3, length.out = n)
  rows <- rep(seq_len(n), times)
  fit_replicated <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items[rows, ], X = binary[rows, ]),
    model,
    CONTROL = set_control("QN")
  )
  fit_weighted <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(
      items,
      X = binary,
      WEIGHTS = times,
      WEIGHTS_TYPE = "frequency"
    ),
    model,
    CONTROL = set_control("QN")
  )
  expect_equal(fit_weighted$parameter_vector, fit_replicated$parameter_vector, tolerance = 1e-5)
})

#### Inputs are validated ####
local({
  expect_error(lavaan.pl:::lavaanpl_data(items, X = binary[1:10, ]))
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      lavaan.pl:::lavaanpl_data(items),
      lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(NA, 1, 2)),
      CONTROL = set_control("QN")
    ),
    "one column per covariate"
  )
  expect_error(lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(NA, 2, 2)))
})

#### SA with covariates agrees with PML ####
local({
  se <- sqrt(lavaan.pl:::compute_var(fit)$asymptotic_variance / n)
  distance <- function(x) {
    max(abs(lavaan.pl:::transform_parameters(x$model, x$parameter_vector) - lavaan.pl:::transform_parameters(fit$model, fit$parameter_vector)) / se)
  }
  fit_sa <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_control("SA", pairs_per_iteration = 4))
  expect_true(fit_sa$stoFit$converged)
  expect_true(distance(fit_sa) < 1)
  vars_sa <- lavaan.pl:::compute_var(fit_sa)
  expect_true(all(vars_sa$optimisation_noise >= 0))
})

#### Frequency tables and individual observations give the same likelihood, gradient, H, J and pair variance ####
local({
  simulated_missing <- items
  set.seed(9)
  simulated_missing[sample(length(items), 300)] <- NA
  gamma <- matrix(NA, 1, 2)
  kappa <- matrix(0, p, 2)
  kappa[2, 1] <- NA
  cases <- list(
    list(data = lavaan.pl:::lavaanpl_data(items), model = lavaan.pl:::lavaanpl_model(loadings_structure)),
    list(data = data, model = lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = gamma, KAPPA = kappa)),
    list(
      data = lavaan.pl:::lavaanpl_data(simulated_missing, MISSING = "available.cases", X = binary),
      model = lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), GAMMA = gamma)
    ),
    list(
      data = lavaan.pl:::lavaanpl_data(
        items,
        X = binary,
        WEIGHTS = rep(1:3, length.out = n),
        CLUSTER = rep(seq_len(n / 10), each = 10)
      ),
      model = lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = gamma)
    ),
    list(
      data = lavaan.pl:::lavaanpl_data(items, X = binary, GROUP = rep(c("a", "b"), length.out = n)),
      model = lavaan.pl:::lavaanpl_groups(list(
        lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = gamma),
        lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), KAPPA = kappa)
      ))
    )
  )
  for (case in cases) {
    start <- suppressWarnings(lavaan.pl:::fit_lavaanpl(
      case$data,
      case$model,
      CONTROL = set_control("QN", maxeval = 1)
    ))
    theta <- start$init + 0.05
    with_tables <- lavaan.pl:::cpp_inputs(case$data, start$model, TRUE)
    with_units <- lavaan.pl:::cpp_inputs(case$data, start$model, FALSE)
    pl_tables <- lavaan.pl:::cpp_pl(theta, with_tables$model, with_tables$data, 1L)
    pl_units <- lavaan.pl:::cpp_pl(theta, with_units$model, with_units$data, 1L)
    hj_tables <- lavaan.pl:::cpp_sandwich(theta, with_tables$model, with_tables$data, TRUE)
    hj_units <- lavaan.pl:::cpp_sandwich(theta, with_units$model, with_units$data, TRUE)
    expect_equal(pl_units$pl, pl_tables$pl, tolerance = 1e-10)
    expect_equal(pl_units$gradient, pl_tables$gradient, tolerance = 1e-8)
    expect_equal(hj_units$H, hj_tables$H, tolerance = 1e-8)
    expect_equal(hj_units$J, hj_tables$J, tolerance = 1e-8)
    expect_equal(hj_units$pairvar, hj_tables$pairvar, tolerance = 1e-8)
  }
})

#### The route is the one with fewer bivariate probabilities ####
local({
  use_tables <- lavaan.pl:::use_tables
  corners <- function(y, categories) {
    pairs <- combn(ncol(y), 2)
    sum(apply(pairs, 2, function(pair) {
      k <- pair[1]
      l <- pair[2]
      sum(mapply(function(a, b) {
        if (is.na(a) || is.na(b)) return(0)
        length(intersect(c(a, a + 1), seq_len(categories[k] - 1))) *
          length(intersect(c(b, b + 1), seq_len(categories[l] - 1)))
      }, y[, k], y[, l]))
    }))
  }
  set.seed(31)
  categories <- c(5, 3, 2, 4)
  y <- sapply(categories, function(m) sample(0:(m - 1), 200, TRUE))
  y[sample(length(y), 40)] <- NA
  a <- 2 - (y == 0) - sweep(y, 2, categories - 1, "==")
  a[is.na(a)] <- 0
  observations <- corners(y, categories)
  expect_equal(sum((rowSums(a)^2 - rowSums(a^2)) / 2), observations)
  grid <- sum(combn(categories - 1, 2, prod))
  crossover <- observations / grid
  stratum <- seq_len(nrow(y))
  expect_true(use_tables(floor(crossover) - 1, y, categories, stratum))
  expect_false(use_tables(ceiling(crossover) + 1, y, categories, stratum))
  expect_true(use_tables(1, y, categories, rep(0, nrow(y))))
  binary <- sapply(1:3, function(j) rep(0:1, 100))
  expect_true(use_tables(50, cbind(1:200, binary)[, -1], rep(2, 3), 1:200))
  expect_false(use_tables(2, binary, rep(2, 3), rep(0:1, 100)))
  expect_false(use_tables(200, binary, rep(2, 3), 1:200))
  expect_true(use_tables(10, cbind(rnorm(200), binary[, 1]), c(1, 2), rep(0:9, 20)))
})

#### Continuous covariates use the individual observations ####
local({
  set.seed(11)
  continuous <- cbind(rnorm(n), rbinom(n, 1, .5))
  items_continuous <- simulate_mimic(n, p, continuous, c(.5, -.3), 12)
  data_continuous <- lavaan.pl:::lavaanpl_data(items_continuous, X = continuous)
  expect_false(data_continuous$tables)
  expect_true(data$tables)
  expect_equal(length(data_continuous$freq), 0)
  model_continuous <- lavaan.pl:::lavaanpl_model(
    loadings_markers,
    PSI = matrix(NA),
    GAMMA = matrix(NA, 1, 2)
  )
  fit_continuous <- lavaan.pl:::fit_lavaanpl(data_continuous, model_continuous, CONTROL = set_control("QN"))
  vars <- lavaan.pl:::compute_var(fit_continuous)
  se <- sqrt(vars$asymptotic_variance / n)
  gamma <- as.vector(lavaan.pl:::vector_to_matrices(fit_continuous$model, fit_continuous$parameter_vector)[[1]]$gamma)
  expect_true(all(gamma / tail(se, 2) * sign(c(.5, -.3)) > 3))
  fit_sa <- lavaan.pl:::fit_lavaanpl(
    data_continuous,
    model_continuous,
    CONTROL = set_control("SA", pairs_per_iteration = 4)
  )
  expect_true(fit_sa$stoFit$converged)
  expect_true(max(abs(
    lavaan.pl:::transform_parameters(fit_sa$model, fit_sa$parameter_vector) - lavaan.pl:::transform_parameters(fit_continuous$model, fit_continuous$parameter_vector)
  ) / se) < 1)
})

#### The available-cases term is the weighted univariate likelihood given the covariates ####
local({
  with_missing <- items
  set.seed(13)
  with_missing[sample(length(items), 400)] <- NA
  model_gamma <- lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(NA, 1, 2))
  data_available <- lavaan.pl:::lavaanpl_data(with_missing, MISSING = "available.cases", X = binary)
  data_pairwise <- lavaan.pl:::lavaanpl_data(with_missing, MISSING = "pairwise", X = binary)
  theta <- suppressWarnings(lavaan.pl:::fit_lavaanpl(
    data_available,
    model_gamma,
    CONTROL = set_control("QN", maxeval = 1)
  ))$init
  pl <- function(x) {
    inputs <- lavaan.pl:::cpp_inputs(x, lavaan.pl:::lavaanpl_groups(list(model_gamma)))
    lavaan.pl:::cpp_pl(theta, inputs$model, inputs$data)$pl
  }
  par <- lavaan.pl:::vector_to_matrices(model_gamma, theta)
  slopes <- par$lambda %*% par$gamma
  cuts <- split(par$tau, rep(seq_len(p), data_available$cat - 1))
  univariate <- sum(sapply(seq_len(p), function(j) {
    observed <- !is.na(with_missing[, j])
    shift <- as.vector(binary[observed, , drop = FALSE] %*% slopes[j, ])
    bounds <- c(-Inf, cuts[[j]], Inf)
    y <- with_missing[observed, j]
    probability <- pnorm(bounds[y + 2] - shift) - pnorm(bounds[y + 1] - shift) + 1e-8
    sum(data_available$uweights[observed] * log(probability))
  }))
  expect_equal(pl(data_available) - pl(data_pairwise), univariate, tolerance = 1e-8)
  fit_available <- lavaan.pl:::fit_lavaanpl(data_available, model_gamma, CONTROL = set_control("QN"))
  expect_true(max(abs(
    lavaan.pl:::cpp_pl(fit_available$parameter_vector, fit_available$cpp_inputs$model, fit_available$cpp_inputs$data, 1L)$gradient
  )) / n < 1e-4)
})

#### Tests against the saturated and nested models ####
local({
  vars <- lavaan.pl:::compute_var(fit)
  test <- lavaan.pl:::compute_plrt(fit, vars)
  expect_equal(test$model$df, sum(data$cat - 1) + p * 2 + p * (p - 1) / 2 - fit$dims$npar)
  expect_equal(test$baseline$df, p * (p - 1) / 2)
  expect_true(test$saturated$logl >= lavaan.pl:::cpp_pl(fit$parameter_vector, fit$cpp_inputs$model, fit$cpp_inputs$data)$pl)
  gamma_fixed <- matrix(c(NA, 0), 1, 2)
  fit_nested <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA), GAMMA = gamma_fixed),
    CONTROL = set_control("QN")
  )
  difference <- lavaan.pl:::compute_plrt(fit_nested, H1 = fit)$model
  pl <- function(x) lavaan.pl:::cpp_pl(x$parameter_vector, x$cpp_inputs$model, x$cpp_inputs$data)$pl
  expect_equal(difference$stat, 2 * (pl(fit) - pl(fit_nested)))
  expect_equal(difference$df, 1)
  expect_true(difference$pvalue.scaled < .05)
  fit_without <- lavaan.pl:::fit_lavaanpl(data, lavaan.pl:::lavaanpl_model(loadings_markers, PSI = matrix(NA)), CONTROL = set_control("QN"))
  expect_equal(
    lavaan.pl:::compute_plrt(fit_without, H1 = fit_nested)$model$df,
    1
  )
  data_units <- lavaan.pl:::lavaanpl_data(items, X = cbind(binary, rnorm(n)))
  fit_units <- lavaan.pl:::fit_lavaanpl(
    data_units,
    lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(NA, 1, 3)),
    CONTROL = set_control("QN")
  )
  expect_false(data_units$tables)
  test_units <- lavaan.pl:::compute_plrt(fit_units, lavaan.pl:::compute_var(fit_units))
  expect_equal(test_units$model$df, sum(data$cat - 1) + p * 3 + p * (p - 1) / 2 - fit_units$dims$npar)
  expect_true(is.finite(test_units$model$stat.scaled))
  expect_true(is.finite(test_units$baseline$stat.scaled))
})

#### The lavaan interface ####
local({
  mimic <- as.data.frame(items)
  names(mimic) <- paste0("y", 1:p)
  for (item in names(mimic)) mimic[[item]] <- ordered(mimic[[item]])
  mimic$x1 <- binary[, 1]
  mimic$x2 <- binary[, 2]
  syntax <- "f =~ y1 + y2 + y3 + y4 + y5 + y6\nf ~ x1 + x2"
  fit_sem <- pmlsem(syntax, mimic, control = set_control("QN"))
  expect_equal(fit_sem@external$lavaanpl_internal$dims$ngam, 2)
  expect_equal(
    unname(coef(fit_sem)[c("f~x1", "f~x2")]),
    as.vector(lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]$gamma),
    tolerance = 1e-3
  )
  fit_direct <- pmlsem(paste0(syntax, "\ny2 ~ x1"), mimic)
  expect_equal(fit_direct@external$lavaanpl_internal$dims$nkap, 1)
  expect_equal(fit_direct@external$lavaanpl_internal$dims$q, 1)
  expect_true(is.finite(coef(fit_direct)["f=~y2"]))
  expect_error(pmlsem(syntax, mimic, conditional.x = FALSE), "conditional.x")
  expect_error(pmlsem(syntax, mimic, fixed.x = FALSE), "fixed.x")
  expect_error(
    pmlsem("f =~ y1 + y2 + y3 + y4 + y5 + y6\ny5 ~ y6", mimic),
    "exogenous covariates"
  )
  fit_sa <- pmlsem(syntax, mimic, control = set_control("SA", pairs_per_iteration = 4))
  expect_true(lavaan::lavInspect(fit_sa, "converged"))
  with_missing <- mimic
  with_missing$y2[1:100] <- NA
  fit_available <- pmlsem(syntax, with_missing, missing = "available.cases")
  expect_true(is.finite(fit_available@optim$logl))
  expect_equal(
    fit_available@optim$logl,
    lavaan.pl:::cpp_pl(
      fit_available@external$lavaanpl_internal$parameter_vector,
      fit_available@external$lavaanpl_internal$cpp_inputs$model,
      fit_available@external$lavaanpl_internal$cpp_inputs$data
    )$pl
  )
})

#### The location and scale of the covariates only change their effects and the thresholds ####
local({
  set.seed(21)
  age <- round(rnorm(n, 40, 10))
  items_age <- simulate_mimic(n, p, cbind(age - 40), -.05, 22)
  model_age <- lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(NA, 1, 1))
  fit_raw <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_age, X = cbind(age)), model_age, CONTROL = set_control("QN"))
  fit_centred <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_age, X = cbind(age - 40)), model_age, CONTROL = set_control("QN"))
  fit_scaled <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_age, X = cbind(1000 * age)), model_age, CONTROL = set_control("QN"))
  par_raw <- lavaan.pl:::vector_to_matrices(fit_raw$model, fit_raw$parameter_vector)[[1]]
  par_centred <- lavaan.pl:::vector_to_matrices(fit_centred$model, fit_centred$parameter_vector)[[1]]
  expect_equal(abs(par_raw$lambda), abs(par_centred$lambda), tolerance = 1e-4)
  expect_equal(abs(par_raw$gamma), abs(par_centred$gamma), tolerance = 1e-4)
  expect_equal(1000 * abs(lavaan.pl:::vector_to_matrices(fit_scaled$model, fit_scaled$parameter_vector)[[1]]$gamma), abs(par_raw$gamma), tolerance = 1e-4)
  slopes <- as.vector(par_raw$lambda %*% par_raw$gamma)
  expect_equal(
    par_raw$tau - rep(slopes, each = 3) * 40,
    par_centred$tau,
    tolerance = 1e-4
  )
  pl <- function(x) lavaan.pl:::cpp_pl(x$parameter_vector, x$cpp_inputs$model, x$cpp_inputs$data, 1L)
  expect_equal(pl(fit_raw)$pl, pl(fit_centred)$pl, tolerance = 1e-8)
  expect_true(max(abs(pl(fit_raw)$gradient)) / n < 1e-3)
  fit_started <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items_age, X = cbind(age)),
    model_age,
    CONTROL = set_control()
  )
  expect_equal(pl(fit_started)$pl, pl(fit_centred)$pl, tolerance = 1e-8)
  fit_sa <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items_age, X = cbind(age)),
    model_age,
    CONTROL = set_control("SA", pairs_per_iteration = 4)
  )
  expect_true(fit_sa$stoFit$converged)
  expect_equal(pl(fit_sa)$pl, pl(fit_centred)$pl, tolerance = 1e-4)
  fixed <- lavaan.pl:::lavaanpl_model(loadings_structure, GAMMA = matrix(0.01, 1, 1))
  fit_fixed <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_age, X = cbind(age)), fixed, CONTROL = set_control("QN"))
  expect_equal(as.vector(lavaan.pl:::vector_to_matrices(fit_fixed$model, fit_fixed$parameter_vector)[[1]]$gamma), 0.01)
  expect_true(max(abs(pl(fit_fixed)$gradient)) / n < 1e-3)
  expect_message(
    suppressWarnings(lavaan.pl:::fit_lavaanpl(
      lavaan.pl:::lavaanpl_data(items_age, X = cbind(age)),
      lavaan.pl:::lavaanpl_model(
        loadings_structure,
        GAMMA = matrix(NA, 1, 1),
        CONSTRAINTS = list(A = matrix(replace(numeric(25), 1, 1), 1), b = -1)
      ),
      CONTROL = set_control("QN")
    )),
    "consider centring"
  )
})

#### ucminf warns when it stops with a large gradient ####
local({
  expect_warning(
    lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_control("QN", maxeval = 1)),
    "may not be at an optimum"
  )
})
