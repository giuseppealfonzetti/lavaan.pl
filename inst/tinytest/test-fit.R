source("fixtures/simulate.R")
#### Data and models ####
three_factors <- build_constr_mat(9, 3, "simple")
latent_correlations <- matrix(.3, 3, 3)
diag(latent_correlations) <- 1
items <- sim_data(
  600,
  gen_loadings(three_factors, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
data <- lavaan.pl:::lavaanpl_data(items)
lsat_data <- lavaan.pl:::lavaanpl_data(as.matrix(data.frame(lapply(LSAT, as.numeric))) - 1)
one_factor <- lavaan.pl:::lavaanpl_model(matrix(NA, 5, 1))
items_two_factors <- items[1:300, 4:9]
two_factors <- lavaan.pl:::lavaanpl_model(build_constr_mat(6, 2, "simple"))

#### lavaan.pl:::fit_lavaanpl(, CONTROL = set_qn_control()) validates its inputs ####
local({
  expect_error(lavaan.pl:::fit_lavaanpl(items, lavaan.pl:::lavaanpl_model(three_factors), CONTROL = set_qn_control()), "inherits")
  expect_error(
    lavaan.pl:::fit_lavaanpl(data, lavaan.pl:::lavaanpl_model(matrix(NA, 5, 1)), CONTROL = set_qn_control()),
    "DATA\\$p == MODEL\\$p"
  )
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      data,
      lavaan.pl:::lavaanpl_model(three_factors),
      START = rep(0, 3),
      CONTROL = set_qn_control()
    ),
    "length"
  )
  expect_error(
    lavaan.pl:::fit_lavaanpl(lsat_data, one_factor, CONTROL = list(MAXEVAL = 5)),
    "inherits"
  )
})

#### lavaan.pl:::fit_lavaanpl(, CONTROL = set_qn_control()) stores its inputs and lavaan.pl:::compute_var() covers the free parameters ####
local({
  model <- lavaan.pl:::lavaanpl_model(three_factors)
  fit <- lavaan.pl:::fit_lavaanpl(data, model, CONTROL = set_qn_control())
  vars <- lavaan.pl:::compute_var(fit)
  expect_identical(fit$data, data)
  expect_identical(fit$model$MODELS[[1]], model)
  expect_inherits(fit$control, "lavaanpl_qn_control")
  expect_identical(fit$method, "QN")
  expect_equal(dim(vars$H), c(fit$dims$npar, fit$dims$npar))
  expect_true(all(is.finite(vars$asymptotic_variance)))
})

#### The control object selects the optimiser ####
local({
  fit_sa <- lavaan.pl:::fit_lavaanpl(
    lsat_data,
    one_factor,
    CONTROL = set_sa_control(average_epochs = 1, pairs_per_iteration = 2)
  )
  fit_short <- lavaan.pl:::fit_lavaanpl(
    lsat_data,
    one_factor,
    CONTROL = set_qn_control(maxeval = 5)
  )
  expect_identical(fit_sa$method, "SA")
  expect_true(fit_short$numFit$info[["neval"]] <= 5)
})

#### Numerical and analytic Hessians agree ####
local({
  fit <- lavaan.pl:::fit_lavaanpl(lsat_data, one_factor, CONTROL = set_qn_control())
  expect_equal(
    lavaan.pl:::compute_var(fit, NUMDERIV = TRUE)$H,
    lavaan.pl:::compute_var(fit)$H,
    tolerance = 1e-5
  )
})

#### Linear constraints hold at the estimate, where the projected gradient vanishes ####
local({
  names_free <- lavaan.pl:::parameter_names(data, lavaan.pl:::lavaanpl_model(three_factors))
  constraints <- rbind(
    replace(
      numeric(length(names_free)),
      match(c("lambda[2,1]", "lambda[3,1]"), names_free),
      c(1, -1)
    ),
    replace(numeric(length(names_free)), match("tau[4,2]", names_free), 1)
  )
  fit <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(
      three_factors,
      CONSTRAINTS = list(A = constraints, b = c(0, 0.25))
    ),
    CONTROL = set_qn_control()
  )
  free_to_full <- fit$lincon$K
  negative_gradient <- function(theta) {
    -lavaan.pl:::cpp_pl(theta, fit$engine$model, fit$engine$data, 1L)$gradient /
      data$n
  }
  hessian_numeric <- numDeriv::jacobian(
    function(x) {
      as.vector(crossprod(
        free_to_full,
        negative_gradient(as.vector(free_to_full %*% x + fit$lincon$c))
      ))
    },
    fit$parameter_vector[fit$lincon$free_idx]
  )
  expect_equal(fit$dims$npar, length(names_free) - 2)
  expect_equal(as.vector(constraints %*% fit$parameter_vector), c(0, 0.25))
  expect_true(
    max(abs(crossprod(free_to_full, negative_gradient(fit$parameter_vector)))) < 1e-5
  )
  expect_equal(
    lavaan.pl:::compute_var(fit)$H,
    (hessian_numeric + t(hessian_numeric)) / 2,
    tolerance = 1e-4
  )
})

#### Invalid linear constraints are rejected ####
local({
  names_free <- lavaan.pl:::parameter_names(data, lavaan.pl:::lavaanpl_model(three_factors))
  constraint <- function(positions, values) {
    matrix(replace(numeric(length(names_free)), positions, values), 1)
  }
  equal_loadings <- constraint(
    match(c("lambda[2,1]", "lambda[3,1]"), names_free),
    c(1, -1)
  )
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      data,
      lavaan.pl:::lavaanpl_model(
        three_factors,
        CONSTRAINTS = list(
          A = constraint(match("eta[1,2]", names_free), 1),
          b = 0
        )
      ),
      CONTROL = set_qn_control()
    ),
    "LATCOV"
  )
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      data,
      lavaan.pl:::lavaanpl_model(
        three_factors,
        CONSTRAINTS = list(A = constraint(c(1, 28), c(1, -1)), b = 0)
      ),
      CONTROL = set_qn_control()
    ),
    "either thresholds"
  )
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      data,
      lavaan.pl:::lavaanpl_model(
        three_factors,
        CONSTRAINTS = list(
          A = rbind(equal_loadings, equal_loadings),
          b = c(0, 1)
        )
      ),
      CONTROL = set_qn_control()
    ),
    "inconsistent"
  )
})

#### Weights: constant weights leave the estimates unchanged, frequency weights equal replicated rows ####
local({
  set.seed(2)
  weights <- sample(1:3, nrow(items_two_factors), replace = TRUE)
  fit <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_two_factors), two_factors, CONTROL = set_qn_control())
  fit_constant <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(
      items_two_factors,
      WEIGHTS = rep(2, nrow(items_two_factors))
    ),
    two_factors,
    CONTROL = set_qn_control()
  )
  fit_weighted <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(
      items_two_factors,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency"
    ),
    two_factors,
    CONTROL = set_qn_control()
  )
  fit_replicated <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items_two_factors[
      rep(seq_len(nrow(items_two_factors)), weights),
    ]),
    two_factors,
    CONTROL = set_qn_control()
  )
  expect_equal(fit_constant$parameter_vector, fit$parameter_vector)
  expect_equal(fit_weighted$parameter_vector, fit_replicated$parameter_vector, tolerance = 1e-6)
})

#### Complete data give the same estimates under every missing-data option ####
local({
  fit_listwise <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_two_factors), two_factors, CONTROL = set_qn_control())
  fit_pairwise <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items_two_factors, MISSING = "pairwise"),
    two_factors,
    CONTROL = set_qn_control()
  )
  fit_available <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(items_two_factors, MISSING = "available.cases"),
    two_factors,
    CONTROL = set_qn_control()
  )
  expect_equal(fit_pairwise$parameter_vector, fit_listwise$parameter_vector)
  expect_equal(fit_available$parameter_vector, fit_listwise$parameter_vector)
})

#### Two unconstrained groups reproduce the separate fits ####
local({
  items_first <- items[1:150, 4:9]
  items_second <- items[151:300, 4:9]
  fit_first <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_first), two_factors, CONTROL = set_qn_control())
  fit_second <- lavaan.pl:::fit_lavaanpl(lavaan.pl:::lavaanpl_data(items_second), two_factors, CONTROL = set_qn_control())
  data_two_groups <- lavaan.pl:::lavaanpl_data(
    rbind(items_first, items_second),
    GROUP = rep(c("b", "a"), each = 150)
  )
  fit_two_groups <- lavaan.pl:::fit_lavaanpl(
    data_two_groups,
    lavaan.pl:::lavaanpl_groups(list(two_factors, two_factors)),
    CONTROL = set_qn_control()
  )
  pl <- function(fit) {
    lavaan.pl:::cpp_pl(fit$parameter_vector, fit$engine$model, fit$engine$data)$pl
  }
  expect_identical(data_two_groups$group, rep(c(0, 1), each = 150))
  expect_equal(
    fit_two_groups$parameter_vector,
    c(fit_first$parameter_vector, fit_second$parameter_vector),
    tolerance = 1e-5
  )
  expect_equal(pl(fit_two_groups), pl(fit_first) + pl(fit_second))
  expect_equal(
    lavaan.pl:::compute_var(fit_two_groups)$asymptotic_variance / 300,
    c(
      lavaan.pl:::compute_var(fit_first)$asymptotic_variance,
      lavaan.pl:::compute_var(fit_second)$asymptotic_variance
    ) /
      150,
    tolerance = 1e-4
  )
  expect_length(lavaan.pl:::vector_to_matrices(fit_two_groups$model, fit_two_groups$parameter_vector), 2)
  expect_equal(
    lavaan.pl:::vector_to_matrices(fit_two_groups$model, fit_two_groups$parameter_vector)[[2]],
    lavaan.pl:::vector_to_matrices(two_factors, fit_second$parameter_vector),
    tolerance = 1e-4
  )
  expect_equal(
    lavaan.pl:::matrices_to_vector(fit_two_groups$model, lavaan.pl:::vector_to_matrices(fit_two_groups$model, fit_two_groups$parameter_vector)),
    fit_two_groups$parameter_vector
  )
})

#### An equality constraint across groups holds and removes one parameter ####
local({
  data_two_groups <- lavaan.pl:::lavaanpl_data(
    items_two_factors,
    GROUP = rep(1:2, each = 150)
  )
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(two_factors, two_factors))
  names_free <- lavaan.pl:::parameter_names(data_two_groups, model_two_groups)
  equal_loadings <- matrix(0, 1, length(names_free))
  equal_loadings[1, match(c("lambda[1,1]", "lambda[1,1].g2"), names_free)] <- c(
    1,
    -1
  )
  fit_free <- lavaan.pl:::fit_lavaanpl(data_two_groups, model_two_groups, CONTROL = set_qn_control())
  fit_equal <- lavaan.pl:::fit_lavaanpl(
    data_two_groups,
    lavaan.pl:::lavaanpl_groups(
      list(two_factors, two_factors),
      list(A = equal_loadings, b = 0)
    ),
    CONTROL = set_qn_control()
  )
  expect_equal(
    fit_equal$parameter_vector[match("lambda[1,1]", names_free)],
    fit_equal$parameter_vector[match("lambda[1,1].g2", names_free)]
  )
  expect_equal(fit_equal$dims$npar, fit_free$dims$npar - 1)
})

#### The number of groups must match between data and model ####
local({
  data_two_groups <- lavaan.pl:::lavaanpl_data(
    items_two_factors,
    GROUP = rep(1:2, each = 150)
  )
  expect_error(lavaan.pl:::fit_lavaanpl(data_two_groups, two_factors, CONTROL = set_qn_control()), "groups")
  expect_error(
    lavaan.pl:::fit_lavaanpl(
      lavaan.pl:::lavaanpl_data(items_two_factors),
      lavaan.pl:::lavaanpl_groups(list(two_factors, two_factors)),
      CONTROL = set_qn_control()
    ),
    "groups"
  )
})

#### Frequency weights equal replicated rows within groups ####
local({
  set.seed(2)
  weights <- sample(1:3, 300, replace = TRUE)
  group <- rep(1:2, each = 150)
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(two_factors, two_factors))
  fit_weighted <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(
      items_two_factors,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency",
      GROUP = group
    ),
    model_two_groups,
    CONTROL = set_qn_control()
  )
  fit_replicated <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(
      items_two_factors[rep(1:300, weights), ],
      GROUP = group[rep(1:300, weights)]
    ),
    model_two_groups,
    CONTROL = set_qn_control()
  )
  expect_equal(fit_weighted$parameter_vector, fit_replicated$parameter_vector, tolerance = 1e-6)
})

#### Weights normalised within groups do not depend on the scale of each group's weights ####
local({
  group <- rep(1:2, each = 150)
  set.seed(3)
  weights <- runif(300, .2, 3)
  rescaled <- weights * ifelse(group == 2, 5, 1)
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(two_factors, two_factors))
  names_free <- lavaan.pl:::parameter_names(
    lavaan.pl:::lavaanpl_data(items_two_factors, GROUP = group),
    model_two_groups
  )
  equal_loadings <- matrix(0, 1, length(names_free))
  equal_loadings[
    1,
    match(c("lambda[1,1]", "lambda[1,1].g2"), names_free)
  ] <- c(1, -1)
  model_equal <- lavaan.pl:::lavaanpl_groups(
    list(two_factors, two_factors),
    list(A = equal_loadings, b = 0)
  )
  estimates <- function(weights, normalization) {
    data <- lavaan.pl:::lavaanpl_data(
      items_two_factors,
      WEIGHTS = weights,
      GROUP = group,
      WEIGHTS_NORMALIZATION = normalization
    )
    lavaan.pl:::fit_lavaanpl(data, model_equal, CONTROL = set_qn_control())$parameter_vector
  }
  expect_equal(estimates(rescaled, "group"), estimates(weights, "group"))
  expect_false(isTRUE(all.equal(
    estimates(rescaled, "total"),
    estimates(weights, "total")
  )))
})
