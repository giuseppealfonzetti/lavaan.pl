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
data <- set_latpl_data(items)
lsat_data <- set_latpl_data(as.matrix(data.frame(lapply(LSAT, as.numeric))) - 1)
one_factor <- set_latpl_model(matrix(NA, 5, 1))
items_two_factors <- items[1:300, 4:9]
two_factors <- set_latpl_model(build_constr_mat(6, 2, "simple"))

#### fit_latpl() validates its inputs ####
local({
  expect_error(fit_latpl(items, set_latpl_model(three_factors)), "inherits")
  expect_error(
    fit_latpl(data, set_latpl_model(matrix(NA, 5, 1))),
    "DATA\\$p == MODEL\\$p"
  )
  expect_error(
    fit_latpl(
      data,
      set_latpl_model(three_factors),
      INIT = set_init(THETA = rep(0, 3))
    ),
    "length"
  )
  expect_error(
    fit_latpl(lsat_data, one_factor, CONTROL = list(MAXEVAL = 5)),
    "inherits"
  )
})

#### fit_latpl() stores its inputs and compute_var() covers the free parameters ####
local({
  model <- set_latpl_model(three_factors)
  fit <- fit_latpl(data, model)
  vars <- compute_var(fit)
  expect_identical(fit$data, data)
  expect_identical(fit$model, model)
  expect_inherits(fit$control, "latpl_ucminf_control")
  expect_identical(fit$method, "ucminf")
  expect_equal(dim(vars$H), c(fit$dims$npar, fit$dims$npar))
  expect_true(all(is.finite(vars$asymptotic_variance)))
})

#### The control object selects the optimiser ####
local({
  fit_sa <- fit_latpl(
    lsat_data,
    one_factor,
    CONTROL = set_sa_control(EPOCHS = 1, PAIRS_PER_ITERATION = 2)
  )
  fit_short <- fit_latpl(
    lsat_data,
    one_factor,
    CONTROL = set_ucminf_control(MAXEVAL = 5)
  )
  expect_identical(fit_sa$method, "SA")
  expect_true(fit_short$numFit$info[["neval"]] <= 5)
})

#### Initialisation by SA reaches the same estimate from different starting values ####
local({
  fit_sa_init <- fit_latpl(
    lsat_data,
    one_factor,
    INIT = set_init(
      METHOD = "SA",
      CONTROL = set_sa_init_control(ADAM_EPOCHS = 100)
    )
  )
  fit_standard_init <- fit_latpl(lsat_data, one_factor)
  expect_false(identical(fit_sa_init$init, fit_standard_init$init))
  expect_equal(fit_sa_init$theta, fit_standard_init$theta, tolerance = 1e-4)
})

#### Numerical and analytic Hessians agree ####
local({
  fit <- fit_latpl(lsat_data, one_factor)
  expect_equal(
    compute_var(fit, NUMDERIV = TRUE)$invH,
    compute_var(fit)$invH,
    tolerance = 1e-5
  )
})

#### Linear constraints hold at the estimate, where the projected gradient vanishes ####
local({
  names_free <- theta_names(data, set_latpl_model(three_factors))
  constraints <- rbind(
    replace(
      numeric(length(names_free)),
      match(c("lambda[2,1]", "lambda[3,1]"), names_free),
      c(1, -1)
    ),
    replace(numeric(length(names_free)), match("tau[4,2]", names_free), 1)
  )
  fit <- fit_latpl(
    data,
    set_latpl_model(
      three_factors,
      CONSTRAINTS = list(A = constraints, b = c(0, 0.25))
    )
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
    fit$theta[fit$lincon$free_idx]
  )
  expect_equal(fit$dims$npar, length(names_free) - 2)
  expect_equal(as.vector(constraints %*% fit$theta), c(0, 0.25))
  expect_true(
    max(abs(crossprod(free_to_full, negative_gradient(fit$theta)))) < 1e-5
  )
  expect_equal(
    compute_var(fit)$H,
    (hessian_numeric + t(hessian_numeric)) / 2,
    tolerance = 1e-4
  )
})

#### Invalid linear constraints are rejected ####
local({
  names_free <- theta_names(data, set_latpl_model(three_factors))
  constraint <- function(positions, values) {
    matrix(replace(numeric(length(names_free)), positions, values), 1)
  }
  equal_loadings <- constraint(
    match(c("lambda[2,1]", "lambda[3,1]"), names_free),
    c(1, -1)
  )
  expect_error(
    fit_latpl(
      data,
      set_latpl_model(
        three_factors,
        CONSTRAINTS = list(
          A = constraint(match("eta[1,2]", names_free), 1),
          b = 0
        )
      )
    ),
    "LATCOV"
  )
  expect_error(
    fit_latpl(
      data,
      set_latpl_model(
        three_factors,
        CONSTRAINTS = list(A = constraint(c(1, 28), c(1, -1)), b = 0)
      )
    ),
    "either thresholds"
  )
  expect_error(
    fit_latpl(
      data,
      set_latpl_model(
        three_factors,
        CONSTRAINTS = list(
          A = rbind(equal_loadings, equal_loadings),
          b = c(0, 1)
        )
      )
    ),
    "inconsistent"
  )
})

#### Weights: constant weights leave the estimates unchanged, frequency weights equal replicated rows ####
local({
  set.seed(2)
  weights <- sample(1:3, nrow(items_two_factors), replace = TRUE)
  fit <- fit_latpl(set_latpl_data(items_two_factors), two_factors)
  fit_constant <- fit_latpl(
    set_latpl_data(
      items_two_factors,
      WEIGHTS = rep(2, nrow(items_two_factors))
    ),
    two_factors
  )
  fit_weighted <- fit_latpl(
    set_latpl_data(
      items_two_factors,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency"
    ),
    two_factors
  )
  fit_replicated <- fit_latpl(
    set_latpl_data(items_two_factors[
      rep(seq_len(nrow(items_two_factors)), weights),
    ]),
    two_factors
  )
  expect_equal(fit_constant$theta, fit$theta)
  expect_equal(fit_weighted$theta, fit_replicated$theta, tolerance = 1e-6)
})

#### Complete data give the same estimates under every missing-data option ####
local({
  fit_listwise <- fit_latpl(set_latpl_data(items_two_factors), two_factors)
  fit_pairwise <- fit_latpl(
    set_latpl_data(items_two_factors, MISSING = "pairwise"),
    two_factors
  )
  fit_available <- fit_latpl(
    set_latpl_data(items_two_factors, MISSING = "available.cases"),
    two_factors
  )
  expect_equal(fit_pairwise$theta, fit_listwise$theta)
  expect_equal(fit_available$theta, fit_listwise$theta)
})

#### Two unconstrained groups reproduce the separate fits ####
local({
  items_first <- items[1:150, 4:9]
  items_second <- items[151:300, 4:9]
  fit_first <- fit_latpl(set_latpl_data(items_first), two_factors)
  fit_second <- fit_latpl(set_latpl_data(items_second), two_factors)
  data_two_groups <- set_latpl_data(
    rbind(items_first, items_second),
    GROUP = rep(c("b", "a"), each = 150)
  )
  fit_two_groups <- fit_latpl(
    data_two_groups,
    set_latpl_groups(list(two_factors, two_factors))
  )
  pl <- function(fit) {
    lavaan.pl:::cpp_pl(fit$theta, fit$engine$model, fit$engine$data)$pl
  }
  expect_identical(data_two_groups$group, rep(c(0, 1), each = 150))
  expect_equal(
    fit_two_groups$theta,
    c(fit_first$theta, fit_second$theta),
    tolerance = 1e-5
  )
  expect_equal(pl(fit_two_groups), pl(fit_first) + pl(fit_second))
  expect_equal(
    compute_var(fit_two_groups)$asymptotic_variance / 300,
    c(
      compute_var(fit_first)$asymptotic_variance,
      compute_var(fit_second)$asymptotic_variance
    ) /
      150,
    tolerance = 1e-4
  )
  expect_length(get_par(fit_two_groups), 2)
  expect_equal(
    get_par(fit_two_groups)[[2]],
    theta2par(two_factors, fit_second$theta),
    tolerance = 1e-4
  )
  expect_equal(
    par2theta(fit_two_groups$model, get_par(fit_two_groups)),
    fit_two_groups$theta
  )
})

#### An equality constraint across groups holds and removes one parameter ####
local({
  data_two_groups <- set_latpl_data(
    items_two_factors,
    GROUP = rep(1:2, each = 150)
  )
  model_two_groups <- set_latpl_groups(list(two_factors, two_factors))
  names_free <- theta_names(data_two_groups, model_two_groups)
  equal_loadings <- matrix(0, 1, length(names_free))
  equal_loadings[1, match(c("lambda[1,1]", "lambda[1,1].g2"), names_free)] <- c(
    1,
    -1
  )
  fit_free <- fit_latpl(data_two_groups, model_two_groups)
  fit_equal <- fit_latpl(
    data_two_groups,
    set_latpl_groups(
      list(two_factors, two_factors),
      list(A = equal_loadings, b = 0)
    )
  )
  expect_equal(
    fit_equal$theta[match("lambda[1,1]", names_free)],
    fit_equal$theta[match("lambda[1,1].g2", names_free)]
  )
  expect_equal(fit_equal$dims$npar, fit_free$dims$npar - 1)
})

#### The number of groups must match between data and model ####
local({
  data_two_groups <- set_latpl_data(
    items_two_factors,
    GROUP = rep(1:2, each = 150)
  )
  expect_error(fit_latpl(data_two_groups, two_factors), "groups")
  expect_error(
    fit_latpl(
      set_latpl_data(items_two_factors),
      set_latpl_groups(list(two_factors, two_factors))
    ),
    "groups"
  )
})

#### Frequency weights equal replicated rows within groups ####
local({
  set.seed(2)
  weights <- sample(1:3, 300, replace = TRUE)
  group <- rep(1:2, each = 150)
  model_two_groups <- set_latpl_groups(list(two_factors, two_factors))
  fit_weighted <- fit_latpl(
    set_latpl_data(
      items_two_factors,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency",
      GROUP = group
    ),
    model_two_groups
  )
  fit_replicated <- fit_latpl(
    set_latpl_data(
      items_two_factors[rep(1:300, weights), ],
      GROUP = group[rep(1:300, weights)]
    ),
    model_two_groups
  )
  expect_equal(fit_weighted$theta, fit_replicated$theta, tolerance = 1e-6)
})

#### Weights normalised within groups do not depend on the scale of each group's weights ####
local({
  group <- rep(1:2, each = 150)
  set.seed(3)
  weights <- runif(300, .2, 3)
  rescaled <- weights * ifelse(group == 2, 5, 1)
  model_two_groups <- set_latpl_groups(list(two_factors, two_factors))
  names_free <- theta_names(
    set_latpl_data(items_two_factors, GROUP = group),
    model_two_groups
  )
  equal_loadings <- matrix(0, 1, length(names_free))
  equal_loadings[
    1,
    match(c("lambda[1,1]", "lambda[1,1].g2"), names_free)
  ] <- c(1, -1)
  model_equal <- set_latpl_groups(
    list(two_factors, two_factors),
    list(A = equal_loadings, b = 0)
  )
  estimates <- function(weights, normalization) {
    data <- set_latpl_data(
      items_two_factors,
      WEIGHTS = weights,
      GROUP = group,
      WEIGHTS_NORMALIZATION = normalization
    )
    fit_latpl(data, model_equal)$theta
  }
  expect_equal(estimates(rescaled, "group"), estimates(weights, "group"))
  expect_false(isTRUE(all.equal(
    estimates(rescaled, "total"),
    estimates(weights, "total")
  )))
})
