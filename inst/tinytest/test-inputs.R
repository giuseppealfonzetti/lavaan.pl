source("fixtures/simulate.R")
#### Data and models ####
bfi_items <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1
two_factors <- matrix(0, 6, 2)
two_factors[1:3, 1] <- NA
two_factors[4:6, 2] <- NA

#### Control constructors return classed lists with defaults ####
local({
  ucminf_control <- set_control("QN")
  sa_control <- set_control("SA")
  expect_inherits(ucminf_control, "lavaanpl_qn_control")
  expect_inherits(sa_control, "lavaanpl_sa_control")
  expect_equal(
    unclass(ucminf_control),
    list(grtol = 1e-6, xtol = 1e-12, stepmax = 1, maxeval = 500L)
  )
  expect_identical(sa_control$average_epochs, 50L)
  expect_identical(sa_control$gain, "auto")
})

#### Control constructors validate their arguments ####
local({
  expect_error(set_control("QN", maxeval = 0))
  expect_error(set_control("QN", grtol = -1))
  expect_error(set_control("SA", precondition_epochs = -1))
  expect_error(set_control("SA", max_precondition_epochs = -1))
  expect_error(set_control(maxeval = 0))
  expect_error(set_control(alpha = 0))
  expect_error(set_control(nstarts = 0))
  expect_identical(set_control(nstarts = 2)$nstarts, 2L)
  expect_identical(lavaan.pl:::control_method(set_control()), "SA-QN")
  expect_error(set_control("SA", decay = 0.4))
  expect_error(set_control("SA", avg_from = 1))
  expect_error(set_control("SA", gain = "full"))
  for (method in c("SA-QN", "SA", "QN")) expect_identical(lavaan.pl:::control_method(set_control(method)), method)
  expect_error(set_control("nope"))
  expect_error(set_control("QN", average_epochs = 5), "unused argument")
})

#### lavaan.pl:::lavaanpl_data() computes the frequencies and validates the data ####
local({
  data <- lavaan.pl:::lavaanpl_data(bfi_items)
  expect_inherits(data, "lavaanpl_data")
  expect_identical(c(data$n, data$p), dim(bfi_items))
  expect_equal(unname(data$cat), rep(6, 6))
  expect_equal(data$freq, lavaan.pl:::cpp_pairs_freq(bfi_items, data$cat))
  expect_stdout(print(data), "<lavaanpl_data>")
  expect_error(lavaan.pl:::lavaanpl_data(bfi_items, MISSING = "complete"))
})

#### lavaan.pl:::lavaanpl_data() drops incomplete rows under listwise deletion ####
local({
  items_missing <- bfi_items
  set.seed(4)
  items_missing[sample(length(items_missing), 150)] <- NA
  expect_message(
    data <- lavaan.pl:::lavaanpl_data(items_missing),
    "incomplete observations removed"
  )
  expect_equal(data$n, sum(stats::complete.cases(items_missing)))
})

#### lavaan.pl:::lavaanpl_data() counts one univariate term per missing and observed item pair ####
local({
  items_missing <- bfi_items
  set.seed(4)
  items_missing[sample(length(items_missing), 150)] <- NA
  data <- lavaan.pl:::lavaanpl_data(items_missing, MISSING = "available.cases")
  missing_per_row <- rowSums(is.na(items_missing))
  expect_equal(
    sum(data$ufreq),
    sum(missing_per_row * (ncol(items_missing) - missing_per_row))
  )
})

#### lavaan.pl:::lavaanpl_data() normalises the weights within groups ####
local({
  group <- rep(1:2, c(1000, nrow(bfi_items) - 1000))
  set.seed(1)
  weights <- runif(nrow(bfi_items), .2, 3)
  within_groups <- lavaan.pl:::lavaanpl_data(
    bfi_items,
    WEIGHTS = weights,
    GROUP = group,
    WEIGHTS_NORMALIZATION = "group"
  )
  total <- lavaan.pl:::lavaanpl_data(bfi_items, WEIGHTS = weights, GROUP = group)
  expect_equal(
    as.vector(tapply(within_groups$weights, group, sum)),
    as.vector(table(group))
  )
  expect_equal(sum(total$weights), nrow(bfi_items))
})

#### lavaan.pl:::lavaanpl_data() numbers clusters within groups ####
local({
  group <- rep(1:2, each = 100)
  clusters <- rep(1:10, length.out = 200)
  data <- lavaan.pl:::lavaanpl_data(bfi_items[1:200, ], CLUSTER = clusters, GROUP = group)
  expect_equal(length(unique(data$cluster)), 20)
  expect_error(
    lavaan.pl:::lavaanpl_data(
      bfi_items[1:200, ],
      CLUSTER = rep(1:2, each = 100),
      GROUP = group
    ),
    "two clusters"
  )
})

#### lavaan.pl:::lavaanpl_model() builds the model and validates it ####
local({
  model <- lavaan.pl:::lavaanpl_model(two_factors)
  expect_inherits(model, "lavaanpl_model")
  expect_identical(c(model$p, model$q), c(6L, 2L))
  expect_identical(diag(model$PSI), c(1, 1))
  expect_stdout(print(model), "<lavaanpl_model>")
  expect_inherits(lavaan.pl:::lavaanpl_model(matrix(NA, 6, 1)), "lavaanpl_model")
  expect_inherits(
    lavaan.pl:::lavaanpl_model(
      two_factors,
      CONSTRAINTS = list(A = matrix(c(1, -1), 1, 2), b = 0)
    ),
    "lavaanpl_model"
  )
  expect_error(lavaan.pl:::lavaanpl_model(
    two_factors,
    CONSTRAINTS = list(A = matrix(1, 1, 2), b = c(0, 1))
  ))
  expect_error(lavaan.pl:::lavaanpl_model(two_factors, CONSTRAINTS = list(A = 1, b = 0)))
})

#### lavaan.pl:::lavaanpl_model() takes the loadings and latent variances as given, with unit variances by default ####
local({
  markers <- replace(two_factors, cbind(c(1, 4), 1:2), 1)
  model <- lavaan.pl:::lavaanpl_model(markers, PSI = matrix(NA, 2, 2))
  expect_identical(model$LAMBDA, markers)
  expect_true(all(is.na(diag(model$PSI))))
  expect_identical(diag(lavaan.pl:::lavaanpl_model(two_factors)$PSI), c(1, 1))
})

#### lavaan.pl:::lavaanpl_model() rejects inconsistent structures ####
local({
  three_factors <- build_constr_mat(9, 3, "simple")
  residual_covariances <- matrix(0, 9, 9)
  residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
  expect_error(lavaan.pl:::lavaanpl_model(three_factors, LATCOV = "orthogonal"))
  expect_error(lavaan.pl:::lavaanpl_model(three_factors, BETA = diag(3)), "diagonal")
  expect_error(
    lavaan.pl:::lavaanpl_model(
      three_factors,
      THETA = replace(residual_covariances, 2, 0.2)
    ),
    "symmetric"
  )
})

#### lavaan.pl:::cpp_pairs_freq() counts every pair once and skips missing responses ####
local({
  set.seed(123)
  p <- 10
  n <- 10
  binary <- matrix(as.numeric(rbinom(n * p, 1, .5)), n, p)
  frequencies <- lavaan.pl:::cpp_pairs_freq(binary, rep(2, p))
  pairs_per_row <- p * (p - 1) / 2
  expect_equal(length(frequencies), 4 * pairs_per_row)
  expect_equal(sum(frequencies), n * pairs_per_row)
  binary[1, 1] <- NA
  expect_equal(sum(lavaan.pl:::cpp_pairs_freq(binary, rep(2, p))), n * pairs_per_row - (p - 1))
  expect_equal(
    sum(lavaan.pl:::cpp_pairs_freq(binary, rep(2, p), WEIGHTS = rep(2, n))),
    2 * (n * pairs_per_row - (p - 1))
  )
})
