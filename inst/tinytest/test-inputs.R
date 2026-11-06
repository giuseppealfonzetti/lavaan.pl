#### Data and models ####
bfi_items <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1
two_factors <- matrix(0, 6, 2)
two_factors[1:3, 1] <- NA
two_factors[4:6, 2] <- NA

#### Control constructors return classed lists with defaults ####
local({
  ucminf_control <- set_qn_control()
  sa_control <- set_sa_control()
  expect_inherits(ucminf_control, "latpl_qn_control")
  expect_inherits(sa_control, "latpl_sa_control")
  expect_equal(
    unclass(ucminf_control),
    list(GRTOL = 1e-6, XTOL = 1e-12, STEPMAX = 1, MAXEVAL = 500L)
  )
  expect_identical(sa_control$AVERAGE_EPOCHS, 50L)
  expect_identical(sa_control$GAIN, "auto")
})

#### Control constructors validate their arguments ####
local({
  expect_error(set_qn_control(MAXEVAL = 0))
  expect_error(set_qn_control(GRTOL = -1))
  expect_error(set_sa_control(PRECONDITION_EPOCHS = -1))
  expect_error(set_sa_control(MAX_PRECONDITION_EPOCHS = -1))
  expect_error(set_saqn_control(MAXEVAL = 0))
  expect_error(set_saqn_control(ALPHA = 0))
  expect_error(set_saqn_control(NSTARTS = 0))
  expect_identical(set_saqn_control(NSTARTS = 2)$NSTARTS, 2L)
  expect_identical(lavaan.pl:::control_method(set_saqn_control()), "SA-QN")
  expect_error(set_sa_control(DECAY = 0.4))
  expect_error(set_sa_control(AVG_FROM = 1))
  expect_error(set_sa_control(GAIN = "full"))
})

#### lavaan.pl:::set_latpl_data() computes the frequencies and validates the data ####
local({
  data <- lavaan.pl:::set_latpl_data(bfi_items)
  expect_inherits(data, "latpl_data")
  expect_identical(c(data$n, data$p), dim(bfi_items))
  expect_equal(unname(data$cat), rep(6, 6))
  expect_equal(data$freq, lavaan.pl:::cpp_pairs_freq(bfi_items, data$cat))
  expect_stdout(print(data), "<latpl_data>")
  expect_error(lavaan.pl:::set_latpl_data(bfi_items, MISSING = "complete"))
})

#### lavaan.pl:::set_latpl_data() drops incomplete rows under listwise deletion ####
local({
  items_missing <- bfi_items
  set.seed(4)
  items_missing[sample(length(items_missing), 150)] <- NA
  expect_message(
    data <- lavaan.pl:::set_latpl_data(items_missing),
    "incomplete observations removed"
  )
  expect_equal(data$n, sum(stats::complete.cases(items_missing)))
})

#### lavaan.pl:::set_latpl_data() counts one univariate term per missing and observed item pair ####
local({
  items_missing <- bfi_items
  set.seed(4)
  items_missing[sample(length(items_missing), 150)] <- NA
  data <- lavaan.pl:::set_latpl_data(items_missing, MISSING = "available.cases")
  missing_per_row <- rowSums(is.na(items_missing))
  expect_equal(
    sum(data$ufreq),
    sum(missing_per_row * (ncol(items_missing) - missing_per_row))
  )
})

#### lavaan.pl:::set_latpl_data() normalises the weights within groups ####
local({
  group <- rep(1:2, c(1000, nrow(bfi_items) - 1000))
  set.seed(1)
  weights <- runif(nrow(bfi_items), .2, 3)
  within_groups <- lavaan.pl:::set_latpl_data(
    bfi_items,
    WEIGHTS = weights,
    GROUP = group,
    WEIGHTS_NORMALIZATION = "group"
  )
  total <- lavaan.pl:::set_latpl_data(bfi_items, WEIGHTS = weights, GROUP = group)
  expect_equal(
    as.vector(tapply(within_groups$weights, group, sum)),
    as.vector(table(group))
  )
  expect_equal(sum(total$weights), nrow(bfi_items))
})

#### lavaan.pl:::set_latpl_data() numbers clusters within groups ####
local({
  group <- rep(1:2, each = 100)
  clusters <- rep(1:10, length.out = 200)
  data <- lavaan.pl:::set_latpl_data(bfi_items[1:200, ], CLUSTER = clusters, GROUP = group)
  expect_equal(length(unique(data$cluster)), 20)
  expect_error(
    lavaan.pl:::set_latpl_data(
      bfi_items[1:200, ],
      CLUSTER = rep(1:2, each = 100),
      GROUP = group
    ),
    "two clusters"
  )
})

#### lavaan.pl:::set_latpl_model() builds the model and validates it ####
local({
  model <- lavaan.pl:::set_latpl_model(two_factors)
  expect_inherits(model, "latpl_model")
  expect_identical(c(model$p, model$q), c(6L, 2L))
  expect_identical(model$CONSTRLOGSD, c(0, 0))
  expect_stdout(print(model), "<latpl_model>")
  expect_inherits(lavaan.pl:::set_latpl_model(matrix(NA, 6, 1)), "latpl_model")
  expect_inherits(
    lavaan.pl:::set_latpl_model(
      two_factors,
      CONSTRAINTS = list(A = matrix(c(1, -1), 1, 2), b = 0)
    ),
    "latpl_model"
  )
  expect_error(lavaan.pl:::set_latpl_model(
    two_factors,
    CONSTRAINTS = list(A = matrix(1, 1, 2), b = c(0, 1))
  ))
  expect_error(lavaan.pl:::set_latpl_model(two_factors, CONSTRAINTS = list(A = 1, b = 0)))
})

#### lavaan.pl:::set_latpl_model() fixes the first loading of each factor under marker identification ####
local({
  marker <- lavaan.pl:::set_latpl_model(two_factors, STDLV = FALSE)
  expect_identical(marker$CONSTRMAT[c(1, 4), ], diag(2) * 1)
  expect_true(all(is.na(marker$CONSTRLOGSD)))
})

#### lavaan.pl:::set_latpl_model() rejects inconsistent structures ####
local({
  three_factors <- lavaan.pl:::build_constr_mat(9, 3, "simple")
  residual_covariances <- matrix(0, 9, 9)
  residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
  expect_error(
    lavaan.pl:::set_latpl_model(three_factors, CORRFLAG = FALSE, LATCOV = "direct"),
    "CORRFLAG"
  )
  expect_error(lavaan.pl:::set_latpl_model(three_factors, REGMAT = diag(3)), "diagonal")
  expect_error(
    lavaan.pl:::set_latpl_model(
      three_factors,
      RESCOV = replace(residual_covariances, 2, 0.2)
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

#### lavaan.pl:::build_constr_mat() builds triangular and simple structures ####
local({
  p <- 10
  q <- 2
  free_triangular <- is.na(lavaan.pl:::build_constr_mat(
    P = p,
    Q = q,
    STRUCT = "triangular"
  ))
  expect_equal(sum(free_triangular), (p * q) - (q * (q - 1) / 2))
  expect_equal(
    sum(free_triangular[lower.tri(free_triangular, diag = TRUE)]),
    (p * q) - (q * (q - 1) / 2)
  )
  expect_equal(
    sum(free_triangular[upper.tri(free_triangular, diag = FALSE)]),
    0
  )
  expect_equal(sum(free_triangular %in% c(0, 1)), p * q)
  free_simple <- is.na(lavaan.pl:::build_constr_mat(P = p, Q = q, STRUCT = "simple"))
  expect_equal(sum(free_simple), p)
  expect_equal(sum(rowSums(free_simple) == 1), p)
})

#### lavaan.pl:::gen_loadings() rejects invalid constraint matrices ####
local({
  expect_error(lavaan.pl:::gen_loadings(CONSTRMAT = 2))
  expect_error(lavaan.pl:::gen_loadings(CONSTRMAT = NA))
})
