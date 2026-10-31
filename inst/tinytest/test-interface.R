D <- set_latpl_data(as.matrix(data.frame(lapply(LSAT, as.numeric))) - 1)
cl <- set_latpl_model(matrix(NA, 5, 1))

# control constructors returns classed lists with defaults
local({
  uc <- set_ucminf_control()
  sa <- set_sa_control()
  init <- set_sa_init_control()
  expect_inherits(uc, "latpl_ucminf_control")
  expect_inherits(sa, "latpl_sa_control")
  expect_inherits(init, "latpl_sa_control")
  expect_equal(
    unclass(uc),
    list(GRTOL = 1e-6, XTOL = 1e-12, STEPMAX = 1, MAXEVAL = 500L, ROBUST = Inf)
  )
  expect_identical(sa$EPOCHS, 20L)
  expect_identical(sa$GAIN, "auto")
  expect_identical(init$EPOCHS, 0L)
  expect_identical(init$BURN_EPOCHS, 5L)
  expect_identical(names(sa), names(init))
})

# control constructors validate their arguments
local({
  expect_error(set_ucminf_control(MAXEVAL = 0))
  expect_error(set_ucminf_control(GRTOL = -1))
  expect_error(set_sa_control(BURN_EPOCHS = -1))
  expect_error(set_sa_control(DECAY = 0.4))
  expect_error(set_sa_control(AVG_FROM = 1))
  expect_error(set_sa_control(GAIN = "full"))
  expect_error(set_sa_init_control(ALPHA = 0))
})

# fit_latpl() checks control objects
local({
  expect_error(fit_latpl(D, cl, CONTROL = list(MAXEVAL = 5)), "inherits")
  fit <- fit_latpl(
    D,
    cl,
    CONTROL = set_sa_control(EPOCHS = 1, PAIRS_PER_ITERATION = 2)
  )
  expect_identical(fit$method, "SA")
  fit <- fit_latpl(D, cl, CONTROL = set_ucminf_control(MAXEVAL = 5))
  expect_true(fit$numFit$info[["neval"]] <= 5)
})

# cfa() checks  control objects
local({
  mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
  expect_error(
    cfa(mod, LSAT, std.lv = TRUE, control = list(ctrl = list(maxeval = 5))),
    "set_ucminf_control"
  )
  expect_error(
    cfa(mod, LSAT, std.lv = TRUE, estimator.args = list(method = "SA")),
    "set_sa_control"
  )
  expect_error(
    cfa(
      mod,
      LSAT,
      std.lv = TRUE,
      estimator.args = list(cpp_control_init = list(MAXE = 2))
    ),
    "set_init"
  )
  expect_error(
    cfa(
      mod,
      LSAT,
      std.lv = TRUE,
      estimator.args = list(init_method = "standard")
    ),
    "set_init"
  )
  expect_error(
    cfa(
      mod,
      LSAT,
      std.lv = TRUE,
      estimator.args = list(init = set_init(THETA = rep(0, 10)))
    ),
    "start"
  )
  fit <- cfa(
    mod,
    LSAT,
    std.lv = TRUE,
    se = "none",
    control = set_ucminf_control(MAXEVAL = 5)
  )
  expect_true(fit@external$lavaan.pl$numFit$info[["neval"]] <= 5)
})

Y <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1

# set_latpl_data() validates  data and computes frequencies
local({
  d <- set_latpl_data(Y)
  expect_inherits(d, "latpl_data")
  expect_identical(c(d$n, d$p), dim(Y))
  expect_equal(unname(d$cat), rep(6, 6))
  expect_equal(d$freq, pairs_freq(Y, d$cat))
  expect_null(d$valfreq)
  expect_false(is.null(set_latpl_data(Y, VALDATA = Y[1:100, ])$valfreq))
  expect_stdout(print(d), "<latpl_data>")

  Yna <- Y
  Yna[1, 1] <- NA
  expect_message(
    d_na <- set_latpl_data(Yna),
    "1 incomplete observations removed"
  )
  expect_equal(d_na$n, nrow(Y) - 1)
  expect_error(set_latpl_data(Yna, MISSING = "complete"))
  expect_error(set_latpl_data(Y, VALDATA = Y[, 1:5]))
  expect_error(set_latpl_data(Y, VALDATA = Y + 1))
})

# set_latpl_model() validates model
local({
  A <- matrix(0, 6, 2)
  A[1:3, 1] <- NA
  A[4:6, 2] <- NA
  m <- set_latpl_model(A)
  expect_inherits(m, "latpl_model")
  expect_identical(c(m$p, m$q), c(6L, 2L))
  expect_identical(m$CONSTRLOGSD, c(0, 0))
  expect_stdout(print(m), "<latpl_model>")

  marker <- set_latpl_model(A, STDLV = FALSE)
  expect_identical(marker$CONSTRMAT[c(1, 4), ], diag(2) * 1)
  expect_true(all(is.na(marker$CONSTRLOGSD)))

  expect_inherits(set_latpl_model(matrix(NA, 6, 1)), "latpl_model")
  expect_inherits(
    set_latpl_model(A, CONSTRAINTS = list(A = matrix(c(1, -1), 1, 2), b = 0)),
    "latpl_model"
  )
  expect_error(set_latpl_model(
    A,
    CONSTRAINTS = list(A = matrix(1, 1, 2), b = c(0, 1))
  ))
  expect_error(set_latpl_model(A, CONSTRAINTS = list(A = 1, b = 0)))
})

# set_init()
local({
  expect_inherits(set_init(), "latpl_init")
  expect_identical(set_init()$METHOD, "standard")
  expect_error(set_init(METHOD = "custom"))
  expect_error(set_init(THETA = c(1, Inf)))
  expect_error(set_init(CONTROL = set_ucminf_control()))
})

# fit_latpl() checks and stores args, and compute_var() recycle
local({
  d <- set_latpl_data(Y)
  m <- set_latpl_model(matrix(NA, 6, 1))
  expect_error(fit_latpl(Y, m), "inherits")
  expect_error(
    fit_latpl(d, set_latpl_model(matrix(NA, 5, 1))),
    "DATA\\$p == MODEL\\$p"
  )
  expect_error(fit_latpl(d, m, INIT = set_init(THETA = rep(0, 3))), "length")

  fit <- fit_latpl(d, m)
  expect_identical(fit$data, d)
  expect_identical(fit$model, m)
  expect_inherits(fit$control, "latpl_ucminf_control")
  expect_identical(fit$method, "ucminf")

  v <- compute_var(fit)
  expect_equal(dim(v$H), c(fit$dims$npar, fit$dims$npar))
  expect_true(all(is.finite(v$asymptotic_variance)))
})

# pairs_freq() skip missing responses
local({
  set.seed(123)
  p <- 10
  n <- 10
  y <- matrix(as.numeric(rbinom(n * p, 1, .5)), n, p)
  f <- pairs_freq(y, rep(2, p))
  expect_equal(length(f), 4 * p * (p - 1) / 2)
  expect_equal(sum(f), n * p * (p - 1) / 2)
  y[1, 1] <- NA
  expect_equal(sum(pairs_freq(y, rep(2, p))), n * p * (p - 1) / 2 - (p - 1))
  expect_equal(
    sum(pairs_freq(y, rep(2, p), WEIGHTS = rep(2, n))),
    2 * (n * p * (p - 1) / 2 - (p - 1))
  )
})

# get_par() returns pars
local({
  fit <- fit_latpl(
    set_latpl_data(Y),
    set_latpl_model(build_constr_mat(6, 2, "simple"), STDLV = FALSE)
  )
  expect_identical(get_par(fit, OPTION = "raw"), fit$theta)
  par <- get_par(fit)
  expect_identical(
    names(par),
    c(
      "thresholds",
      "loadings",
      "latent_correlations",
      "regressions",
      "residual_covariances"
    )
  )
  expect_equal(par, theta2par(fit$model, fit$theta))
  tr <- get_par(fit, OPTION = "transformed")
  expect_equal(tr[length(tr) - 2], cov2cor(par$latent_correlations)[1, 2])
  expect_equal(tail(tr, 2), sqrt(diag(par$latent_correlations)))
  expect_stdout(print(fit), "Free parameters")
})

# build_constr_mat() output
local({
  p <- 10
  q <- 2
  A <- is.na(build_constr_mat(P = p, Q = q, STRUCT = 'triangular'))
  expect_equal(sum(A), (p * q) - (q * (q - 1) / 2))
  expect_equal(sum(A[lower.tri(A, diag = T)]), (p * q) - (q * (q - 1) / 2))
  expect_equal(sum(A[upper.tri(A, diag = F)]), 0)
  expect_equal(sum(A %in% c(0, 1)), p * q)
  A <- is.na(build_constr_mat(P = p, Q = q, STRUCT = 'simple'))
  expect_equal(sum(A), p)
  expect_equal(sum(rowSums(A) == 1), p)
})


# gen_loadings() errors
local({
  p <- 10
  q <- 2
  A <- build_constr_mat(P = p, Q = q, STRUCT = 'triangular')
  expect_error(gen_loadings(CONSTRMAT = 2))
  expect_error(gen_loadings(CONSTRMAT = NA))
})
