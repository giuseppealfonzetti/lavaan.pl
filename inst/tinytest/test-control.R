D <- set_plfa_data(as.matrix(data.frame(lapply(LSAT, as.numeric))) - 1)
cl <- set_plfa_model(matrix(NA, 5, 1))

# control constructors return classed lists with the documented defaults
local({
  uc <- set_ucminf_control()
  sa <- set_sa_control()
  init <- set_sa_init_control()
  expect_inherits(uc, "plfa_ucminf_control")
  expect_inherits(sa, "plfa_sa_control")
  expect_inherits(init, "plfa_sa_control")
  expect_equal(unclass(uc), list(GRTOL = 1e-6, XTOL = 1e-12, STEPMAX = 1, MAXEVAL = 500L))
  expect_identical(sa$MAXE, 30L)
  expect_identical(sa$STEP0, 1)
  expect_identical(init$MAXE, 5L)
  expect_null(init$STEP0)
  expect_identical(names(sa), names(init))
})

# control constructors validate their arguments
local({
  expect_error(set_ucminf_control(MAXEVAL = 0))
  expect_error(set_ucminf_control(GRTOL = -1))
  expect_error(set_sa_control(BURNE = -1))
  expect_error(set_sa_control(BURNE = 30, MAXE = 30))
  expect_error(set_sa_control(CHECK_TOL = 1))
  expect_error(set_sa_init_control(STEP0 = 0))
})

# fit_plFA() checks and honours its control objects
local({
  expect_error(fit_plFA(D, cl, CONTROL = list(MAXEVAL = 5)), "inherits")
  expect_warning(fit <- fit_plFA(D, cl, CONTROL = set_sa_control(MAXE = 1, BURNE = 0, PAIRS_PER_ITERATION = 2)), "at least 40 items")
  expect_identical(fit$method, "SA")
  fit <- fit_plFA(D, cl, CONTROL = set_ucminf_control(MAXEVAL = 5))
  expect_true(fit$numFit$info[["neval"]] <= 5)
})

# cfa() checks and honours its control objects
local({
  mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
  expect_error(cfa(mod, LSAT, std.lv = TRUE, control = list(ctrl = list(maxeval = 5))), "set_ucminf_control")
  expect_error(cfa(mod, LSAT, std.lv = TRUE, estimator.args = list(method = "SA")), "set_sa_control")
  expect_error(cfa(mod, LSAT, std.lv = TRUE, estimator.args = list(cpp_control_init = list(MAXE = 2))), "set_init")
  expect_error(cfa(mod, LSAT, std.lv = TRUE, estimator.args = list(init_method = "standard")), "set_init")
  expect_error(cfa(mod, LSAT, std.lv = TRUE, estimator.args = list(init = set_init(THETA = rep(0, 10)))), "start")
  fit <- cfa(mod, LSAT, std.lv = TRUE, se = "none", control = set_ucminf_control(MAXEVAL = 5))
  expect_true(fit@external$lavaan.pl$numFit$info[["neval"]] <= 5)
})
