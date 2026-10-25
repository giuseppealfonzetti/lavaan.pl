Y <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1

# set_plfa_data() validates the data and computes the frequencies once
local({
  d <- set_plfa_data(Y)
  expect_inherits(d, "plfa_data")
  expect_identical(c(d$n, d$p), dim(Y))
  expect_equal(unname(d$cat), rep(6, 6))
  expect_equal(d$freq, pairs_freq(Y, d$cat))
  expect_null(d$valfreq)
  expect_false(is.null(set_plfa_data(Y, VALDATA = Y[1:100, ])$valfreq))
  expect_stdout(print(d), "<plfa_data>")

  Yna <- Y
  Yna[1, 1] <- NA
  expect_message(d_na <- set_plfa_data(Yna), "1 incomplete observations removed")
  expect_equal(d_na$n, nrow(Y) - 1)
  expect_error(set_plfa_data(Yna, MISSING = "complete"))
  expect_error(set_plfa_data(Y, VALDATA = Y[, 1:5]))
  expect_error(set_plfa_data(Y, VALDATA = Y + 1))
})

# set_plfa_model() validates the model and applies the identification
local({
  A <- matrix(0, 6, 2)
  A[1:3, 1] <- NA
  A[4:6, 2] <- NA
  m <- set_plfa_model(A)
  expect_inherits(m, "plfa_model")
  expect_identical(c(m$p, m$q), c(6L, 2L))
  expect_identical(m$CONSTRLOGSD, c(0, 0))
  expect_stdout(print(m), "<plfa_model>")

  marker <- set_plfa_model(A, STDLV = FALSE)
  expect_identical(marker$CONSTRMAT[c(1, 4), ], diag(2) * 1)
  expect_true(all(is.na(marker$CONSTRLOGSD)))

  expect_inherits(set_plfa_model(matrix(NA, 6, 1)), "plfa_model")
  expect_inherits(set_plfa_model(A, LLC = list(list(c(2, 1), c(0.5, 3, 1)))), "plfa_model")
  expect_error(set_plfa_model(A, LLC = list(list(c(7, 1), c(0.5, 3, 1)))))
  expect_error(set_plfa_model(A, LLC = list(list(c(2, 1), c(3, 1)))))
})

# set_init() validates its arguments
local({
  expect_inherits(set_init(), "plfa_init")
  expect_identical(set_init()$METHOD, "standard")
  expect_error(set_init(METHOD = "custom"))
  expect_error(set_init(THETA = c(1, Inf)))
  expect_error(set_init(CONTROL = set_ucminf_control()))
})

# fit_plFA() checks its inputs, stores them, and computeVar() reuses the data
local({
  d <- set_plfa_data(Y)
  m <- set_plfa_model(matrix(NA, 6, 1))
  expect_error(fit_plFA(Y, m), "inherits")
  expect_error(fit_plFA(d, set_plfa_model(matrix(NA, 5, 1))), "DATA\\$p == MODEL\\$p")
  expect_error(fit_plFA(d, m, INIT = set_init(THETA = rep(0, 3))), "length")

  fit <- fit_plFA(d, m)
  expect_identical(fit$data, d)
  expect_identical(fit$model, m)
  expect_inherits(fit$control, "plfa_ucminf_control")
  expect_identical(fit$method, "ucminf")

  v <- computeVar(fit)
  expect_equal(dim(v$H), c(fit$dims$npar, fit$dims$npar))
  expect_true(all(is.finite(v$asymptotic_variance)))
})
