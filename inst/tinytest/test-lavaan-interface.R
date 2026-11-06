#### Data and models ####
lsat_model <- "eta =~ y1 + y2 + y3 + y4 + y5"
three_factors <- "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nE =~ E1 + E2 + E3"
twelve_items <- c(paste0("A", 1:4), paste0("C", 1:4), paste0("E", 1:4))
bfi_sample <- as.data.frame(bfi)
for (item in twelve_items) {
  bfi_sample[[item]] <- ordered(bfi_sample[[item]])
}
set.seed(1)
bfi_sample <- bfi_sample[complete.cases(bfi_sample[, twelve_items]), ][
  sample(2000, 800),
]
three_factors_four_items <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4\nE =~ E1 + E2 + E3 + E4"

#### pmlcfa() returns a lavaan.pl object for binary and ordinal items ####
local({
  for (fit in list(
    pmlcfa(lsat_model, LSAT, std.lv = TRUE),
    pmlcfa(three_factors, bfi, std.lv = TRUE)
  )) {
    expect_inherits(fit, "lavaan.pl")
    expect_inherits(coef(fit), "lavaan.vector")
    expect_false(is.null(fit@external$lavaanpl_internal))
  }
})

#### Integer items are ordinal by default, and items that look continuous are rejected ####
local({
  integers <- as.data.frame(lapply(LSAT, function(x) as.integer(x) - 1L))
  expect_equal(coef(pmlcfa(lsat_model, integers, std.lv = TRUE)), coef(pmlcfa(lsat_model, LSAT, std.lv = TRUE)))
  continuous <- replace(integers, "y1", list(integers$y1 + 0.5))
  expect_error(pmlcfa(lsat_model, continuous, std.lv = TRUE), "non-integer values or more than 10 categories")
  many <- replace(integers, "y1", list(rep(0:11, length.out = nrow(integers))))
  expect_error(pmlcfa(lsat_model, many, std.lv = TRUE), "y1 would be ordinal")
  expect_inherits(pmlcfa(lsat_model, continuous, std.lv = TRUE, ordered = c("y2", "y3", "y4", "y5")), "lavaan.pl")
})

#### Estimators other than PML are left to lavaan ####
local({
  expect_error(pmlcfa(lsat_model, LSAT, estimator = "DWLS"), "lavaan::cfa")
  expect_error(pmlsem(lsat_model, LSAT, estimator = "DWLS"), "lavaan::sem")
  expect_false(is.null(pmlcfa(lsat_model, LSAT, estimator = "PML", se = "none")@external$lavaanpl_internal))
})

#### std.lv gives the same fit ####
local({
  fit_marker <- pmlcfa(three_factors, bfi, std.lv = FALSE, control = set_qn_control())
  fit_stdlv <- pmlcfa(three_factors, bfi, std.lv = TRUE, control = set_qn_control())
  expect_equal(
    fit_marker@external$lavaanpl_internal$numFit$value,
    fit_stdlv@external$lavaanpl_internal$numFit$value
  )
  expect_equal(fitted(fit_marker), fitted(fit_stdlv), tolerance = 1e-3)
})

#### SA is selected through control ####
local({
  fit <- suppressWarnings(pmlcfa(
    three_factors,
    bfi,
    std.lv = TRUE,
    control = set_sa_control(),
    se = "none"
  ))
  expect_inherits(fit, "lavaan.pl")
  expect_identical(fit@external$lavaanpl_internal$method, "SA")
})

#### Fixed loadings and equality constraints are honoured ####
local({
  fit <- pmlcfa(
    "A =~ -0.4*A1 + A2 + a*A3 + b*A4 + c*A5\nC =~ 0.6*C1 + C2 + C3 + C4 + C5\na == 0.5*b + 0.25*c",
    bfi,
    std.lv = TRUE
  )
  loadings <- lavaan::lavInspect(fit, "est")$lambda
  expect_equal(loadings[1, 1], -0.4)
  expect_equal(loadings[6, 2], 0.6)
  expect_equal(loadings[3, 1], 0.5 * loadings[4, 1] + 0.25 * loadings[5, 1])
})

#### update() refits with lavaan.pl ####
local({
  fit <- pmlcfa(lsat_model, LSAT, std.lv = TRUE)
  fit_no_se <- update(fit, se = "none")
  fit_four_items <- update(fit, model = "eta =~ y1 + y2 + y3 + y4")
  expect_identical(fit_no_se@Options$se, "none")
  expect_false(is.null(fit_no_se@external$lavaanpl_internal))
  expect_equal(length(coef(fit_four_items)), 8)
  expect_false(is.null(fit_four_items@external$lavaanpl_internal))
})

#### Numerical and analytic standard errors agree ####
local({
  fit_analytic <- pmlcfa(lsat_model, LSAT, std.lv = TRUE)
  fit_numeric <- pmlcfa(
    lsat_model,
    LSAT,
    std.lv = TRUE,
    estimator.args = list(computevar_numderiv = TRUE)
  )
  expect_equal(
    attr(fit_analytic@vcov$vcov, "E.inv"),
    attr(fit_numeric@vcov$vcov, "E.inv"),
    tolerance = 1e-5
  )
})

#### pmlcfa() validates the control arguments ####
local({
  expect_error(
    pmlcfa(
      lsat_model,
      LSAT,
      std.lv = TRUE,
      control = list(ctrl = list(maxeval = 5))
    ),
    "set_qn_control"
  )
  expect_error(
    pmlcfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(method = "SA")),
    "optim.method"
  )
  expect_error(
    pmlcfa(
      lsat_model,
      LSAT,
      std.lv = TRUE,
      estimator.args = list(init = 1)
    ),
    "Unsupported `estimator.args`: init"
  )
  fit <- pmlcfa(
    lsat_model,
    LSAT,
    std.lv = TRUE,
    se = "none",
    control = set_qn_control(maxeval = 5)
  )
  expect_true(fit@external$lavaanpl_internal$numFit$info[["neval"]] <= 5)
})

#### Unsupported models are rejected ####
local({
  bfi_weighted <- data.frame(bfi, w = 1)
  expect_error(
    pmlcfa("A =~ A1 + a*A2 + A3\nC =~ C1 + C2 + C3\na > 2", bfi, std.lv = TRUE),
    "Inequality constraints"
  )
  expect_error(
    pmlcfa("A =~ A1 + a*A2 + b*A3\nC =~ C1 + C2 + C3\na > b", bfi, std.lv = TRUE),
    "Inequality constraints"
  )
  expect_error(
    pmlcfa(three_factors, bfi, std.lv = TRUE, bounds = "standard"),
    "bounds"
  )
  expect_error(pmlcfa(
    three_factors,
    bfi,
    std.lv = TRUE,
    missing = "doubly.robust"
  ))
  expect_error(
    pmlcfa(
      three_factors,
      bfi_weighted,
      std.lv = TRUE,
      sampling.weights = "w",
      sampling.weights.normalization = "none"
    ),
    "normalization"
  )
  expect_error(
    pmlcfa(
      "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v*A\nA ~~ c*C\nc == v",
      bfi,
      std.lv = FALSE
    ),
    "covariances and variances"
  )
  expect_error(
    pmlcfa(
      "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v1*A\nC ~~ v2*C\nv1 == v2 + 0.5",
      bfi,
      std.lv = FALSE
    ),
    "variances"
  )
  expect_error(
    pmlcfa("A =~ A1 + a*A2 + A3\nA1 | t*t1\na == t", bfi, std.lv = TRUE),
    "either thresholds"
  )
  expect_error(
    pmlsem("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA1 ~ C", bfi, std.lv = TRUE),
    "observed variables"
  )
})

#### A just-identified higher-order factor reproduces the correlated factors ####
local({
  fit_correlated <- pmlcfa(
    three_factors_four_items,
    bfi_sample,
    std.lv = TRUE,
    control = set_qn_control()
  )
  fit_higher_order <- pmlsem(
    paste(three_factors_four_items, "\nG =~ A + C + E"),
    bfi_sample,
    std.lv = TRUE,
    control = set_qn_control()
  )
  expect_equal(
    fit_higher_order@optim$logl,
    fit_correlated@optim$logl,
    tolerance = 1e-8
  )
})

#### A residual covariance improves the fit in a nested test ####
local({
  fit_correlated <- pmlcfa(
    three_factors_four_items,
    bfi_sample,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  fit_residual <- pmlcfa(
    paste(three_factors_four_items, "\nA1 ~~ A2"),
    bfi_sample,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  difference <- lavTestLRT(fit_residual, fit_correlated)
  expect_true(difference[2, "Df diff"] > 0)
  expect_true(difference[2, "Chisq diff"] > 0)
})

#### optim.method selects the algorithm, and is shown by lavaan ####
local({
  fits <- list(
    `SA-QN` = pmlcfa(lsat_model, LSAT, std.lv = TRUE),
    SA = pmlcfa(lsat_model, LSAT, std.lv = TRUE, optim.method = "SA"),
    QN = pmlcfa(lsat_model, LSAT, std.lv = TRUE, optim.method = "QN")
  )
  for (method in names(fits)) {
    expect_identical(fits[[method]]@external$lavaanpl_internal$method, method)
    expect_identical(fits[[method]]@Options$optim.method, method)
    expect_true(any(grepl(method, capture.output(summary(fits[[method]])), fixed = TRUE)))
  }
  expect_equal(unname(coef(fits$`SA-QN`)), unname(coef(fits$QN)), tolerance = 1e-4)
  expect_identical(
    pmlcfa(lsat_model, LSAT, std.lv = TRUE, optim.method = "SA", control = set_sa_control(average_epochs = 5))@external$lavaanpl_internal$control$average_epochs,
    5L
  )
  expect_identical(
    pmlcfa(lsat_model, LSAT, std.lv = TRUE, control = set_sa_control())@external$lavaanpl_internal$method,
    "SA"
  )
  expect_error(pmlcfa(lsat_model, LSAT, std.lv = TRUE, optim.method = "nlminb"), "must be")
  expect_error(pmlcfa(lsat_model, LSAT, std.lv = TRUE, optim.method = "QN", control = set_sa_control()), "is for SA")
  expect_error(pmlcfa(lsat_model, LSAT, std.lv = TRUE, control = list()), "must be created")
})

#### lavaan's slots come from the engine, not from a second lavaan fit ####
local({
  set.seed(2)
  frame <- as.data.frame(bfi)[1:600, ]
  frame$x1 <- rnorm(600)
  frame$w <- runif(600, .5, 2)
  syntax <- "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nC ~ x1"
  fit <- pmlsem(syntax, frame, std.lv = TRUE, sampling.weights = "w")
  engine_fit <- fit@external$lavaanpl_internal
  logl <- lavaan.pl:::cpp_pl(engine_fit$parameter_vector, engine_fit$engine$model, engine_fit$engine$data)$pl
  expect_identical(fit@optim$logl, logl)
  expect_equal(fit@optim$fx, -logl / nobs(fit))
  expect_equal(unname(vcov(fit)[1, 1]), parameterEstimates(fit)$se[1]^2)
  expect_false(isTRUE(all.equal(
    lavaan::lavInspect(fit, "start")$lambda,
    lavaan::lavInspect(fit, "est")$lambda,
    check.attributes = FALSE
  )))
  fit_ordinal <- pmlcfa(lsat_model, LSAT, std.lv = TRUE)
  expect_equal(fit_ordinal@optim$npar, length(coef(fit_ordinal)))
  expect_true(all(is.finite(parameterEstimates(pmlcfa(lsat_model, LSAT, std.lv = TRUE, se = "none"))$est)))
  expect_warning(
    suppressMessages(pmlcfa("A =~ A1 + A2\nC =~ C1 + C2 + C3", frame, std.lv = TRUE, se = "none")),
    "ov variances are negative"
  )
})

#### lavaan's snake_case options are understood ####
local({
  fit <- pmlcfa(lsat_model, LSAT, std_lv = TRUE, optim_method = "QN")
  expect_true(fit@Options$std.lv)
  expect_identical(fit@external$lavaanpl_internal$method, "QN")
  expect_identical(fit@Options$optim.method, "QN")
  expect_false(fit@Options$ceq.simple)
  expect_error(
    pmlcfa("eta =~ y1 + a*y2 + a*y3 + y4 + y5", LSAT, std.lv = TRUE, ceq.simple = TRUE),
    "ceq.simple"
  )
})

#### samplestats = FALSE gives the same fit ####
local({
  for (group in list(NULL, "gender")) {
    full <- pmlcfa(three_factors, bfi_sample, std.lv = TRUE, group = group)
    fast <- pmlcfa(three_factors, bfi_sample, std.lv = TRUE, group = group, samplestats = FALSE)
    expect_equal(coef(fast), coef(full))
    expect_equal(vcov(fast), vcov(full))
    expect_equal(logLik(fast), logLik(full))
  }
})
