#### Data and models ####
lsat_model <- "eta =~ y1 + y2 + y3 + y4 + y5"
three_factors <- "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nE =~ E1 + E2 + E3"
twelve_items <- c(paste0("A", 1:4), paste0("C", 1:4), paste0("E", 1:4))
bfi_sample <- as.data.frame(bfi)
for (item in twelve_items) {
  bfi_sample[[item]] <- ordered(bfi_sample[[item]])
}
set.seed(1)
bfi_sample <- bfi_sample[complete.cases(bfi_sample[, twelve_items]), ][sample(2000, 800), ]
three_factors_four_items <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4\nE =~ E1 + E2 + E3 + E4"

#### cfa() returns a lavaan object for binary and ordinal items ####
local({
  for (fit in list(cfa(lsat_model, LSAT, std.lv = TRUE), cfa(three_factors, bfi, std.lv = TRUE))) {
    expect_identical(as.character(class(fit)), "lavaan")
    expect_inherits(coef(fit), "lavaan.vector")
    expect_false(is.null(fit@external$lavaan.pl))
  }
})

#### Estimators other than PML are fitted by lavaan ####
local({
  expect_null(cfa(lsat_model, LSAT, estimator = "DWLS")@external$lavaan.pl)
})

#### std.lv gives the same fit ####
local({
  fit_marker <- cfa(three_factors, bfi, std.lv = FALSE)
  fit_stdlv <- cfa(three_factors, bfi, std.lv = TRUE)
  expect_equal(fit_marker@external$lavaan.pl$numFit$value, fit_stdlv@external$lavaan.pl$numFit$value)
  expect_equal(fitted(fit_marker), fitted(fit_stdlv), tolerance = 1e-3)
})

#### SA is selected through control ####
local({
  fit <- suppressWarnings(cfa(three_factors, bfi, std.lv = TRUE, control = set_sa_control(), se = "none"))
  expect_identical(as.character(class(fit)), "lavaan")
  expect_identical(fit@external$lavaan.pl$method, "SA")
})

#### Fixed loadings and equality constraints are honoured ####
local({
  fit <- cfa(
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
  fit <- cfa(lsat_model, LSAT, std.lv = TRUE)
  fit_no_se <- update(fit, se = "none")
  fit_four_items <- update(fit, model = "eta =~ y1 + y2 + y3 + y4")
  expect_identical(fit_no_se@Options$se, "none")
  expect_false(is.null(fit_no_se@external$lavaan.pl))
  expect_equal(length(coef(fit_four_items)), 8)
  expect_false(is.null(fit_four_items@external$lavaan.pl))
})

#### Numerical and analytic standard errors agree ####
local({
  fit_analytic <- cfa(lsat_model, LSAT, std.lv = TRUE)
  fit_numeric <- cfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(computevar_numderiv = TRUE))
  expect_equal(attr(fit_analytic@vcov$vcov, "E.inv"), attr(fit_numeric@vcov$vcov, "E.inv"), tolerance = 1e-5)
})

#### cfa() validates the control arguments ####
local({
  expect_error(cfa(lsat_model, LSAT, std.lv = TRUE, control = list(ctrl = list(maxeval = 5))), "set_ucminf_control")
  expect_error(cfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(method = "SA")), "set_sa_control")
  expect_error(
    cfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(cpp_control_init = list(MAXE = 2))),
    "set_init"
  )
  expect_error(cfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(init_method = "standard")), "set_init")
  expect_error(
    cfa(lsat_model, LSAT, std.lv = TRUE, estimator.args = list(init = set_init(THETA = rep(0, 10)))),
    "start"
  )
  fit <- cfa(lsat_model, LSAT, std.lv = TRUE, se = "none", control = set_ucminf_control(MAXEVAL = 5))
  expect_true(fit@external$lavaan.pl$numFit$info[["neval"]] <= 5)
})

#### Unsupported models are rejected ####
local({
  bfi_weighted <- data.frame(bfi, w = 1)
  expect_error(cfa("A =~ A1 + a*A2 + A3\nC =~ C1 + C2 + C3\na > 2", bfi, std.lv = TRUE), "Inequality constraints")
  expect_error(cfa("A =~ A1 + a*A2 + b*A3\nC =~ C1 + C2 + C3\na > b", bfi, std.lv = TRUE), "Inequality constraints")
  expect_error(cfa(three_factors, bfi, std.lv = TRUE, bounds = "standard"), "bounds")
  expect_error(cfa(three_factors, bfi, std.lv = TRUE, missing = "doubly.robust"))
  expect_error(
    cfa(three_factors, bfi_weighted, std.lv = TRUE, sampling.weights = "w", sampling.weights.normalization = "none"),
    "normalization"
  )
  expect_error(
    cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v*A\nA ~~ c*C\nc == v", bfi, std.lv = FALSE),
    "covariances and variances"
  )
  expect_error(
    cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v1*A\nC ~~ v2*C\nv1 == v2 + 0.5", bfi, std.lv = FALSE),
    "variances"
  )
  expect_error(cfa("A =~ A1 + a*A2 + A3\nA1 | t*t1\na == t", bfi, std.lv = TRUE), "either thresholds")
  expect_error(sem("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA1 ~ C", bfi, std.lv = TRUE), "observed variables")
})

#### A just-identified higher-order factor reproduces the correlated factors ####
local({
  fit_correlated <- cfa(three_factors_four_items, bfi_sample, std.lv = TRUE)
  fit_higher_order <- sem(paste(three_factors_four_items, "\nG =~ A + C + E"), bfi_sample, std.lv = TRUE)
  expect_equal(fit_higher_order@optim$logl, fit_correlated@optim$logl, tolerance = 1e-8)
})

#### A residual covariance improves the fit in a nested test ####
local({
  fit_correlated <- cfa(three_factors_four_items, bfi_sample, std.lv = TRUE, test = "mean.var.adjusted")
  fit_residual <- cfa(paste(three_factors_four_items, "\nA1 ~~ A2"), bfi_sample, std.lv = TRUE,
                      test = "mean.var.adjusted")
  difference <- lavTestLRT(fit_residual, fit_correlated)
  expect_true(difference[2, "Df diff"] > 0)
  expect_true(difference[2, "Chisq diff"] > 0)
})
