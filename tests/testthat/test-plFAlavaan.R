test_that("[binary] cfa function works", {

  # Generate some data
  set.seed(456)
  mod <- "
  f1 =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5
  f2 =~ 0.7*y6 + 0.7*y7 + 0.7*y8 + 0.7*y9 + 0.7*y10
  f1 ~~ 0.3*f2
  f1 ~~ 1*f1
  f2 ~~ 1*f2

  y1 | 0.5*t1
  y2 | 0.6*t1
  y3 | 0.7*t1
  y4 | 0.8*t1
  y5 | 0.9*t1
  y6 | 0.5*t1
  y7 | 0.6*t1
  y8 | 0.7*t1
  y9 | 0.8*t1
  y10 | 0.9*t1
  "
  dat <- lavaan::simulateData(mod, sample.nobs = 2000)
  dat <- as.data.frame(lapply(dat, ordered))

  # lavaan syntax fit
  mod <- "
    eta1 =~ y1 + y2 + y3 + y4 + y5
    eta2 =~ y6 + y7 + y8 + y9 + y10
  "
  fit <- cfa(mod, dat, std.lv = TRUE)

  expect_s4_class(fit, "lavaan")
  expect_s4_class(fit, "plFAlavaan")
  expect_true(inherits(coef(fit), "lavaan.vector"))
})

test_that("Revert to lavaan works", {
  fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, estimator = "DWLS")
  expect_false(inherits(fit, "plFAlavaan"))
})

test_that("[ordinal] cfa function works", {

  # Generate some data
  set.seed(456)
  mod <- "
  f1 =~ 0.7*y1 + 0.7*y2 + 0.7*y3
  f2 =~ 0.7*y4 + 0.7*y5 + 0.7*y6
  f3 =~ 0.7*y7 + 0.7*y8 + 0.7*y9
  f1 ~~ 0.3*f2
  f1 ~~ 0.2*f3
  f2 ~~ 0.4*f3
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f3 ~~ 1*f3

  y1 | -0.5*t1 + 0*t2 + 0.5*t3
  y2 | -0.5*t1 + 0*t2 + 0.5*t3
  y3 | -0.5*t1 + 0*t2 + 0.5*t3
  y4 | -0.5*t1 + 0*t2 + 0.5*t3
  y5 | -0.5*t1 + 0*t2 + 0.5*t3
  y6 | -0.5*t1 + 0*t2 + 0.5*t3
  y7 | -0.5*t1 + 0*t2 + 0.5*t3
  y8 | -0.5*t1 + 0*t2 + 0.5*t3
  y9 | -0.5*t1 + 0*t2 + 0.5*t3
  "
  dat <- lavaan::simulateData(mod, sample.nobs = 2000)
  dat <- as.data.frame(lapply(dat, ordered))

  # lavaan syntax fit
  mod <- "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
    eta3 =~ y7 + y8 + y9
  "
  fit <- cfa(mod, dat, std.lv = TRUE)

  expect_s4_class(fit, "lavaan")
  expect_s4_class(fit, "plFAlavaan")
  expect_true(inherits(coef(fit), "lavaan.vector"))
})

test_that("`std.lv` argument works", {
  mod <- "
    A =~ A1 + A2 + A3
    C =~ C1 + C2 + C3
    E =~ E1 + E2 + E3
  "
  fit1 <- cfa(mod, bfi, std.lv = FALSE)
  fit2 <- cfa(mod, bfi, std.lv = TRUE)
  expect_equal(
    fit1@external$lavaan.pl@numFit$value,
    fit2@external$lavaan.pl@numFit$value
  )
  expect_equal(fitted(fit1), fitted(fit2), tolerance = 1e-3)
})

test_that("Method = 'SA' works", {
  mod <- "
    A =~ A1 + A2 + A3 + A4 + A5
    C =~ C1 + C2 + C3 + C4 + C5
    E =~ E1 + E2 + E3 + E4 + E5
    N =~ N1 + N2 + N3 + N4 + N5
    O =~ O1 + O2 + O3 + O4 + O5
  "

  fit <- cfa(mod, bfi, std.lv = TRUE, estimator.args = list(method = "SA"),
             se = "none")

  expect_s4_class(fit, "lavaan")
  expect_s4_class(fit, "plFAlavaan")
  expect_true(inherits(coef(fit), "lavaan.vector"))
})

test_that("Constraints work", {
  mod <- "
    A =~ -0.4*A1 + A2 + a*A3 + b*A4 + c*A5
    C =~ 0.6*C1 + C2 + C3 + C4 + C5
    A ~~ -0.4*C
    a == 0.5*b + 0.25*c
  "
  fit <- cfa(mod, bfi, std.lv = TRUE)
  Lambda <- lavaan::inspect(fit, "est")$lambda
  Psi <- lavaan::inspect(fit, "est")$psi
  expect_equal(Lambda[1, 1], -0.4)
  expect_equal(Lambda[6, 2], 0.6)
  expect_equal(Psi[1, 2], -0.4)

})

# Test for starting values
#
# fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT)
# fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, start = coef(fit))

test_that("Against lavaan::cfa", {
  mod <- "eta1 =~ y1 + y2 + y3 + y4 + y5"

  fit1 <- cfa(mod, LSAT, std.lv = TRUE)
  # fit2 <- cfa(mod, bfi, std.lv = TRUE, estimator.args = list(method = "SA"))
  fit3 <- lavaan::cfa(mod, LSAT, std.lv = TRUE, estimator = "WLSMV")
  fit4 <- lavaan::cfa(mod, LSAT, std.lv = TRUE, estimator = "PML")

  expect_equal(coef(fit1), coef(fit3), tolerance = 1e-2)
  expect_equal(coef(fit3), coef(fit4), tolerance = 1e-3)

})


test_that("Hinv is the same", {

  information_matrix <- function(fit, type) {
    lavargs <- list(
      lavmodel = fit@Model,
      lavsamplestats = fit@SampleStats,
      lavdata = fit@Data,
      lavoptions = fit@Options,
      lavcache = fit@Cache
    )

    if (type == "observed")
      out <- do.call("lav_model_info_observed", lavargs,
                     envir = asNamespace("lavaan"))
    else if (type == "expected")
      out <- do.call("lav_model_info_expected", lavargs,
                     envir = asNamespace("lavaan"))
    else stop("Invalid type")
    return(out)
  }

  # From lavaan.pl
  mod <- "eta1 =~ y1 + y2 + y3 + y4 + y5"

  fit <- cfa(mod, LSAT, std.lv = TRUE)
  Hinv1 <- with(fit@external, computeVar(lavaan.pl, D)$invH[idx_plFA2lav, idx_plFA2lav])

  # Getting Hinv via numDeriv
  fit <- cfa(mod, LSAT, std.lv = TRUE,
             estimator.args = list(computevar_numderiv = TRUE))
  Hinv2 <- with(fit@external, vars$invH[idx_plFA2lav, idx_plFA2lav])

  # Gtting Hinv from lavaan, which also uses numDeriv
  Hinv3 <- solve(information_matrix(fit, "observed"))

  expect_equal(Hinv1, Hinv2, tolerance = 1e-5)
  expect_equal(Hinv1, Hinv3, tolerance = 1e-5)

})

test_that("Jmat is the same", {
  mod <- "eta1 =~ y1 + y2 + y3 + y4 + y5"
  fit1 <- cfa(mod, LSAT, std.lv = TRUE)
  fit2 <- lavaan::cfa(mod, LSAT, std.lv = TRUE, estimator = "PML")

  J1 <- with(fit1@external, computeVar(lavaan.pl, D)$J[idx_plFA2lav, idx_plFA2lav])
  J2 <- lavaan:::lav_model_info_firstorder(
    lavmodel =fit2@Model,
    lavsamplestats = fit2@SampleStats,
    lavdata = fit2@Data,
    lavoptions = fit2@Options,
    lavcache = fit2@Cache
  )
  expect_equal(J1, J2, tolerance = 1e-5)
})




