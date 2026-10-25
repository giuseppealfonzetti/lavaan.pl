if (!at_home()) exit_file("Not run on CRAN")

grid <- list(
  "binary, one factor" = list(
    model = "eta =~ y1 + y2 + y3 + y4 + y5", data = LSAT, std.lv = TRUE),
  "ordinal, std.lv" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3", data = bfi, std.lv = TRUE),
  "ordinal, marker" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3", data = bfi, std.lv = FALSE),
  "covariance written first" = list(
    model = "A ~~ C\nA =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3", data = bfi, std.lv = FALSE),
  "orthogonal factors" = list(
    model = "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ 0*C", data = bfi, std.lv = TRUE),
  "fixed loading" = list(
    model = "A =~ A1 + A2 + A3\nC =~ 0.6*C1 + C2 + C3", data = bfi, std.lv = TRUE),
  "cross-loading" = list(
    model = "A =~ A1 + A2 + A3 + C1\nC =~ C1 + C2 + C3", data = bfi, std.lv = TRUE),
  "equality constraints" = list(
    model = "A =~ A1 + a*A2 + b*A3 + c*A4 + A5\nC =~ C1 + d*C2 + d*C3\na == 0.5*b + 0.25*c",
    data = bfi, std.lv = TRUE),
  "defined parameter" = list(
    model = "A =~ A1 + a*A2 + b*A3\nC =~ C1 + C2 + C3\nab := a*b", data = bfi, std.lv = FALSE)
)

for (nm in names(grid)) {
  g <- grid[[nm]]
  ours <- cfa(g$model, g$data, std.lv = g$std.lv)
  ref <- lavaan::cfa(g$model, g$data, std.lv = g$std.lv, estimator = "PML")

  local({
    expect_identical(as.character(class(ours)), "lavaan")
    expect_true(lavaan::lavInspect(ours, "converged"))
    pe <- parameterEstimates(ours)
    pr <- parameterEstimates(ref)
    expect_identical(pe[, 1:3], pr[, 1:3])
    expect_equal(pe$est, pr$est, tolerance = 1e-4)
    expect_equal(pe$se, pr$se, tolerance = 1e-3)
    expect_equal(fitted(ours), fitted(ref), tolerance = 1e-4)
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
  })
}

# mean.var.adjusted test and fit measures match lavaan PML
local({
  fm <- c("chisq", "df", "chisq.scaled", "df.scaled", "pvalue.scaled", "chisq.scaling.factor",
          "baseline.chisq", "baseline.df", "baseline.chisq.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled")
  for (nm in c("binary, one factor", "ordinal, std.lv", "ordinal, marker", "cross-loading", "equality constraints")) {
    g <- grid[[nm]]
    ours <- cfa(g$model, g$data, std.lv = g$std.lv, test = "mean.var.adjusted")
    ref <- lavaan::cfa(g$model, g$data, std.lv = g$std.lv, estimator = "PML", test = "mean.var.adjusted")
    expect_equal(ours@test[[1]][c("stat", "df")], ref@test[[1]][c("stat", "df")], tolerance = 1e-3)
    expect_equal(ours@test[[2]][c("stat", "df", "pvalue", "scaling.factor")],
                 ref@test[[2]][c("stat", "df", "pvalue", "scaling.factor")], tolerance = 1e-3)
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
  }
})

# Nested test matches lavaan's PLRT with the exact restriction
local({
  m1 <- "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3 + C4 + C5"
  f1 <- cfa(m1, bfi, std.lv = TRUE, test = "mean.var.adjusted")
  h0 <- list(
    cfa(m1, bfi, std.lv = TRUE, orthogonal = TRUE, test = "mean.var.adjusted"),
    cfa("A =~ A1 + a*A2 + a*A3 + A4 + A5\nC =~ C1 + b*C2 + b*C3 + C4 + C5", bfi, std.lv = TRUE, test = "mean.var.adjusted")
  )
  for (f0 in h0) {
    ours <- lavTestLRT(f1, f0)
    expect_equal(lavaan.pl::anova(f1, f0), ours)
    ref <- suppressWarnings(lavaan::lavTestLRT(f1, f0))
    m10 <- lavaan:::lav_test_diff_m10(m1 = f1, m0 = f0, test = FALSE)
    iH <- solve(lavaan::lavTech(m10, "information.observed"))
    J <- lavaan::lavTech(m10, "first.order")
    A <- lavaan:::lav_test_diff_a(m1 = f1, m0 = f0, method = "exact")
    M <- (A %*% iH %*% J %*% iH %*% t(A)) %*% solve(A %*% iH %*% t(A))
    plrt <- 2 * (f1@optim$logl - f0@optim$logl)
    expect_equal(ours[2, "Chisq diff"], sum(diag(M)) / sum(diag(M %*% M)) * plrt, tolerance = 1e-4)
    expect_equal(ours[2, "Df diff"], sum(diag(M))^2 / sum(diag(M %*% M)), tolerance = 1e-4)
    cols <- c("Df", "PL_AIC", "PL_BIC", "Chisq")
    expect_equal(unlist(ours[, cols]), unlist(ref[, cols]), tolerance = 1e-6)
  }
})

# unsupported models are rejected
local({
  expect_error(cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ 0.3*C", bfi, std.lv = TRUE), "covariances")
  expect_error(cfa("A =~ A1 + A2 + A3\nA1 | -1*t1", bfi, std.lv = TRUE), "threshold")
  expect_error(cfa("A =~ A1 + A2 + A3\nA1 ~~ A2", bfi, std.lv = TRUE), "residual")
  expect_error(cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nC ~ A", bfi, std.lv = TRUE), "regression")
  expect_error(cfa("A =~ A1 + A2 + A3", bfi, std.lv = TRUE, parameterization = "theta"), "parameterization")
  expect_error(cfa("A =~ A1 + A2 + A3", bfi, std.lv = TRUE, group = "gender"), "group")
})
