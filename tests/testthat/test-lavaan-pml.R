skip_on_cran()

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

  test_that(paste("matches lavaan PML:", nm), {
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

test_that("mean.var.adjusted test matches lavaan PML", {
  for (nm in c("ordinal, std.lv", "equality constraints")) {
    g <- grid[[nm]]
    ours <- cfa(g$model, g$data, std.lv = g$std.lv, test = "mean.var.adjusted")
    ref <- lavaan::cfa(g$model, g$data, std.lv = g$std.lv, estimator = "PML", test = "mean.var.adjusted")
    expect_equal(ours@test[[1]][c("stat", "df")], ref@test[[1]][c("stat", "df")], tolerance = 1e-3)
    expect_equal(ours@test[[2]][c("stat", "df")], ref@test[[2]][c("stat", "df")], tolerance = 1e-3)
  }
})

test_that("unsupported models are rejected", {
  expect_error(cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ 0.3*C", bfi, std.lv = TRUE), "covariances")
  expect_error(cfa("A =~ A1 + A2 + A3\nA1 | -1*t1", bfi, std.lv = TRUE), "threshold")
  expect_error(cfa("A =~ A1 + A2 + A3\nA1 ~~ A2", bfi, std.lv = TRUE), "residual")
  expect_error(cfa("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nC ~ A", bfi, std.lv = TRUE), "regression")
  expect_error(cfa("A =~ A1 + A2 + A3", bfi, std.lv = TRUE, parameterization = "theta"), "parameterization")
  expect_error(cfa("A =~ A1 + A2 + A3", bfi, std.lv = TRUE, group = "gender"), "group")
})
