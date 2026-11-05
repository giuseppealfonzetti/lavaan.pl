#### parameterEstimates() and standardizedSolution() keep the standard errors ####
local({
  fit <- pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
  plain <- as(fit, "lavaan")
  expect_equal(parameterEstimates(fit), lavaan::parameterEstimates(plain))
  expect_equal(parameterestimates(fit), lavaan::parameterEstimates(plain))
  expect_equal(lavParameterEstimates(fit), lavaan::parameterEstimates(plain))
  expect_true(all(is.finite(parameterEstimates(fit)$se)))
  expect_equal(standardizedSolution(fit), lavaan::standardizedSolution(plain))
  expect_equal(standardizedsolution(fit), lavaan::standardizedSolution(plain))
  expect_null(lavaan::parameterEstimates(fit)$se)
  lavaan_fit <- lavaan::cfa("eta =~ y1 + y2 + y3", LSAT)
  expect_equal(parameterEstimates(lavaan_fit), lavaan::parameterEstimates(lavaan_fit))
})
