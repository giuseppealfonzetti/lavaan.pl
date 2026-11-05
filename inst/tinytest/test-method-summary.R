#### summary() keeps lavaan's inference ####
local({
  fit <- pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
  out <- capture.output(summary(fit, standardized = TRUE))
  expect_true(startsWith(out[1], "lavaan.pl"))
  expect_true(any(grepl("Std.Err", out, fixed = TRUE)))
  expect_equal(out[-1], capture.output(lavaan::summary(as(fit, "lavaan"), standardized = TRUE)))
})
