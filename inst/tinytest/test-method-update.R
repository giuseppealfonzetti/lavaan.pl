#### update() refits from the caller's frame ####
local({
  items <- LSAT[1:500, ]
  syntax <- "eta =~ y1 + y2 + y3 + y4 + y5"
  fit <- pmlcfa(syntax, items, std.lv = TRUE)
  refit <- update(fit, se = "none")
  expect_inherits(refit, "lavaan.pl")
  expect_identical(refit@Options$se, "none")
  expect_identical(refit@Options$std.lv, TRUE)
  expect_equal(lavaan::nobs(refit), 500L)
  expect_equal(length(coef(update(fit, model = "eta =~ y1 + y2 + y3 + y4"))), 8)
  added <- update(fit, add = "y1 ~~ y2")
  expect_equal(coef(added), coef(pmlcfa(paste(syntax, "y1 ~~ y2", sep = "\n"), items, std.lv = TRUE)))
})
