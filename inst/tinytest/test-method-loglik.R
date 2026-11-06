#### logLik() is the pairwise log-likelihood ####
local({
  fit <- pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
  engine_fit <- get_lavaanpl_internal(fit)
  vars <- fit@external$vars
  ll <- logLik(fit)
  expect_inherits(ll, "logLik")
  expect_equal(as.numeric(ll), fit@optim$logl)
  expect_equal(attr(ll, "df"), sum(diag(solve(vars$H, vars$J))))
  expect_equal(AIC(fit), -2 * fit@optim$logl + 2 * attr(ll, "df"))
  expect_equal(BIC(fit), -2 * fit@optim$logl + log(engine_fit$data$n) * attr(ll, "df"))
  expect_equal(logLik(pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE, se = "none")), ll, tolerance = 1e-6)
})
