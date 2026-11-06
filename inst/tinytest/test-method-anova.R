bfi_sample <- as.data.frame(bfi)
three_factors <- "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nE =~ E1 + E2 + E3"

#### anova() and lavTestLRT() use lavaan.pl:::compute_plrt() ####
local({
  fit_h0 <- pmlcfa(three_factors, bfi_sample, std.lv = TRUE, test = "mean.var.adjusted")
  fit_h1 <- pmlcfa(paste(three_factors, "\nA1 ~~ C1"), bfi_sample, std.lv = TRUE, test = "mean.var.adjusted")
  ours <- anova(fit_h1, fit_h0)
  expect_true(grepl("lavaan.pl", attr(ours, "heading"), fixed = TRUE))
  expect_equal(lavTestLRT(fit_h1, fit_h0), ours)
  expect_equal(rownames(ours), c("fit_h1", "fit_h0"))
  expect_equal(ours$PL_AIC, c(AIC(fit_h1), AIC(fit_h0)))
  expect_equal(ours$PL_BIC, c(BIC(fit_h1), BIC(fit_h0)))
  expect_equal(lavaan::anova(fit_h1, fit_h0), ours)
  expect_false(grepl("lavaan.pl", attr(lavaan::lavTestLRT(fit_h1, fit_h0), "heading"), fixed = TRUE))
})
