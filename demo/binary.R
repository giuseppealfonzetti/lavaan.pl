library(lavaan.pl)
data("LSAT")
mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
fit <- pmlcfa(mod, LSAT, std.lv = TRUE, estimator = "PML", se = "none")
summary(fit)

# summary(update(fit, se = "robust.huber.white"))
# summary(update(fit, se = "robust.huber.white", test = "mean.var.adjusted"))
# summary(update(fit, se = TRUE))
