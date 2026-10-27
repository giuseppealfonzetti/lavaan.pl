if (!at_home()) exit_file("Not run on CRAN")

p <- 40; q <- 4; n <- 1000
cm <- build_constr_mat(p, q, "simple")
S <- matrix(.3, q, q); diag(S) <- 1
Y <- sim_data(n, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1.2, 0, 1.2), S, SEED = 2)
D <- set_latpl_data(Y)
M <- set_latpl_model(cm)

pml <- fit_latpl(D, M)
sa <- fit_latpl(D, M, CONTROL = set_sa_control())

# SA estimates close to PML
expect_equal(sa$stoFit$convergence$convergence_full, 1)
expect_true(max(abs(get_par(sa, OPTION = "transformed") - get_par(pml, OPTION = "transformed"))) < 0.02)

# Noise term uses the number of averaged iterations
vs <- compute_var(sa)
it <- sa$stoFit$last_iter - sa$stoFit$burnt
P <- sa$dims$pairs; ppi <- sa$control$PAIRS_PER_ITERATION
a1 <- P * (P - ppi) / (ppi * (P - 1))
a2 <- (P - ppi) / (ppi * (P - 1))
invH <- solve(vs$H)
noise <- invH %*% (a1 * vs$H - a2 * vs$J) %*% invH / (n * it)
eig <- eigen((noise + t(noise)) / 2, symmetric = TRUE)
expect_equal(vs$noise, eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors)))
expect_equal(vs$optimisation_noise, diag(vs$trJacob %*% vs$noise %*% t(vs$trJacob)))
expect_null(compute_var(pml)$noise)

# Fixed budget: exactly MAXE epochs, averaging after BURNE
fb <- fit_latpl(D, M, CONTROL = set_sa_control(CHECK_TOL = FALSE, BURNE = 2, MAXE = 4))
maxtpe <- P %/% ppi
expect_equal(fb$stoFit$last_iter, 4 * maxtpe + 1)
expect_equal(fb$stoFit$burnt, 3 * maxtpe + 1)
expect_equal(fb$stoFit$convergence$convergence_full, 0)

# cfa() with SA: noise term added to the standard errors
colnames(Y) <- paste0("y", 1:p)
dat <- as.data.frame(lapply(as.data.frame(Y), ordered))
mod <- paste(sapply(1:q, function(h) paste0("f", h, " =~ ", paste0("y", which(is.na(cm[, h])), collapse = " + "))), collapse = "\n")
fit <- cfa(mod, dat, std.lv = TRUE, control = set_sa_control())
template <- lavaan::cfa(mod, dat, std.lv = TRUE, estimator = "PML", se = "robust.huber.white", information = "observed", test = "none", do.fit = FALSE)
plain <- as_lavaan(fit@external$lavaan.pl, modifyList(fit@external$vars, list(noise = NULL)), template)
expect_true(lavaan::lavInspect(fit, "converged"))
expect_true(all(is.finite(lavaan::parameterEstimates(fit)$se)))
expect_true(all(diag(lavaan::vcov(fit)) >= diag(lavaan::vcov(plain))))
expect_true(any(diag(lavaan::vcov(fit)) > diag(lavaan::vcov(plain))))
expect_equal(attr(fit@vcov$vcov, "E.inv"), attr(plain@vcov$vcov, "E.inv"))
expect_true(lavaan::lavInspect(as_lavaan(fb, NULL, template), "converged"))

# SA with available cases lands close to the PML fit with available cases
Ym <- Y
set.seed(6)
for (h in 1:q) {
  idx <- which(is.na(cm[, h]))
  for (j in idx[-1]) Ym[runif(n) < c(.1, .25, .25, .5)[Y[, idx[1]] + 1], j] <- NA
}
Dm <- set_latpl_data(unname(Ym), MISSING = "available.cases")
pml_m <- fit_latpl(Dm, M)
sa_m <- fit_latpl(Dm, M, CONTROL = set_sa_control())
expect_true(max(abs(get_par(sa_m, OPTION = "transformed") - get_par(pml_m, OPTION = "transformed"))) < 0.02)
