if (!at_home()) {
  exit_file("Not run on CRAN")
}

p <- 40
q <- 4
n <- 1000
cm <- build_constr_mat(p, q, "simple")
S <- matrix(.3, q, q)
diag(S) <- 1
Y <- sim_data(
  n,
  gen_loadings(cm, SEED = 1, LB = .4, UB = .8),
  c(-1.2, 0, 1.2),
  S,
  SEED = 2
)
D <- set_latpl_data(Y)
# distance between SA and PML estimates, in PML standard errors
dist <- function(sa, pml) {
  se <- sqrt(compute_var(pml)$asymptotic_variance / pml$dims$n)
  max((abs(
    get_par(sa, OPTION = "transformed") - get_par(pml, OPTION = "transformed")
  ) /
    se)[se > 0])
}

# SA estimates within one standard error of PML
for (stdlv in c(TRUE, FALSE)) {
  M <- set_latpl_model(cm, STDLV = stdlv)
  pml <- fit_latpl(D, M)
  sa <- fit_latpl(D, M, CONTROL = set_sa_control())
  expect_true(sa$stoFit$converged)
  expect_identical(sa$stoFit$gain, "hessian")
  expect_true(dist(sa, pml) < 1)
}
M <- set_latpl_model(cm)
pml <- fit_latpl(D, M)
sd <- fit_latpl(D, M, CONTROL = set_sa_control(GAIN = "diagonal"))
expect_identical(sd$stoFit$gain, "diagonal")
expect_identical(sd$stoFit$newton_epochs, 0L)
expect_true(dist(sd, pml) < 1)

# Reproducible with a fixed seed, different with another seed
sa <- fit_latpl(D, M, CONTROL = set_sa_control())
expect_identical(fit_latpl(D, M, CONTROL = set_sa_control())$theta, sa$theta)
expect_false(identical(
  fit_latpl(D, M, CONTROL = set_sa_control(SEED = 7))$theta,
  sa$theta
))

# Averaged phase lasts EPOCHS epochs, averaging after AVG_FROM of it
E <- sa$stoFit$epoch_length
expect_equal(E, ceiling(sa$dims$pairs / 16))
expect_equal(sa$stoFit$averaged, 20 * E - floor(0.25 * 20 * E))
expect_equal(
  sa$stoFit$iterations,
  (sa$stoFit$adam_epochs + sa$stoFit$newton_epochs + 20) * E
)
expect_equal(tail(sa$stoFit$path$phase, 1), 3)

# Noise term uses the number of averaged iterations
vs <- compute_var(sa)
P <- sa$dims$pairs
ppi <- sa$control$PAIRS_PER_ITERATION
a1 <- P * (P - ppi) / (ppi * (P - 1))
a2 <- (P - ppi) / (ppi * (P - 1))
invH <- solve(vs$H)
noise <- invH %*% (a1 * vs$H - a2 * vs$J) %*% invH / (n * sa$stoFit$averaged)
eig <- eigen((noise + t(noise)) / 2, symmetric = TRUE)
expect_equal(vs$noise, eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors)))
expect_equal(
  vs$optimisation_noise,
  diag(vs$trJacob %*% vs$noise %*% t(vs$trJacob))
)
expect_null(compute_var(pml)$noise)

# Linear constraints hold exactly
nm <- theta_names(D, M)
A <- rbind(
  replace(
    numeric(length(nm)),
    match(c("lambda[2,1]", "lambda[3,1]"), nm),
    c(1, -1)
  ),
  replace(numeric(length(nm)), match("tau[4,2]", nm), 1)
)
Mc <- set_latpl_model(cm, CONSTRAINTS = list(A = A, b = c(0, 0.25)))
sc <- fit_latpl(D, Mc, CONTROL = set_sa_control())
expect_equal(as.vector(A %*% sc$theta), c(0, 0.25))
expect_true(dist(sc, fit_latpl(D, Mc)) < 1)
B <- matrix(0, q, q)
B[2, 1] <- B[3, 2] <- NA
RC <- matrix(0, p, p)
RC[1, 2] <- RC[2, 1] <- NA
for (Ms in list(
  set_latpl_model(cm, LATCOV = "direct"),
  set_latpl_model(cm, CORRFLAG = FALSE, REGMAT = B, RESCOV = RC),
  set_latpl_model(cm, PARAMETERIZATION = "theta")
)) {
  expect_true(
    dist(fit_latpl(D, Ms, CONTROL = set_sa_control()), fit_latpl(D, Ms)) < 1
  )
}

# Starting values outside the parameter space are rejected
expect_error(
  fit_latpl(
    D,
    M,
    CONTROL = set_sa_control(),
    INIT = set_init(THETA = replace(pml$theta, 121, 5))
  ),
  "outside the parameter space"
)

# SA Adam phase as initialisation of ucminf
fi <- fit_latpl(D, M, INIT = set_init(METHOD = "SA"))
expect_equal(fi$theta, pml$theta, tolerance = 1e-4)

# Adam phase from several starts: SA continues from the best, counting the epochs of all
nll <- function(theta) -lavaan.pl:::cpp_pl(theta, sa$engine$model, sa$engine$data)$pl
starts <- list(sa$init, pml$theta)
one <- lavaan.pl:::run_sa(starts[1], D, sa$engine, sa$lincon, set_sa_init_control(), TRUE, FALSE)
two <- lavaan.pl:::run_sa(starts, D, sa$engine, sa$lincon, set_sa_init_control(), TRUE, FALSE)
expect_true(nll(two$theta) <= nll(one$theta))
expect_equal(sum(two$path$phase == 0), 2)
full <- lavaan.pl:::run_sa(starts, D, sa$engine, sa$lincon, set_sa_control(), FALSE, FALSE)
expect_equal(full$iterations, (full$adam_epochs + full$newton_epochs + 20) * E)

# cfa() with SA
colnames(Y) <- paste0("y", 1:p)
dat <- as.data.frame(lapply(as.data.frame(Y), ordered))
mod <- paste(
  sapply(1:q, function(h) {
    paste0("f", h, " =~ ", paste0("y", which(is.na(cm[, h])), collapse = " + "))
  }),
  collapse = "\n"
)
fit <- cfa(mod, dat, std.lv = TRUE, control = set_sa_control())
template <- lavaan::cfa(
  mod,
  dat,
  std.lv = TRUE,
  estimator = "PML",
  se = "robust.huber.white",
  information = "observed",
  test = "none",
  do.fit = FALSE
)
plain <- as_lavaan(
  fit@external$lavaan.pl,
  modifyList(fit@external$vars, list(noise = NULL)),
  template
)
expect_true(lavaan::lavInspect(fit, "converged"))
expect_true(all(is.finite(lavaan::parameterEstimates(fit)$se)))
expect_true(all(diag(lavaan::vcov(fit)) >= diag(lavaan::vcov(plain))))
expect_true(any(diag(lavaan::vcov(fit)) > diag(lavaan::vcov(plain))))
expect_equal(attr(fit@vcov$vcov, "E.inv"), attr(plain@vcov$vcov, "E.inv"))

# SA with available cases
Ym <- unname(Y)
set.seed(6)
for (h in 1:q) {
  idx <- which(is.na(cm[, h]))
  for (j in idx[-1]) {
    Ym[runif(n) < c(.1, .25, .25, .5)[Y[, idx[1]] + 1], j] <- NA
  }
}
Dm <- set_latpl_data(Ym, MISSING = "available.cases")
expect_true(
  dist(fit_latpl(Dm, M, CONTROL = set_sa_control()), fit_latpl(Dm, M)) < 1
)
