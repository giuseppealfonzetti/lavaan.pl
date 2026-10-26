cm <- build_constrMat(9, 3, "simple")
S <- matrix(.3, 3, 3); diag(S) <- 1
Y <- sim_data(600, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
D <- set_plfa_data(Y)

# Same unconstrained model: direct and LKJ parametrisations reach the same optimum
for (stdlv in c(TRUE, FALSE)) {
  fl <- fit_plFA(D, set_plfa_model(cm, STDLV = stdlv))
  fd <- fit_plFA(D, set_plfa_model(cm, STDLV = stdlv, LATCOV = "direct"))
  expect_equal(fd$numFit$value, fl$numFit$value, tolerance = 1e-8)
  expect_equal(getPar(fd, OPTION = "transformed"), getPar(fl, OPTION = "transformed"), tolerance = 1e-4)
  vl <- computeVar(fl)
  vd <- computeVar(fd)
  expect_equal(vd$asymptotic_variance, vl$asymptotic_variance, tolerance = 1e-3)
  expect_equal(computePLRT(fd, vd)$model[c("stat", "df", "stat.scaled", "df.scaled")],
               computePLRT(fl, vl)$model[c("stat", "df", "stat.scaled", "df.scaled")], tolerance = 1e-3)
}

# Analytic gradient and H of the direct parametrisation
for (stdlv in c(TRUE, FALSE)) {
  M <- set_plfa_model(cm, STDLV = stdlv, LATCOV = "direct")
  dims <- lavaan.pl:::check_dims(D$Y, M)
  th <- fit_plFA(D, M)$theta + 0.02
  ngr <- function(theta, gr) {
    out <- lavaan.pl:::cpp_multiThread_completePairwise(
      N = dims$n, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD, THETA = theta,
      FREQ = D$freq, CORRFLAG = 2L, NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar, GRFLAG = gr
    )
    if (gr) out$iter_ngradient / dims$n else out$iter_nll / dims$n
  }
  hj <- lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = th, FREQ = D$freq, DATA = D$Y, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD,
    N = dims$n, CORRFLAG = 2L, NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar
  )
  expect_equal(ngr(th, 1), numDeriv::grad(function(x) ngr(x, 0), th), tolerance = 1e-5)
  Hnum <- numDeriv::jacobian(function(x) ngr(x, 1), th)
  expect_equal(hj$H, (Hnum + t(Hnum)) / 2, tolerance = 1e-4)
}

# Names, constraints and guards
Md <- set_plfa_model(cm, LATCOV = "direct")
nm <- theta_names(D, Md)
expect_equal(tail(nm, 3), c("psi[1,2]", "psi[1,3]", "psi[2,3]"))
A <- matrix(replace(numeric(length(nm)), match("psi[1,3]", nm), 1), 1)
f0 <- fit_plFA(D, set_plfa_model(cm, LATCOV = "direct", CONSTRAINTS = list(A = A, b = 0)))
expect_equal(getPar(f0, OPTION = "list")$latent_correlations[1, 3], 0)
expect_error(set_plfa_model(cm, CORRFLAG = FALSE, LATCOV = "direct"), "CORRFLAG")
expect_error(fit_plFA(D, set_plfa_model(cm, CONSTRAINTS = list(A = A, b = 0))), "LATCOV")
expect_error(suppressWarnings(fit_plFA(D, Md, CONTROL = set_sa_control())), "SA")
