cm <- build_constrMat(9, 3, "simple")
B <- matrix(0, 3, 3); B[2, 1] <- .5; B[3, 1] <- .3; B[3, 2] <- .4
A <- solve(diag(3) - B)
S <- A %*% diag(c(1, .7, .6)) %*% t(A)
S <- stats::cov2cor(S)
Y <- sim_data(600, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
D <- set_plfa_data(Y)
R <- matrix(0, 3, 3); R[2, 1] <- NA; R[3, 1] <- NA; R[3, 2] <- NA
RC <- matrix(0, 9, 9); RC[1, 2] <- RC[2, 1] <- NA; RC[4, 7] <- RC[7, 4] <- NA; RC[5, 6] <- RC[6, 5] <- 0.1

# Analytic gradient and H with regressions and residual covariances, on the
# orthogonal, LKJ and direct paths
for (code in 0:2) for (stdlv in c(TRUE, FALSE)) {
  M <- set_plfa_model(cm, STDLV = stdlv, CORRFLAG = code > 0, LATCOV = if (code == 2) "direct" else "lkj", REGMAT = R, RESCOV = RC)
  dims <- lavaan.pl:::check_dims(D$Y, M)
  th <- c(rep(c(-1, .3, 1), 9), seq(.4, .7, length.out = dims$nload), if (code == 2) c(.1, .15, .05) else if (code == 1) c(.1, .2, .1),
          rep(-.8, dims$nvar), c(.3, .2, .25), c(.05, .1))
  f <- function(theta, gr) {
    out <- lavaan.pl:::cpp_multiThread_completePairwise(
      N = dims$n, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD, THETA = theta, FREQ = D$freq,
      CORRFLAG = code, NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar, GRFLAG = gr,
      REGMAT = M$REGMAT, RESCOV = M$RESCOV
    )
    if (gr) out$iter_ngradient / dims$n else out$iter_nll / dims$n
  }
  hj <- lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = th, FREQ = D$freq, DATA = D$Y, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD,
    N = dims$n, CORRFLAG = code, NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar,
    REGMAT = M$REGMAT, RESCOV = M$RESCOV
  )
  expect_equal(f(th, 1), numDeriv::grad(function(x) f(x, 0), th), tolerance = 1e-5)
  Hnum <- numDeriv::jacobian(function(x) f(x, 1), th)
  expect_equal(hj$H, (Hnum + t(Hnum)) / 2, tolerance = 1e-4)
}

# A saturated structural model is equivalent to the correlated-factor model
fc <- fit_plFA(D, set_plfa_model(cm, STDLV = FALSE))
fs <- fit_plFA(D, set_plfa_model(cm, STDLV = FALSE, CORRFLAG = FALSE, REGMAT = R))
expect_equal(fs$numFit$value, fc$numFit$value, tolerance = 1e-6)
gs <- getPar(fs, OPTION = "list")
As <- solve(diag(3) - gs$regressions)
expect_equal(As %*% gs$latent_correlations %*% t(As), getPar(fc, OPTION = "list")$latent_correlations, tolerance = 1e-3)
expect_equal(computePLRT(fs)$model$stat, computePLRT(fc)$model$stat, tolerance = 1e-4)
f0 <- fit_plFA(D, set_plfa_model(cm, STDLV = FALSE, CORRFLAG = FALSE, REGMAT = replace(R, 3, 0)))
expect_equal(computePLRT(f0, H1 = fc)$model[c("stat", "df")], computePLRT(f0, H1 = fs)$model[c("stat", "df")], tolerance = 1e-4)

# Fixed residual covariances, names and guards
RF <- matrix(0, 9, 9); RF[1, 2] <- RF[2, 1] <- 0.1; RF[4, 7] <- RF[7, 4] <- NA
Mr <- set_plfa_model(cm, REGMAT = R, RESCOV = RF)
fr <- fit_plFA(D, Mr)
expect_equal(getPar(fr, OPTION = "list")$residual_covariances[2, 1], 0.1)
expect_equal(tail(theta_names(D, Mr), 4), c("beta[2,1]", "beta[3,1]", "beta[3,2]", "theta[4,7]"))
expect_identical(lavaan.pl:::check_init_method(list(p = 20), "SA", NULL, Mr), "standard")
expect_identical(lavaan.pl:::check_init_method(list(p = 20), "SA", NULL, set_plfa_model(cm)), "SA")
expect_error(suppressWarnings(fit_plFA(D, Mr, CONTROL = set_sa_control())), "SA")
expect_error(set_plfa_model(cm, REGMAT = diag(3)), "diagonal")
expect_error(set_plfa_model(cm, RESCOV = replace(RC, 2, 0.2)), "symmetric")
