cm <- build_constrMat(6, 2, "simple")
S <- matrix(.3, 2, 2); diag(S) <- 1
Y <- sim_data(400, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
D <- set_plfa_data(Y)
M0 <- set_plfa_model(cm)
nm <- theta_names(D, M0)
d <- length(nm)
f0 <- fit_plFA(D, M0)

# Parameter names follow the parameter vector
expect_equal(length(nm), length(f0$theta))
expect_equal(nm[c(1, 19, 24, 25)], c("tau[1,1]", "lambda[1,1]", "lambda[6,2]", "eta[1,2]"))
expect_equal(tail(theta_names(D, set_plfa_model(cm, STDLV = FALSE)), 3), c("eta[1,2]", "logsd[1]", "logsd[2]"))
expect_equal(f0$dims$npar, d)

# Equal loadings and a fixed threshold, through the developer interface
A <- rbind(replace(numeric(d), match(c("lambda[2,1]", "lambda[3,1]"), nm), c(1, -1)),
           replace(numeric(d), match("tau[4,2]", nm), 1))
M <- set_plfa_model(cm, CONSTRAINTS = list(A = A, b = c(0, 0.25)))
f <- fit_plFA(D, M)
expect_equal(f$dims$npar, d - 2)
expect_equal(as.vector(A %*% f$theta), c(0, 0.25))
dims <- f$dims
ngr <- function(theta) {
  lavaan.pl:::cpp_multiThread_completePairwise(
    N = dims$n, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD, THETA = theta,
    FREQ = D$freq, CORRFLAG = M$CORRFLAG, NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr,
    NVAR = dims$nvar, GRFLAG = 1
  )$iter_ngradient / dims$n
}
K <- f$lincon$K
expect_true(max(abs(crossprod(K, ngr(f$theta)))) < 1e-5)

# Free-coordinate H against the numerical Jacobian of the constrained gradient
v <- computeVar(f)
Hnum <- numDeriv::jacobian(function(x) as.vector(crossprod(K, ngr(as.vector(K %*% x + f$lincon$c)))), f$theta[f$lincon$free_idx])
expect_equal(v$H, (Hnum + t(Hnum)) / 2, tolerance = 1e-4)

# Unsupported constraints
expect_error(fit_plFA(D, set_plfa_model(cm, CONSTRAINTS = list(A = matrix(replace(numeric(d), match("eta[1,2]", nm), 1), 1), b = 0))), "correlations")
expect_error(fit_plFA(D, set_plfa_model(cm, CONSTRAINTS = list(A = matrix(replace(numeric(d), c(1, 19), c(1, -1)), 1), b = 0))), "either thresholds")
expect_error(fit_plFA(D, set_plfa_model(cm, CONSTRAINTS = list(A = rbind(A[1, ], A[1, ]), b = c(0, 1)))), "inconsistent")
expect_error(fit_plFA(D, M, CONTROL = set_sa_control()), "SA")
