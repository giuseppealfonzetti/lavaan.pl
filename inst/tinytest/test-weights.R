cm <- build_constrMat(6, 2, "simple")
S <- matrix(.3, 2, 2); diag(S) <- 1
Y <- sim_data(300, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
M <- set_plfa_model(cm)
set.seed(2)
w <- sample(1:3, nrow(Y), replace = TRUE)

hj <- function(D, theta) {
  dims <- lavaan.pl:::check_dims(D$Y, M)
  lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = theta, FREQ = D$freq, DATA = D$Y, C_VEC = D$cat, CONSTRMAT = M$CONSTRMAT,
    CONSTRLOGSD = M$CONSTRLOGSD, N = D$n, CORRFLAG = M$CORRFLAG, NTHR = dims$nthr,
    NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar, UFREQ = D$ufreq, UW = D$uweights,
    GW = D$weights, JW = D$jweights, CL = D$cluster
  )
}

# Unit weights reproduce the unweighted fit
f <- fit_plFA(set_plfa_data(Y), M)
expect_equal(fit_plFA(set_plfa_data(Y, WEIGHTS = rep(2, nrow(Y))), M)$theta, f$theta)

# Integer frequency weights reproduce the replicated data
Dw <- set_plfa_data(Y, WEIGHTS = w, WEIGHTS_TYPE = "frequency")
Dr <- set_plfa_data(Y[rep(seq_len(nrow(Y)), w), ])
fw <- fit_plFA(Dw, M)
expect_equal(fw$theta, fit_plFA(Dr, M)$theta, tolerance = 1e-6)
th <- fw$theta + 0.05
hw <- hj(Dw, th)
hr <- hj(Dr, th)
expect_equal(hw$H, hr$H)
expect_equal(hw$J, hr$J)

# Weighted scores sum to the weighted gradient
Dd <- set_plfa_data(Y, WEIGHTS = w)
hd <- hj(Dd, th)
expect_equal(hd$gradJ, hd$gradH)

# One cluster per unit: J of the centred cluster totals
n <- nrow(Y)
for (D in list(set_plfa_data(Y), Dd)) {
  h0 <- hj(D, th)
  D$cluster <- as.numeric(seq_len(n) - 1)
  h1 <- hj(D, th)
  zbar <- h0$gradJ / n
  expect_equal(h1$J, n / (n - 1) * (h0$J - tcrossprod(zbar)))
}
