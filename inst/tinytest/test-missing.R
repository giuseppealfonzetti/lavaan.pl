cm <- build_constrMat(8, 2, "simple")
S <- matrix(.4, 2, 2); diag(S) <- 1
Yc <- sim_data(500, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
Y <- Yc
set.seed(4)
for (j in c(2:4, 6:8)) Y[runif(500) < c(.1, .25, .25, .6)[Y[, if (j < 5) 1 else 5] + 1], j] <- NA
M <- set_plfa_model(cm)

# Without missing values all strategies give the complete-data fit
f_l <- fit_plFA(set_plfa_data(Yc), M)
expect_equal(fit_plFA(set_plfa_data(Yc, MISSING = "pairwise"), M)$theta, f_l$theta)
expect_equal(fit_plFA(set_plfa_data(Yc, MISSING = "available.cases"), M)$theta, f_l$theta)

# Listwise deletion
expect_message(d_l <- set_plfa_data(Y), "incomplete observations removed")
expect_equal(d_l$n, sum(stats::complete.cases(Y)))

# Available-cases univariate frequencies weight each observed response by the missing items of its unit
D <- set_plfa_data(Y, MISSING = "available.cases")
m <- rowSums(is.na(Y))
expect_equal(sum(D$ufreq), sum(m * (ncol(Y) - m)))

# Gradient, H and casewise scores of the complete-pairs and available-cases objectives
for (mi in c("pairwise", "available.cases")) {
  D <- set_plfa_data(Y, MISSING = mi)
  f <- fit_plFA(D, M)
  dims <- f$dims
  obj <- function(theta, gr) {
    out <- lavaan.pl:::cpp_multiThread_completePairwise(
      N = dims$n, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT, CONSTRLOGSD = M$CONSTRLOGSD,
      THETA = theta, FREQ = D$freq, CORRFLAG = M$CORRFLAG, NTHR = dims$nthr, NLOAD = dims$nload,
      NCORR = dims$ncorr, NVAR = dims$nvar, GRFLAG = gr, UFREQ = D$ufreq
    )
    if (gr) out$iter_ngradient else out$iter_nll
  }
  th <- f$theta + 0.05
  hj <- lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = th, FREQ = D$freq, DATA = D$Y, C_VEC = dims$cat, CONSTRMAT = M$CONSTRMAT,
    CONSTRLOGSD = M$CONSTRLOGSD, N = dims$n, CORRFLAG = M$CORRFLAG, NTHR = dims$nthr,
    NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar, UFREQ = D$ufreq, UW = D$uweights
  )
  expect_equal(obj(th, 1), numDeriv::grad(function(x) obj(x, 0), th), tolerance = 1e-4)
  expect_equal(-hj$gradH, obj(th, 1))
  expect_equal(hj$gradJ, hj$gradH)
  expect_equal(hj$H, numDeriv::jacobian(function(x) obj(x, 1), th) / dims$n, tolerance = 1e-5)
}
