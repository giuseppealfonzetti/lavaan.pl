cm <- build_constr_mat(9, 3, "simple")
cm[4, 1] <- NA
S <- matrix(.3, 3, 3); diag(S) <- 1
Y <- sim_data(300, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
Y[, 1:3] <- 1 * (Y[, 1:3] >= 2)
Ym <- Y
set.seed(4)
Ym[sample(length(Ym), 150)] <- NA
D <- set_latpl_data(Ym, MISSING = "available.cases")
n <- D$n
B <- matrix(0, 3, 3); B[2, 1] <- B[3, 1] <- B[3, 2] <- NA
RC <- matrix(0, 9, 9); RC[1, 2] <- RC[2, 1] <- NA; RC[5, 6] <- RC[6, 5] <- 0.1

# Casewise scores of each unit, from the engine on its own response pattern
unit_scores <- function(f, theta) {
  t(sapply(seq_len(n), function(i) {
    y <- D$Y[i, , drop = FALSE]
    um <- D$uweights[i]
    di <- f$engine$data
    di$FREQ <- pairs_freq(y, D$cat)
    di$UFREQ <- unlist(lapply(seq_len(D$p), function(j) um * (seq_len(D$cat[j]) - 1 == y[j] & !is.na(y[j]))))
    lavaan.pl:::cpp_pl(theta, f$engine$model, di, 1L)$gradient
  }))
}

# Gradient, H, J, scores and Jacobian of rho in every latent form, identification and structure,
# with binary and ordinal items and available-cases missing values
for (form in c("orthogonal", "lkj", "direct")) for (stdlv in c(TRUE, FALSE)) for (sem in c(FALSE, TRUE)) {
  M <- set_latpl_model(cm, STDLV = stdlv, CORRFLAG = form != "orthogonal", LATCOV = if (form == "direct") "direct" else "lkj",
                      REGMAT = if (sem) B, RESCOV = if (sem) RC)
  f <- suppressWarnings(fit_latpl(D, M, CONTROL = set_ucminf_control(MAXEVAL = 1)))
  dims <- f$dims
  th <- c(unlist(lapply(D$cat - 1, function(m) seq(-0.8, 0.8, length.out = m))), seq(.4, .6, length.out = dims$nload),
          if (form == "direct") c(.1, .15, .05) else if (form == "lkj") c(.1, .2, .1), rep(-.5, dims$nvar),
          if (sem) c(.3, .2, .25, .05))
  pl <- function(x, o) {
    r <- lavaan.pl:::cpp_pl(x, f$engine$model, f$engine$data, o)
    if (o) r$gradient else r$pl
  }
  g <- pl(th, 1L)
  expect_equal(g, numDeriv::grad(function(x) pl(x, 0L), th), tolerance = 1e-6)
  hj <- lavaan.pl:::cpp_HJ(th, f$engine$model, f$engine$data)
  Hn <- -numDeriv::jacobian(function(x) pl(x, 1L), th) / n
  expect_equal(hj$H, (Hn + t(Hn)) / 2, tolerance = 1e-6)
  expect_equal(hj$gradient, g)
  sc <- unit_scores(f, th)
  expect_equal(colSums(sc), g)
  expect_equal(hj$J, crossprod(sc) / n)
  r <- lavaan.pl:::cpp_rho(th, f$engine$model, f$engine$data)
  expect_equal(r$jacobian, numDeriv::jacobian(function(x) lavaan.pl:::cpp_rho(x, f$engine$model, f$engine$data)$rho, th), tolerance = 1e-6)
}

# Weights: unit weights leave the fit unchanged, integer frequency weights reproduce the
# replicated data, and clusters give the centred cluster totals of the weighted scores
Yc <- Y[, 4:9]
M <- set_latpl_model(build_constr_mat(6, 2, "simple"))
set.seed(2)
w <- sample(1:3, nrow(Yc), replace = TRUE)
f <- fit_latpl(set_latpl_data(Yc), M)
expect_equal(fit_latpl(set_latpl_data(Yc, WEIGHTS = rep(2, nrow(Yc))), M)$theta, f$theta)
fw <- fit_latpl(set_latpl_data(Yc, WEIGHTS = w, WEIGHTS_TYPE = "frequency"), M)
fr <- fit_latpl(set_latpl_data(Yc[rep(seq_len(nrow(Yc)), w), ]), M)
expect_equal(fw$theta, fr$theta, tolerance = 1e-6)
th <- fw$theta + 0.05
hw <- lavaan.pl:::cpp_HJ(th, fw$engine$model, fw$engine$data)
hr <- lavaan.pl:::cpp_HJ(th, fr$engine$model, fr$engine$data)
expect_equal(hw$H, hr$H)
expect_equal(hw$J, hr$J)
cl <- rep(1:30, length.out = nrow(Yc))
fd <- fit_latpl(set_latpl_data(Yc, WEIGHTS = w, CLUSTER = cl), M)
f1 <- fit_latpl(set_latpl_data(Yc), M)
hd <- lavaan.pl:::cpp_HJ(th, fd$engine$model, fd$engine$data)
s1 <- t(sapply(seq_len(nrow(Yc)), function(i) {
  di <- f1$engine$data
  di$FREQ <- pairs_freq(Yc[i, , drop = FALSE], f1$data$cat)
  lavaan.pl:::cpp_pl(th, f1$engine$model, di, 1L)$gradient
}))
z <- rowsum(s1 * fd$data$weights, cl)
z <- sweep(z, 2, colMeans(z))
expect_equal(hd$J, 30 / 29 * crossprod(z) / nrow(Yc))
expect_equal(hd$gradient, colSums(s1 * fd$data$weights))

# Saturated engine: gradient and information blocks against numerical derivatives
Ds <- set_latpl_data(Ym[, 4:9], MISSING = "available.cases")
ds <- list(N = Ds$n, C_VEC = Ds$cat, FREQ = Ds$freq, UFREQ = Ds$ufreq, Y = Ds$Y, UW = Ds$uweights)
nt <- sum(Ds$cat - 1)
kap <- c(unlist(lapply(Ds$cat - 1, function(m) seq(-0.8, 0.8, length.out = m))), seq(.1, .3, length.out = 15))
sat <- function(x, o) lavaan.pl:::cpp_saturated(x[1:nt], x[-(1:nt)], ds, o, FALSE)
s <- sat(kap, 2L)
expect_equal(s$gradient, numDeriv::grad(function(x) sat(x, 0L)$ll, kap), tolerance = 1e-6)
Hs <- -numDeriv::jacobian(function(x) sat(x, 1L)$gradient, kap)
expect_equal(cbind(s$htt, s$htr), Hs[1:nt, ], tolerance = 1e-6)
expect_equal(s$hrr, diag(Hs)[-(1:nt)], tolerance = 1e-6)
expect_true(max(abs(Hs[-(1:nt), -(1:nt)] - diag(diag(Hs)[-(1:nt)]))) < 1e-6)
expect_equal(s$hschur, s$htt - s$htr %*% (t(s$htr) / s$hrr))

# Missing values: without NA all strategies give the complete-data fit; listwise drops;
# available-cases univariate counts weight each observed response by the missing items
f_l <- fit_latpl(set_latpl_data(Yc), M)
expect_equal(fit_latpl(set_latpl_data(Yc, MISSING = "pairwise"), M)$theta, f_l$theta)
expect_equal(fit_latpl(set_latpl_data(Yc, MISSING = "available.cases"), M)$theta, f_l$theta)
expect_message(d_l <- set_latpl_data(Ym), "incomplete observations removed")
expect_equal(d_l$n, sum(stats::complete.cases(Ym)))
m <- rowSums(is.na(Ym))
expect_equal(sum(D$ufreq), sum(m * (ncol(Ym) - m)))

# Robust C-estimator: gradient and H against numerical derivatives with capped cells; with a
# large tuning constant it reproduces the PML fit, with J centred; no PLRT nor available cases
Yr <- Y
set.seed(6)
Yr[1:30, ] <- sample(0:3, 270, replace = TRUE)
Yr[1:30, 1:3] <- 1 * (Yr[1:30, 1:3] >= 2)
Dr <- set_latpl_data(Yr)
Mr <- set_latpl_model(cm, REGMAT = B, CORRFLAG = FALSE, STDLV = FALSE)
f <- suppressWarnings(fit_latpl(Dr, Mr, CONTROL = set_ucminf_control(MAXEVAL = 1, ROBUST = 0.3)))
th <- c(unlist(lapply(Dr$cat - 1, function(m) seq(-0.8, 0.8, length.out = m))), seq(.4, .6, length.out = f$dims$nload),
        rep(-.5, f$dims$nvar), c(.3, .2, .25))
pl <- function(x, o) {
  r <- lavaan.pl:::cpp_pl(x, f$engine$model, f$engine$data, o)
  if (o) r$gradient else r$pl
}
expect_equal(pl(th, 1L), numDeriv::grad(function(x) pl(x, 0L), th), tolerance = 1e-6)
hj <- lavaan.pl:::cpp_HJ(th, f$engine$model, f$engine$data)
Hn <- -numDeriv::jacobian(function(x) pl(x, 1L), th) / Dr$n
expect_equal(hj$H, (Hn + t(Hn)) / 2, tolerance = 1e-6)
fp <- fit_latpl(Dr, Mr)
fr <- fit_latpl(Dr, Mr, CONTROL = set_ucminf_control(ROBUST = 1e8))
expect_equal(fr$theta, fp$theta, tolerance = 1e-6)
sc <- sapply(seq_len(Dr$n), function(i) {
  di <- fp$engine$data
  di$FREQ <- pairs_freq(Dr$Y[i, , drop = FALSE], Dr$cat)
  lavaan.pl:::cpp_pl(fp$theta, fp$engine$model, di, 1L)$gradient
})
sc <- sc - rowMeans(sc)
expect_equal(lavaan.pl:::cpp_HJ(fp$theta, fr$engine$model, fr$engine$data)$J, tcrossprod(sc) / Dr$n)
expect_error(fit_latpl(D, Mr, CONTROL = set_ucminf_control(ROBUST = 0.6)), "available.cases")
expect_error(compute_plrt(fit_latpl(Dr, Mr, CONTROL = set_ucminf_control(ROBUST = 0.6))), "robust")
