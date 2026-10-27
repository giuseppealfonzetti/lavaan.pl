cm <- build_constr_mat(6, 2, "simple")
S <- matrix(.3, 2, 2); diag(S) <- 1
Y <- sim_data(300, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
Ym <- Y
set.seed(4)
for (j in c(2:3, 5:6)) Ym[runif(300) < c(.1, .3, .3, .6)[Y[, if (j < 4) 1 else 4] + 1], j] <- NA

# Available cases with missing values, design weights with clusters, then complete data
set.seed(5)
Dwc <- set_latpl_data(Y, WEIGHTS = runif(300, .2, 3), CLUSTER = rep(1:40, length.out = 300))
for (D in list(set_latpl_data(Ym, MISSING = "available.cases"), Dwc, set_latpl_data(Y))) {
  f <- fit_latpl(D, set_latpl_model(cm))
  v <- compute_var(f)
  r <- compute_plrt(f, v)
  n <- f$dims$n; nt <- f$dims$nthr; psi <- (nt + 1):f$dims$npar
  dat <- list(C_VEC = f$data$cat, FREQ = f$data$freq, UFREQ = f$data$ufreq, Y = f$data$Y, UW = f$data$uweights)
  sat <- function(tau, rho, order, scores) lavaan.pl:::cpp_saturated(tau, rho, dat, order, scores)

  # Saturated fit is stationary
  expect_true(max(abs(sat(r$saturated$tau, r$saturated$rho, 1, FALSE)$gradient)) / n < 1e-6)
  expect_equal(r$model$df, nt + f$dims$pairs - f$dims$npar)
  expect_equal(r$baseline$df, f$dims$pairs)

  # Mean and variance adjustment matches the dense computation
  rho_theta <- function(theta) {
    gp <- get_par(replace(f, "theta", list(theta)), OPTION = "list")
    S <- gp$loadings %*% gp$latent_correlations %*% t(gp$loadings)
    S[upper.tri(S)]
  }
  s <- sat(f$theta[1:nt], rho_theta(f$theta), 2, TRUE)
  Hs <- rbind(cbind(s$htt, s$htr), cbind(t(s$htr), diag(s$hrr))) / n
  Sc <- rbind(s$st, s$sr)
  if (length(f$data$cluster)) {
    Z <- t(rowsum(t(Sc) * f$data$weights, f$data$cluster))
    Z <- Z - rowMeans(Z)
    Js <- ncol(Z) / (ncol(Z) - 1) * tcrossprod(Z) / n
  } else {
    Js <- tcrossprod(Sc) / n
  }
  iHs <- solve(Hs)
  A <- solve(iHs[-(1:nt), -(1:nt)])
  P1 <- A %*% (iHs %*% Js %*% iHs)[-(1:nt), -(1:nt)]
  iH <- solve(v$H)
  G <- (iH %*% v$J %*% iH)[psi, psi]
  P0 <- solve(iH[psi, psi]) %*% G
  delta <- numDeriv::jacobian(rho_theta, f$theta)[, psi]
  E <- sum(diag(P1)) - sum(diag(P0))
  V <- 2 * sum(diag(P1 %*% P1)) + 2 * sum(diag(P0 %*% P0)) - 4 * sum(diag(t(delta) %*% A %*% delta %*% G %*% P0))
  expect_equal(r$model$scaling.factor, (V / 2) / E)
  expect_equal(r$model$df.scaled, E^2 / (V / 2))
  expect_equal(r$model$stat, 2 * (r$saturated$logl - s$ll))

  # Closed-form independence thresholds maximise the independence pairwise likelihood
  expect_true(max(abs(sat(r$independence$tau, rep(0, f$dims$pairs), 1, FALSE)$gradient[1:nt])) / n < 1e-6)
  expect_equal(r$baseline$stat, 2 * (r$saturated$logl - r$independence$logl))
}

# Nested test: orthogonal against correlated factors, exact restriction on the latent correlation
f0 <- fit_latpl(set_latpl_data(Y), set_latpl_model(cm, CORRFLAG = FALSE))
rn <- compute_plrt(f0, H1 = f)$model
hj <- lavaan.pl:::cpp_HJ(c(f0$theta, 0), f$engine$model, f$engine$data)
iH1 <- solve(hj$H)
k <- f$dims$npar
expect_equal(rn$df, 1)
expect_equal(rn$scaling.factor, (iH1 %*% hj$J %*% iH1)[k, k] / iH1[k, k])
expect_equal(rn$df.scaled, 1)
expect_equal(rn$stat.scaled, rn$stat / rn$scaling.factor)
expect_error(compute_plrt(f, H1 = f0), "not nested")

# anova() leaves other model classes to their own methods
l1 <- lm(mpg ~ wt, mtcars)
l2 <- lm(mpg ~ wt + hp, mtcars)
expect_equal(lavaan.pl::anova(l1, l2), stats::anova(l1, l2))

# Constrained thresholds: full-vector adjustment (lavaan's plrt2) matches the dense computation
D <- set_latpl_data(Y)
M0 <- set_latpl_model(cm)
nm <- theta_names(D, M0)
A <- matrix(replace(numeric(length(nm)), match("tau[2,1]", nm), 1), 1)
f <- fit_latpl(D, set_latpl_model(cm, CONSTRAINTS = list(A = A, b = -1)))
v <- compute_var(f)
r <- compute_plrt(f, v)
n <- f$dims$n; nt <- f$dims$nthr; lc <- f$lincon
rho_free <- function(x) {
  gp <- get_par(replace(f, "theta", list(as.vector(lc$K %*% x + lc$c))), OPTION = "list")
  S <- gp$loadings %*% gp$latent_correlations %*% t(gp$loadings)
  S[upper.tri(S)]
}
xf <- f$theta[lc$free_idx]
s <- lavaan.pl:::cpp_saturated(f$theta[1:nt], rho_free(xf), list(C_VEC = f$data$cat, FREQ = f$data$freq, UFREQ = numeric(0),
                                                                 Y = f$data$Y, UW = numeric(0)), 2, TRUE)
Hs <- rbind(cbind(s$htt, s$htr), cbind(t(s$htr), diag(s$hrr))) / n
P1 <- solve(Hs, tcrossprod(rbind(s$st, s$sr)) / n)
P0 <- solve(v$H, v$J)
delta <- rbind(lc$K[1:nt, ], numDeriv::jacobian(rho_free, xf))
E <- sum(diag(P1)) - sum(diag(P0))
V <- 2 * sum(diag(P1 %*% P1)) + 2 * sum(diag(P0 %*% P0)) - 4 * sum(diag(t(delta) %*% Hs %*% delta %*% P0 %*% P0 %*% solve(v$H)))
expect_equal(r$model$df, nt + f$dims$pairs - f$dims$npar)
expect_equal(r$model$scaling.factor, (V / 2) / E)
expect_equal(r$model$df.scaled, E^2 / (V / 2))
