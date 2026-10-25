cm <- build_constrMat(6, 2, "simple")
S <- matrix(.3, 2, 2); diag(S) <- 1
Y <- sim_data(300, gen_loadings(cm, SEED = 1, LB = .4, UB = .8), c(-1, 0.3, 1), S, SEED = 3)
f <- fit_plFA(set_plfa_data(Y), set_plfa_model(cm))
v <- computeVar(f)
r <- computePLRT(f, v)
n <- f$dims$n; nt <- f$dims$nthr; psi <- (nt + 1):f$dims$npar
sat <- function(tau, rho, order, scores) lavaan.pl:::cpp_saturated(tau, rho, f$data$freq, f$data$Y, f$dims$cat, order, scores)

# Saturated fit is stationary
expect_true(max(abs(sat(r$saturated$tau, r$saturated$rho, 1, FALSE)$gradient)) / n < 1e-6)
expect_equal(r$model$df, nt + f$dims$pairs - f$dims$npar)
expect_equal(r$baseline$df, f$dims$pairs)

# Mean and variance adjustment matches the dense computation
rho_theta <- function(theta) {
  gp <- getPar(replace(f, "theta", list(theta)), OPTION = "list")
  S <- gp$loadings %*% gp$latent_correlations %*% t(gp$loadings)
  S[upper.tri(S)]
}
s <- sat(f$theta[1:nt], rho_theta(f$theta), 2, TRUE)
Hs <- rbind(cbind(s$htt, s$htr), cbind(t(s$htr), diag(s$hrr))) / n
Js <- tcrossprod(rbind(s$st, s$sr)) / n
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

# Nested test: orthogonal against correlated factors, exact restriction on the latent correlation
f0 <- fit_plFA(set_plfa_data(Y), set_plfa_model(cm, CORRFLAG = FALSE))
rn <- computePLRT(f0, H1 = f)$model
hj <- lavaan.pl:::cpp_sample_estimators_HJ(
  THETA = c(f0$theta, 0), FREQ = f$data$freq, DATA = f$data$Y, C_VEC = f$dims$cat, CONSTRMAT = f$model$CONSTRMAT,
  CONSTRLOGSD = f$model$CONSTRLOGSD, LLC = f$model$LLC, N = n, CORRFLAG = TRUE, NTHR = nt, NLOAD = f$dims$nload,
  NCORR = 1, NVAR = 0
)
iH1 <- solve(hj$H)
k <- f$dims$npar
expect_equal(rn$df, 1)
expect_equal(rn$scaling.factor, (iH1 %*% hj$J %*% iH1)[k, k] / iH1[k, k])
expect_equal(rn$df.scaled, 1)
expect_equal(rn$stat.scaled, rn$stat / rn$scaling.factor)
expect_error(computePLRT(f, H1 = f0), "not nested")

# anova() leaves other model classes to their own methods
l1 <- lm(mpg ~ wt, mtcars)
l2 <- lm(mpg ~ wt + hp, mtcars)
expect_equal(lavaan.pl::anova(l1, l2), stats::anova(l1, l2))
