cm <- build_constr_mat(9, 3, "simple")
S <- matrix(.3, 3, 3)
diag(S) <- 1
Y <- sim_data(
  600,
  gen_loadings(cm, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  S,
  SEED = 3
)
D <- set_latpl_data(Y)
B <- matrix(0, 3, 3)
B[2, 1] <- B[3, 1] <- B[3, 2] <- NA
RC <- matrix(0, 9, 9)
RC[1, 2] <- RC[2, 1] <- NA
RC[5, 6] <- RC[6, 5] <- 0.1

# theta2par() and par2theta() invert each other
for (form in c("orthogonal", "lkj", "direct")) {
  for (stdlv in c(TRUE, FALSE)) {
    for (sem in c(FALSE, TRUE)) {
      M <- set_latpl_model(
        cm,
        STDLV = stdlv,
        CORRFLAG = form != "orthogonal",
        LATCOV = if (form == "direct") "direct" else "lkj",
        REGMAT = if (sem) B,
        RESCOV = if (sem) RC
      )
      nm <- theta_names(D, M)
      set.seed(1)
      th <- runif(length(nm), -0.3, 0.3)
      par <- theta2par(M, th)
      expect_equal(par2theta(M, par), th)
      expect_equal(par$residual_covariances[6, 5], if (sem) 0.1 else 0)
    }
  }
}
expect_equal(
  tail(
    theta_names(
      D,
      set_latpl_model(cm, LATCOV = "direct", REGMAT = B, RESCOV = RC)
    ),
    7
  ),
  c(
    "psi[1,2]",
    "psi[1,3]",
    "psi[2,3]",
    "beta[2,1]",
    "beta[3,1]",
    "beta[3,2]",
    "theta[1,2]"
  )
)
expect_equal(
  tail(theta_names(D, set_latpl_model(cm, STDLV = FALSE)), 4),
  c("eta[2,3]", "logsd[1]", "logsd[2]", "logsd[3]")
)

# Linear constraints: equal loadings and a fixed threshold, through the developer interface
M0 <- set_latpl_model(cm)
nm <- theta_names(D, M0)
d <- length(nm)
A <- rbind(
  replace(numeric(d), match(c("lambda[2,1]", "lambda[3,1]"), nm), c(1, -1)),
  replace(numeric(d), match("tau[4,2]", nm), 1)
)
f <- fit_latpl(
  D,
  set_latpl_model(cm, CONSTRAINTS = list(A = A, b = c(0, 0.25)))
)
expect_equal(f$dims$npar, d - 2)
expect_equal(as.vector(A %*% f$theta), c(0, 0.25))
K <- f$lincon$K
ngr <- function(theta) {
  -lavaan.pl:::cpp_pl(theta, f$engine$model, f$engine$data, 1L)$gradient / D$n
}
expect_true(max(abs(crossprod(K, ngr(f$theta)))) < 1e-5)
v <- compute_var(f)
Hnum <- numDeriv::jacobian(
  function(x) as.vector(crossprod(K, ngr(as.vector(K %*% x + f$lincon$c)))),
  f$theta[f$lincon$free_idx]
)
expect_equal(v$H, (Hnum + t(Hnum)) / 2, tolerance = 1e-4)
expect_error(
  fit_latpl(
    D,
    set_latpl_model(
      cm,
      CONSTRAINTS = list(
        A = matrix(replace(numeric(d), match("eta[1,2]", nm), 1), 1),
        b = 0
      )
    )
  ),
  "LATCOV"
)
expect_error(
  fit_latpl(
    D,
    set_latpl_model(
      cm,
      CONSTRAINTS = list(
        A = matrix(replace(numeric(d), c(1, 28), c(1, -1)), 1),
        b = 0
      )
    )
  ),
  "either thresholds"
)
expect_error(
  fit_latpl(
    D,
    set_latpl_model(
      cm,
      CONSTRAINTS = list(A = rbind(A[1, ], A[1, ]), b = c(0, 1))
    )
  ),
  "inconsistent"
)

# LKJ and direct parametrisations
for (stdlv in c(TRUE, FALSE)) {
  fl <- fit_latpl(D, set_latpl_model(cm, STDLV = stdlv))
  fd <- fit_latpl(D, set_latpl_model(cm, STDLV = stdlv, LATCOV = "direct"))
  expect_equal(fd$numFit$value, fl$numFit$value, tolerance = 1e-8)
  expect_equal(
    get_par(fd, OPTION = "transformed"),
    get_par(fl, OPTION = "transformed"),
    tolerance = 1e-4
  )
  vl <- compute_var(fl)
  vd <- compute_var(fd)
  expect_equal(vd$asymptotic_variance, vl$asymptotic_variance, tolerance = 1e-3)
  expect_equal(
    compute_plrt(fd, vd)$model[c("stat", "df", "stat.scaled", "df.scaled")],
    compute_plrt(fl, vl)$model[c("stat", "df", "stat.scaled", "df.scaled")],
    tolerance = 1e-3
  )
}
Md <- set_latpl_model(cm, LATCOV = "direct")
nd <- theta_names(D, Md)
f0 <- fit_latpl(
  D,
  set_latpl_model(
    cm,
    LATCOV = "direct",
    CONSTRAINTS = list(
      A = matrix(replace(numeric(d), match("psi[1,3]", nd), 1), 1),
      b = 0
    )
  )
)
expect_equal(get_par(f0, OPTION = "list")$latent_correlations[1, 3], 0)
expect_error(
  set_latpl_model(cm, CORRFLAG = FALSE, LATCOV = "direct"),
  "CORRFLAG"
)

# saturated structural model equivalent to the correlated-factor model
fc <- fit_latpl(D, set_latpl_model(cm, STDLV = FALSE))
fs <- fit_latpl(
  D,
  set_latpl_model(cm, STDLV = FALSE, CORRFLAG = FALSE, REGMAT = B)
)
expect_equal(fs$numFit$value, fc$numFit$value, tolerance = 1e-6)
gs <- get_par(fs, OPTION = "list")
Ts <- solve(diag(3) - gs$regressions)
expect_equal(
  Ts %*% gs$latent_correlations %*% t(Ts),
  get_par(fc, OPTION = "list")$latent_correlations,
  tolerance = 1e-3
)
expect_equal(
  compute_plrt(fs)$model$stat,
  compute_plrt(fc)$model$stat,
  tolerance = 1e-4
)
fn <- fit_latpl(
  D,
  set_latpl_model(
    cm,
    STDLV = FALSE,
    CORRFLAG = FALSE,
    REGMAT = replace(B, 3, 0)
  )
)
expect_equal(
  compute_plrt(fn, H1 = fc)$model[c("stat", "df")],
  compute_plrt(fn, H1 = fs)$model[c("stat", "df")],
  tolerance = 1e-4
)

# Fixed residual covariances and othrs
RF <- matrix(0, 9, 9)
RF[1, 2] <- RF[2, 1] <- 0.1
RF[4, 7] <- RF[7, 4] <- NA
fr <- fit_latpl(D, set_latpl_model(cm, REGMAT = B, RESCOV = RF))
expect_equal(get_par(fr, OPTION = "list")$residual_covariances[2, 1], 0.1)
expect_error(set_latpl_model(cm, REGMAT = diag(3)), "diagonal")
expect_error(set_latpl_model(cm, RESCOV = replace(RC, 2, 0.2)), "symmetric")

# Two unconstrained groups reproduce the separate fits; equality constraints across groups
Y1 <- Y[1:150, 4:9]
Y2 <- Y[151:300, 4:9]
cm2 <- build_constr_mat(6, 2, "simple")
M2 <- set_latpl_model(cm2)
f1 <- fit_latpl(set_latpl_data(Y1), M2)
f2 <- fit_latpl(set_latpl_data(Y2), M2)
Dg <- set_latpl_data(rbind(Y1, Y2), GROUP = rep(c("b", "a"), each = 150))
Mg <- set_latpl_groups(list(M2, M2))
fg <- fit_latpl(Dg, Mg)
expect_identical(Dg$group, rep(c(0, 1), each = 150))
expect_equal(fg$theta, c(f1$theta, f2$theta), tolerance = 1e-5)
expect_equal(
  lavaan.pl:::cpp_pl(fg$theta, fg$engine$model, fg$engine$data)$pl,
  lavaan.pl:::cpp_pl(f1$theta, f1$engine$model, f1$engine$data)$pl + lavaan.pl:::cpp_pl(f2$theta, f2$engine$model, f2$engine$data)$pl
)
expect_equal(
  compute_var(fg)$asymptotic_variance / 300,
  c(compute_var(f1)$asymptotic_variance, compute_var(f2)$asymptotic_variance) / 150,
  tolerance = 1e-4
)
expect_length(get_par(fg), 2)
expect_equal(get_par(fg)[[2]], theta2par(M2, f2$theta), tolerance = 1e-4)
expect_equal(par2theta(Mg, get_par(fg)), fg$theta)
nm <- theta_names(Dg, Mg)
expect_identical(nm, c(theta_names(Dg, M2), paste0(theta_names(Dg, M2), ".g2")))
A <- matrix(0, 1, length(nm))
A[1, match(c("lambda[1,1]", "lambda[1,1].g2"), nm)] <- c(1, -1)
fc <- fit_latpl(Dg, set_latpl_groups(list(M2, M2), list(A = A, b = 0)))
expect_equal(fc$theta[match("lambda[1,1]", nm)], fc$theta[match("lambda[1,1].g2", nm)])
expect_equal(fc$dims$npar, fg$dims$npar - 1)
expect_error(fit_latpl(Dg, M2), "groups")
expect_error(fit_latpl(set_latpl_data(Y1), Mg), "groups")
Yg <- rbind(Y1, Y2)
expect_error(set_latpl_data(Yg, GROUP = rep(1:2, 150), VALDATA = Yg), "VALDATA")
expect_error(set_latpl_data(Yg, GROUP = rep(1:2, 150), CLUSTER = rep(1:30, 10)), "CLUSTER")

# Frequency weights equal replicated rows, within groups
set.seed(2)
w <- sample(1:3, 300, replace = TRUE)
gr <- rep(1:2, each = 150)
fw <- fit_latpl(set_latpl_data(Yg, WEIGHTS = w, WEIGHTS_TYPE = "frequency", GROUP = gr), Mg)
fr <- fit_latpl(set_latpl_data(Yg[rep(1:300, w), ], GROUP = gr[rep(1:300, w)]), Mg)
expect_equal(fw$theta, fr$theta, tolerance = 1e-6)
