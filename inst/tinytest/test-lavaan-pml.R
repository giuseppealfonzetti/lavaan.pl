if (!at_home()) {
  exit_file("Not run on CRAN")
}

grid <- list(
  "binary, one factor" = list(
    model = "eta =~ y1 + y2 + y3 + y4 + y5",
    data = LSAT,
    std.lv = TRUE
  ),
  "ordinal, std.lv" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3",
    data = bfi,
    std.lv = TRUE
  ),
  "ordinal, marker" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3",
    data = bfi,
    std.lv = FALSE
  ),
  "covariance written first" = list(
    model = "A ~~ C\nA =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3",
    data = bfi,
    std.lv = FALSE
  ),
  "orthogonal factors" = list(
    model = "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ 0*C",
    data = bfi,
    std.lv = TRUE
  ),
  "fixed loading" = list(
    model = "A =~ A1 + A2 + A3\nC =~ 0.6*C1 + C2 + C3",
    data = bfi,
    std.lv = TRUE
  ),
  "cross-loading" = list(
    model = "A =~ A1 + A2 + A3 + C1\nC =~ C1 + C2 + C3",
    data = bfi,
    std.lv = TRUE
  ),
  "equality constraints" = list(
    model = "A =~ A1 + a*A2 + b*A3 + c*A4 + A5\nC =~ C1 + d*C2 + d*C3\na == 0.5*b + 0.25*c",
    data = bfi,
    std.lv = TRUE
  ),
  "defined parameter" = list(
    model = "A =~ A1 + a*A2 + b*A3\nC =~ C1 + C2 + C3\nab := a*b",
    data = bfi,
    std.lv = FALSE
  ),
  "theta, std.lv" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3",
    data = bfi,
    std.lv = TRUE,
    parameterization = "theta"
  ),
  "theta, marker, residual covariance" = list(
    model = "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3\nA1 ~~ A2",
    data = bfi,
    std.lv = FALSE,
    parameterization = "theta"
  ),
  "latent mean" = list(
    model = "A =~ A1 + A2 + A3 + A4\nA ~ 1\nA1 | 0*t1 + t2 + t3 + t4 + t5",
    data = bfi,
    std.lv = TRUE
  )
)

for (nm in names(grid)) {
  g <- grid[[nm]]
  pz <- if (is.null(g$parameterization)) "delta" else g$parameterization
  ours <- cfa(g$model, g$data, std.lv = g$std.lv, parameterization = pz)
  ref <- lavaan::cfa(g$model, g$data, std.lv = g$std.lv, estimator = "PML", parameterization = pz)

  local({
    expect_identical(as.character(class(ours)), "lavaan")
    expect_true(lavaan::lavInspect(ours, "converged"))
    pe <- parameterEstimates(ours)
    pr <- parameterEstimates(ref)
    expect_identical(pe[, 1:3], pr[, 1:3])
    expect_equal(pe$est, pr$est, tolerance = 1e-4)
    expect_equal(pe$se, pr$se, tolerance = 1e-3)
    expect_equal(fitted(ours), fitted(ref), tolerance = 1e-4)
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
  })
}

# mean.var.adjusted test and fit measures match lavaan PML
local({
  fm <- c(
    "chisq",
    "df",
    "chisq.scaled",
    "df.scaled",
    "pvalue.scaled",
    "chisq.scaling.factor",
    "baseline.chisq",
    "baseline.df",
    "baseline.chisq.scaled",
    "cfi.scaled",
    "tli.scaled",
    "rmsea.scaled"
  )
  for (nm in c(
    "binary, one factor",
    "ordinal, std.lv",
    "ordinal, marker",
    "cross-loading",
    "equality constraints",
    "theta, std.lv",
    "theta, marker, residual covariance",
    "latent mean"
  )) {
    g <- grid[[nm]]
    pz <- if (is.null(g$parameterization)) "delta" else g$parameterization
    ours <- cfa(g$model, g$data, std.lv = g$std.lv, parameterization = pz, test = "mean.var.adjusted")
    ref <- lavaan::cfa(
      g$model,
      g$data,
      std.lv = g$std.lv,
      parameterization = pz,
      estimator = "PML",
      test = "mean.var.adjusted"
    )
    expect_equal(
      ours@test[[1]][c("stat", "df")],
      ref@test[[1]][c("stat", "df")],
      tolerance = 1e-3
    )
    expect_equal(
      ours@test[[2]][c("stat", "df", "pvalue", "scaling.factor")],
      ref@test[[2]][c("stat", "df", "pvalue", "scaling.factor")],
      tolerance = 1e-3
    )
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
  }
})

# Nested test matches lavaan's PLRT
local({
  m1 <- "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3 + C4 + C5"
  f1 <- cfa(m1, bfi, std.lv = TRUE, test = "mean.var.adjusted")
  h0 <- list(
    cfa(m1, bfi, std.lv = TRUE, orthogonal = TRUE, test = "mean.var.adjusted"),
    cfa(
      "A =~ A1 + a*A2 + a*A3 + A4 + A5\nC =~ C1 + b*C2 + b*C3 + C4 + C5",
      bfi,
      std.lv = TRUE,
      test = "mean.var.adjusted"
    )
  )
  for (f0 in h0) {
    ours <- lavTestLRT(f1, f0)
    expect_equal(lavaan.pl::anova(f1, f0), ours)
    ref <- suppressWarnings(lavaan::lavTestLRT(f1, f0))
    m10 <- lavaan:::lav_test_diff_m10(m1 = f1, m0 = f0, test = FALSE)
    iH <- solve(lavaan::lavTech(m10, "information.observed"))
    J <- lavaan::lavTech(m10, "first.order")
    A <- lavaan:::lav_test_diff_a(m1 = f1, m0 = f0, method = "exact")
    M <- (A %*% iH %*% J %*% iH %*% t(A)) %*% solve(A %*% iH %*% t(A))
    plrt <- 2 * (f1@optim$logl - f0@optim$logl)
    expect_equal(
      ours[2, "Chisq diff"],
      sum(diag(M)) / sum(diag(M %*% M)) * plrt,
      tolerance = 1e-4
    )
    expect_equal(
      ours[2, "Df diff"],
      sum(diag(M))^2 / sum(diag(M %*% M)),
      tolerance = 1e-4
    )
    cols <- c("Df", "PL_AIC", "PL_BIC", "Chisq")
    expect_equal(unlist(ours[, cols]), unlist(ref[, cols]), tolerance = 1e-6)
  }
})

# Missing data match lavaan PML
local({
  cm <- build_constr_mat(12, 2, "simple")
  L <- cm
  L[is.na(L)] <- c(.4, .5, .6, .7, .8, .9, .4, .5, .6, .7, .8, .9)
  Yc <- sim_data(
    1000,
    L,
    c(-1.5, .5, 1.5),
    matrix(c(1, .6, .6, 1), 2),
    SEED = 4
  )
  Y <- Yl <- Yc
  set.seed(5)
  pm <- c(.119, .269, .269, .731)
  for (j in 2:6) {
    Y[runif(1000) < pm[Yc[, 1] + 1], j] <- NA
  }
  for (j in 8:12) {
    Y[runif(1000) < pm[Yc[, 7] + 1], j] <- NA
  }
  Yl[, c(2, 8)] <- Y[, c(2, 8)]
  colnames(Y) <- colnames(Yl) <- paste0("y", 1:12)
  mod <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6\nf2 =~ y7 + y8 + y9 + y10 + y11 + y12"
  fm <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "pvalue.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled",
    "tli.scaled"
  )
  for (mi in c("listwise", "pairwise", "available.cases")) {
    dat <- as.data.frame(lapply(
      as.data.frame(if (mi == "listwise") Yl else Y),
      ordered
    ))
    ours <- suppressWarnings(cfa(
      mod,
      dat,
      std.lv = TRUE,
      missing = mi,
      test = "mean.var.adjusted"
    ))
    ref <- lavaan::cfa(
      mod,
      dat,
      std.lv = TRUE,
      estimator = "PML",
      missing = mi,
      test = "mean.var.adjusted"
    )
    pe <- parameterEstimates(ours)
    pr <- parameterEstimates(ref)
    expect_equal(
      lavaan::lavInspect(ours, "nobs"),
      lavaan::lavInspect(ref, "nobs")
    )
    expect_equal(pe$est, pr$est, tolerance = 1e-4)
    expect_equal(pe$se, pr$se, tolerance = 1e-3)
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
    if (mi == "listwise") {
      next
    }
    f0 <- suppressWarnings(cfa(
      mod,
      dat,
      std.lv = TRUE,
      orthogonal = TRUE,
      missing = mi,
      test = "mean.var.adjusted"
    ))
    m10 <- lavaan:::lav_test_diff_m10(m1 = ours, m0 = f0, test = FALSE)
    iH <- solve(lavaan::lavTech(m10, "information.observed"))
    J <- lavaan::lavTech(m10, "first.order")
    A <- lavaan:::lav_test_diff_a(m1 = ours, m0 = f0, method = "exact")
    M <- (A %*% iH %*% J %*% iH %*% t(A)) %*% solve(A %*% iH %*% t(A))
    expect_equal(
      lavTestLRT(ours, f0)[2, "Chisq diff"],
      sum(diag(M)) / sum(diag(M %*% M)) * 2 * (ours@optim$logl - f0@optim$logl),
      tolerance = 1e-4
    )
  }
  expect_error(cfa(
    mod,
    as.data.frame(lapply(as.data.frame(Y), ordered)),
    std.lv = TRUE,
    missing = "doubly.robust"
  ))
})

# Sampling weights and clusters match lavaan PML
local({
  d <- as.data.frame(bfi)
  items <- c("A1", "A2", "A3", "A4", "C1", "C2", "C3", "C4")
  for (v in items) {
    d[[v]] <- ordered(d[[v]])
  }
  set.seed(1)
  d <- d[complete.cases(d[, items]), ][sample(2000, 400), ]
  d$w <- runif(400, .2, 3)
  d$cl <- rep(1:40, length.out = 400)
  dm <- d
  for (j in c("A2", "A3", "C2", "C3")) {
    dm[[j]][runif(400) < .15] <- NA
  }
  mod <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4"
  fm <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled"
  )
  cases <- list(
    list(data = d, a = list(sampling.weights = "w")),
    list(
      data = d,
      a = list(sampling.weights = "w", sampling.weights.type = "frequency")
    ),
    list(data = d, a = list(cluster = "cl")),
    list(data = d, a = list(sampling.weights = "w", cluster = "cl")),
    list(data = dm, a = list(sampling.weights = "w", missing = "pairwise"))
  )
  for (cs in cases) {
    tst <- if (is.null(cs$a$cluster)) "mean.var.adjusted" else "none"
    ours <- suppressWarnings(do.call(
      cfa,
      c(list(mod, cs$data, std.lv = TRUE, test = tst), cs$a)
    ))
    ref <- suppressWarnings(do.call(
      lavaan::cfa,
      c(list(mod, cs$data, std.lv = TRUE, estimator = "PML", test = tst), cs$a)
    ))
    expect_equal(
      parameterEstimates(ours)$est,
      parameterEstimates(ref)$est,
      tolerance = 1e-4
    )
    expect_equal(
      parameterEstimates(ours)$se,
      parameterEstimates(ref)$se,
      tolerance = 1e-3
    )
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
    if (tst != "none") {
      expect_equal(
        fitMeasures(ours, fm),
        fitMeasures(ref, fm),
        tolerance = 1e-3
      )
    }
  }
  expect_error(
    cfa("A =~ A1 + a*A2 + A3\nC =~ C1 + C2 + C3\na > 2", d, std.lv = TRUE),
    "Inequality constraints"
  )
  expect_error(
    cfa("A =~ A1 + a*A2 + b*A3\nC =~ C1 + C2 + C3\na > b", d, std.lv = TRUE),
    "Inequality constraints"
  )
  expect_error(cfa(mod, d, std.lv = TRUE, bounds = "standard"), "bounds")
  expect_error(
    cfa(
      mod,
      d,
      std.lv = TRUE,
      sampling.weights = "w",
      sampling.weights.normalization = "none"
    ),
    "normalization"
  )
})

# Linear constraints match lavaan PML
local({
  d <- as.data.frame(bfi)
  items <- c("A1", "A2", "A3", "A4", "C1", "C2", "C3", "C4")
  for (v in items) {
    d[[v]] <- ordered(d[[v]])
  }
  set.seed(1)
  d <- d[complete.cases(d[, items]), ][sample(2000, 600), ]
  base <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4"
  fm <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled"
  )
  cases <- list(
    list(
      m = "A =~ A1 + a*A2 + b*A3 + A4\nC =~ C1 + C2 + C3 + C4\na == 0.5*b + 0.1",
      lv = TRUE
    ),
    list(m = paste(base, "\nA1 | -0.4*t1"), lv = TRUE),
    list(m = paste(base, "\nA ~~ v*A\nC ~~ v*C"), lv = FALSE),
    list(m = paste(base, "\nA ~~ v1*A\nC ~~ v2*C\nv1 == 2*v2"), lv = FALSE)
  )
  for (cs in cases) {
    ours <- suppressWarnings(cfa(
      cs$m,
      d,
      std.lv = cs$lv,
      test = "mean.var.adjusted"
    ))
    ref <- suppressWarnings(lavaan::cfa(
      cs$m,
      d,
      std.lv = cs$lv,
      estimator = "PML",
      test = "mean.var.adjusted"
    ))
    expect_equal(
      parameterEstimates(ours)$est,
      parameterEstimates(ref)$est,
      tolerance = 1e-4
    )
    expect_equal(
      parameterEstimates(ours)$se,
      parameterEstimates(ref)$se,
      tolerance = 1e-3
    )
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
  }
  # Equal thresholds
  f <- cfa(paste(base, "\nA3 | t*t3\nA4 | t*t3"), d, std.lv = TRUE)
  g <- lavaan::lavTech(f, "gradient")
  jc <- f@Model@ceq.JAC
  N <- qr.Q(qr(t(jc)), complete = TRUE)[, -seq_len(nrow(jc)), drop = FALSE]
  expect_true(max(abs(crossprod(N, g))) < 1e-4 * max(abs(g)))
})

# Fixed and constrained latent covariances match lavaan PML
local({
  d <- as.data.frame(bfi)
  items <- c(paste0("A", 1:4), paste0("C", 1:4), paste0("E", 1:4))
  for (v in items) {
    d[[v]] <- ordered(d[[v]])
  }
  set.seed(1)
  d <- d[complete.cases(d[, items]), ][sample(2000, 800), ]
  three <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4\nE =~ E1 + E2 + E3 + E4"
  fm <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled"
  )
  cases <- list(
    list(m = paste(three, "\nA ~~ 0*E"), lv = TRUE),
    list(m = paste(three, "\nA ~~ 0.3*C"), lv = TRUE),
    list(m = paste(three, "\nA ~~ r*C\nA ~~ r*E"), lv = TRUE),
    list(m = paste(three, "\nA ~~ 0.1*E"), lv = FALSE),
    list(m = paste(three, "\nA ~~ r*C\nA ~~ r*E"), lv = FALSE)
  )
  for (cs in cases) {
    ours <- suppressWarnings(cfa(
      cs$m,
      d,
      std.lv = cs$lv,
      test = "mean.var.adjusted"
    ))
    ref <- suppressWarnings(lavaan::cfa(
      cs$m,
      d,
      std.lv = cs$lv,
      estimator = "PML",
      test = "mean.var.adjusted"
    ))
    expect_identical(ours@external$lavaan.pl$model$LATCOV, "direct")
    expect_equal(
      parameterEstimates(ours)$est,
      parameterEstimates(ref)$est,
      tolerance = 1e-4
    )
    expect_equal(
      parameterEstimates(ours)$se,
      parameterEstimates(ref)$se,
      tolerance = 1e-3
    )
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
  }
  # Bifactor
  mb <- paste(
    three,
    "\nG =~ A1 + A2 + C1 + C2 + E1 + E2\nG ~~ 0*A\nG ~~ 0*C\nG ~~ 0*E"
  )
  ours <- suppressWarnings(cfa(mb, d, std.lv = TRUE))
  ref <- suppressWarnings(lavaan::cfa(mb, d, std.lv = TRUE, estimator = "PML"))
  expect_true(ours@optim$logl >= ref@optim$logl - 1e-6)
  # Nested test across parametrisations
  f1 <- cfa(three, d, std.lv = TRUE, test = "mean.var.adjusted")
  f0 <- cfa(
    paste(three, "\nA ~~ 0*E"),
    d,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  m10 <- lavaan:::lav_test_diff_m10(m1 = f1, m0 = f0, test = FALSE)
  iH <- solve(lavaan::lavTech(m10, "information.observed"))
  J <- lavaan::lavTech(m10, "first.order")
  A <- lavaan:::lav_test_diff_a(m1 = f1, m0 = f0, method = "exact")
  M <- (A %*% iH %*% J %*% iH %*% t(A)) %*% solve(A %*% iH %*% t(A))
  expect_equal(
    lavTestLRT(f1, f0)[2, "Chisq diff"],
    sum(diag(M)) / sum(diag(M %*% M)) * 2 * (f1@optim$logl - f0@optim$logl),
    tolerance = 1e-4
  )
})

# Residual covs, structural regres and higher-order factors match lavaan PML
local({
  d <- as.data.frame(bfi)
  items <- c(paste0("A", 1:4), paste0("C", 1:4), paste0("E", 1:4))
  for (v in items) {
    d[[v]] <- ordered(d[[v]])
  }
  set.seed(1)
  d <- d[complete.cases(d[, items]), ][sample(2000, 800), ]
  three <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4\nE =~ E1 + E2 + E3 + E4"
  fm <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled"
  )
  cases <- list(
    list(m = paste(three, "\nA1 ~~ A2\nC1 ~~ E1"), lv = TRUE),
    list(m = paste(three, "\nA1 ~~ A2\nC1 ~~ E1"), lv = FALSE),
    list(m = paste(three, "\nC ~ A\nE ~ A + C"), lv = TRUE),
    list(m = paste(three, "\nC ~ A\nE ~ A + C"), lv = FALSE),
    list(m = paste(three, "\nE ~ A + C"), lv = FALSE),
    list(m = paste(three, "\nE ~ b*A + b*C"), lv = TRUE),
    list(m = paste(three, "\nC ~ A\nE ~ A + C\nA1 ~~ A2"), lv = TRUE),
    list(m = paste(three, "\nG =~ A + C + E"), lv = FALSE)
  )
  for (cs in cases) {
    ours <- suppressWarnings(sem(
      cs$m,
      d,
      std.lv = cs$lv,
      test = "mean.var.adjusted"
    ))
    ref <- suppressWarnings(lavaan::sem(
      cs$m,
      d,
      std.lv = cs$lv,
      estimator = "PML",
      test = "mean.var.adjusted"
    ))
    expect_equal(
      parameterEstimates(ours)$est,
      parameterEstimates(ref)$est,
      tolerance = 1e-4
    )
    expect_equal(
      parameterEstimates(ours)$se,
      parameterEstimates(ref)$se,
      tolerance = 1e-3
    )
    expect_equal(ours@optim$logl, ref@optim$logl, tolerance = 1e-6)
    expect_equal(fitMeasures(ours, fm), fitMeasures(ref, fm), tolerance = 1e-3)
  }
  # A just-identified higher-order factor reproduces the correlated-factor model
  fc <- cfa(three, d, std.lv = TRUE, test = "mean.var.adjusted")
  expect_equal(
    sem(paste(three, "\nG =~ A + C + E"), d, std.lv = TRUE)@optim$logl,
    fc@optim$logl,
    tolerance = 1e-8
  )
  # Nested test of equal regressions,
  f1 <- sem(
    paste(three, "\nE ~ A + C"),
    d,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  f0 <- sem(
    paste(three, "\nE ~ b*A + b*C"),
    d,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  m10 <- lavaan:::lav_test_diff_m10(m1 = f1, m0 = f0, test = FALSE)
  iH <- solve(lavaan::lavTech(m10, "information.observed"))
  J <- lavaan::lavTech(m10, "first.order")
  A <- lavaan:::lav_test_diff_a(m1 = f1, m0 = f0, method = "exact")
  M <- (A %*% iH %*% J %*% iH %*% t(A)) %*% solve(A %*% iH %*% t(A))
  expect_equal(
    lavTestLRT(f1, f0)[2, "Chisq diff"],
    sum(diag(M)) / sum(diag(M %*% M)) * 2 * (f1@optim$logl - f0@optim$logl),
    tolerance = 1e-4
  )
  # Residual covariances vs  model witht them
  fr <- cfa(
    paste(three, "\nA1 ~~ A2"),
    d,
    std.lv = TRUE,
    test = "mean.var.adjusted"
  )
  expect_equal(lavTestLRT(fr, fc)[2, "Df diff"] > 0, TRUE)
  expect_equal(lavTestLRT(fr, fc)[2, "Chisq diff"] > 0, TRUE)
})

# unsupported models are rejected
local({
  expect_error(
    cfa(
      "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v*A\nA ~~ c*C\nc == v",
      bfi,
      std.lv = FALSE
    ),
    "covariances and variances"
  )
  expect_error(
    cfa(
      "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA ~~ v1*A\nC ~~ v2*C\nv1 == v2 + 0.5",
      bfi,
      std.lv = FALSE
    ),
    "variances"
  )
  expect_error(
    cfa("A =~ A1 + a*A2 + A3\nA1 | t*t1\na == t", bfi, std.lv = TRUE),
    "either thresholds"
  )
  expect_error(
    sem("A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3\nA1 ~ C", bfi, std.lv = TRUE),
    "observed variables"
  )
  expect_error(
    cfa("A =~ A1 + A2 + A3", bfi, std.lv = TRUE, group = "gender"),
    "group"
  )
})
