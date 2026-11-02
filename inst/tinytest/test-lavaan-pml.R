if (!at_home()) {
  exit_file("Not run on CRAN")
}

#### Data and models ####
source("fixtures/lavaan-pml-cases.R")
reference <- readRDS("fixtures/lavaan-pml-reference.rds")
ours <- lapply(lavaan_pml_cases, function(case) {
  suppressWarnings(case(lavaan_pl_fitters))
})
ours_nested <- nested_tests(ours, lavTestLRT)

#### Helpers ####
cases_in <- function(group) {
  names(ours)[startsWith(names(ours), paste0(group, "/"))]
}

expect_matches_reference <- function(
  name,
  fit_measures = NULL,
  est_tolerance = 1e-4
) {
  fit <- ours[[name]]
  ref <- reference$cases[[name]]
  expect_equal(
    parameterEstimates(fit)$est,
    ref$estimates$est,
    tolerance = est_tolerance
  )
  expect_equal(parameterEstimates(fit)$se, ref$estimates$se, tolerance = 1e-3)
  expect_equal(fit@optim$logl, ref$logl, tolerance = 1e-6)
  if (length(fit_measures)) {
    expect_equal(
      unclass(fitMeasures(fit, fit_measures)),
      ref$fit_measures[fit_measures],
      tolerance = 1e-3
    )
  }
}

expect_exact_a_difference <- function(fit_h1, fit_h0) {
  m10 <- lavaan:::lav_test_diff_m10(m1 = fit_h1, m0 = fit_h0, test = FALSE)
  inverse_information <- solve(lavaan::lavTech(m10, "information.observed"))
  first_order <- lavaan::lavTech(m10, "first.order")
  restriction <- lavaan:::lav_test_diff_a(
    m1 = fit_h1,
    m0 = fit_h0,
    method = "exact"
  )
  m <- (restriction %*%
    inverse_information %*%
    first_order %*%
    inverse_information %*%
    t(restriction)) %*%
    solve(restriction %*% inverse_information %*% t(restriction))
  difference <- lavTestLRT(fit_h1, fit_h0)
  expect_equal(
    difference[2, "Chisq diff"],
    sum(diag(m)) /
      sum(diag(m %*% m)) *
      2 *
      (fit_h1@optim$logl - fit_h0@optim$logl),
    tolerance = 1e-4
  )
  expect_equal(
    difference[2, "Df diff"],
    sum(diag(m))^2 / sum(diag(m %*% m)),
    tolerance = 1e-4
  )
}

fit_measures_full <- c(
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
fit_measures_short <- c(
  "chisq",
  "chisq.scaled",
  "df.scaled",
  "baseline.chisq.scaled",
  "cfi.scaled"
)

#### The references cover every case ####
expect_identical(names(reference$cases), names(lavaan_pml_cases))

#### Estimates, standard errors, covariance matrix and log-likelihood ####
for (name in cases_in("grid")) {
  local({
    fit <- ours[[name]]
    ref <- reference$cases[[name]]
    expect_identical(as.character(class(fit)), "lavaan")
    expect_true(lavaan::lavInspect(fit, "converged"))
    expect_identical(parameterEstimates(fit)[, 1:3], ref$estimates[, 1:3])
    expect_matches_reference(name)
    expect_equal(unclass(vcov(fit)), ref$vcov, tolerance = 1e-4)
    expect_equal(fitted(fit), ref$fitted, tolerance = 1e-4)
  })
}

#### The inverse Hessian matches lavaan's observed information ####
local({
  fit <- ours[["grid/binary, one factor"]]
  observed_information <- lavaan:::lav_model_info_observed(
    lavmodel = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavdata = fit@Data,
    lavoptions = fit@Options,
    lavcache = fit@Cache
  )
  expect_equal(
    attr(fit@vcov$vcov, "E.inv"),
    solve(observed_information),
    tolerance = 1e-5
  )
})

#### Test statistics and fit measures ####
for (name in cases_in("tests")) {
  local({
    fit <- ours[[name]]
    ref <- reference$cases[[name]]
    expect_equal(
      fit@test[[1]][c("stat", "df")],
      ref$test[[1]][c("stat", "df")],
      tolerance = 1e-3
    )
    expect_equal(
      fit@test[[2]][c("stat", "df", "pvalue", "scaling.factor")],
      ref$test[[2]][c("stat", "df", "pvalue", "scaling.factor")],
      tolerance = 1e-3
    )
    expect_equal(
      unclass(fitMeasures(fit, fit_measures_full)),
      ref$fit_measures[fit_measures_full],
      tolerance = 1e-3
    )
  })
}

#### Nested tests match lavaan's information and its exact restriction ####
local({
  for (comparison in c("orthogonal", "equal_loadings")) {
    expect_equal(
      ours_nested[[comparison]],
      reference$nested[[comparison]],
      tolerance = 1e-6
    )
  }
  expect_equal(
    lavaan.pl::anova(ours[["nested/correlated"]], ours[["nested/orthogonal"]]),
    lavTestLRT(ours[["nested/correlated"]], ours[["nested/orthogonal"]])
  )
  expect_exact_a_difference(
    ours[["nested/correlated"]],
    ours[["nested/orthogonal"]]
  )
  expect_exact_a_difference(
    ours[["nested/correlated"]],
    ours[["nested/equal loadings"]]
  )
})

#### Missing data ####
local({
  fit_measures <- c(
    "chisq",
    "chisq.scaled",
    "df.scaled",
    "pvalue.scaled",
    "baseline.chisq.scaled",
    "cfi.scaled",
    "tli.scaled"
  )
  for (missing in c("listwise", "pairwise", "available.cases")) {
    name <- paste0("missing/", missing)
    expect_equal(
      lavaan::lavInspect(ours[[name]], "nobs"),
      reference$cases[[name]]$nobs
    )
    expect_matches_reference(name, fit_measures)
  }
  for (missing in c("pairwise", "available.cases")) {
    expect_exact_a_difference(
      ours[[paste0("missing/", missing)]],
      ours[[paste0("missing/", missing, ", orthogonal")]]
    )
  }
})

#### Sampling weights and clusters ####
local({
  for (name in cases_in("weights")) {
    expect_matches_reference(name, fit_measures_short)
  }
  for (name in cases_in("clusters")) {
    expect_matches_reference(name)
  }
})

#### Linear equality constraints ####
for (name in cases_in("constraints")) {
  expect_matches_reference(name, fit_measures_short)
}

#### Equal thresholds give a stationary point of lavaan's objective ####
local({
  fit <- cfa(
    paste(two_factors_four_items, "\nA3 | t*t3\nA4 | t*t3"),
    bfi_600,
    std.lv = TRUE
  )
  gradient <- lavaan::lavTech(fit, "gradient")
  constraints_jacobian <- fit@Model@ceq.JAC
  null_space <- qr.Q(qr(t(constraints_jacobian)), complete = TRUE)[,
    -seq_len(nrow(constraints_jacobian)),
    drop = FALSE
  ]
  expect_true(
    max(abs(crossprod(null_space, gradient))) < 1e-4 * max(abs(gradient))
  )
})

#### Fixed and constrained latent covariances ####
for (name in setdiff(
  cases_in("latent covariances"),
  "latent covariances/free"
)) {
  expect_identical(ours[[name]]@external$lavaan.pl$model$LATCOV, "direct")
  expect_matches_reference(name, fit_measures_short)
}

#### A bifactor model reaches at least lavaan's log-likelihood ####
expect_true(
  ours[["bifactor"]]@optim$logl >= reference$cases[["bifactor"]]$logl - 1e-6
)

#### Nested test of a zero latent covariance ####
expect_exact_a_difference(
  ours[["latent covariances/free"]],
  ours[["latent covariances/zero"]]
)

#### Residual covariances, structural regressions and higher-order factors ####
for (name in cases_in("structural")) {
  expect_matches_reference(name, fit_measures_short)
}

#### Nested test of equal regressions ####
expect_exact_a_difference(
  ours[["structural/two regressions"]],
  ours[["structural/equal regressions"]]
)

#### Multigroup models ####
for (name in cases_in("multigroup")) {
  local({
    fit <- ours[[name]]
    ref <- reference$cases[[name]]
    fit_measures <- c(
      "chisq",
      "df",
      "chisq.scaled",
      "df.scaled",
      "pvalue.scaled",
      "chisq.scaling.factor",
      "baseline.chisq.scaled",
      "cfi.scaled",
      "tli.scaled"
    )
    expect_identical(
      parameterEstimates(fit)[, c("lhs", "op", "rhs", "group")],
      ref$estimates[, c("lhs", "op", "rhs", "group")]
    )
    expect_matches_reference(name, fit_measures, est_tolerance = 1e-3)
    expect_equal(
      lapply(fit@test, function(x) x[c("stat", "stat.group", "df")]),
      lapply(ref$test, function(x) x[c("stat", "stat.group", "df")]),
      tolerance = 1e-3
    )
  })
}

#### The invariance chain ####
local({
  for (parameterization in c("delta", "theta")) {
    level <- function(name) {
      ours[[paste0("multigroup/", parameterization, "/", name)]]
    }
    chain <- lavTestLRT(
      level("configural"),
      level("thresholds"),
      level("thresholds_loadings"),
      level("thresholds_loadings_intercepts")
    )
    comparison <- paste0("invariance_", parameterization)
    expect_equal(
      ours_nested[[comparison]],
      reference$nested[[comparison]],
      tolerance = 1e-6
    )
    expect_true(all(chain[-1, "Chisq diff"] > 0 & chain[-1, "Df diff"] > 0))
    expect_error(
      lavTestLRT(level("loadings"), level("thresholds_loadings")),
      "not nested"
    )
  }
  expect_exact_a_difference(
    ours[["multigroup/theta/configural"]],
    ours[["multigroup/theta/loadings"]]
  )
})
