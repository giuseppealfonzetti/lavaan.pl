#### Fitting functions ####
lavaan_pl_fitters <- list(cfa = lavaan.pl::cfa, sem = lavaan.pl::sem)
lavaan_pml_fitters <- list(
  cfa = function(...) lavaan::cfa(..., estimator = "PML"),
  sem = function(...) lavaan::sem(..., estimator = "PML")
)

#### Data ####
bfi_sample <- function(items, size) {
  data <- as.data.frame(bfi)
  for (item in items) {
    data[[item]] <- ordered(data[[item]])
  }
  set.seed(1)
  data[complete.cases(data[, items]), ][sample(2000, size), ]
}

eight_items <- c("A1", "A2", "A3", "A4", "C1", "C2", "C3", "C4")
twelve_items <- c(eight_items, "E1", "E2", "E3", "E4")

bfi_weights <- bfi_sample(eight_items, 400)
bfi_weights$w <- runif(400, .2, 3)
bfi_weights$cl <- rep(1:40, length.out = 400)
bfi_weights_missing <- bfi_weights
for (item in c("A2", "A3", "C2", "C3")) {
  bfi_weights_missing[[item]][runif(400) < .15] <- NA
}
bfi_600 <- bfi_sample(eight_items, 600)
bfi_800 <- bfi_sample(twelve_items, 800)

simulated_loadings <- build_constr_mat(12, 2, "simple")
simulated_loadings[is.na(simulated_loadings)] <- c(
  .4,
  .5,
  .6,
  .7,
  .8,
  .9,
  .4,
  .5,
  .6,
  .7,
  .8,
  .9
)
simulated_complete <- sim_data(
  1000,
  simulated_loadings,
  c(-1.5, .5, 1.5),
  matrix(c(1, .6, .6, 1), 2),
  SEED = 4
)
simulated_missing <- simulated_listwise <- simulated_complete
set.seed(5)
missing_probability <- c(.119, .269, .269, .731)
for (j in 2:6) {
  simulated_missing[
    runif(1000) < missing_probability[simulated_complete[, 1] + 1],
    j
  ] <- NA
}
for (j in 8:12) {
  simulated_missing[
    runif(1000) < missing_probability[simulated_complete[, 7] + 1],
    j
  ] <- NA
}
simulated_listwise[, c(2, 8)] <- simulated_missing[, c(2, 8)]
colnames(simulated_missing) <- colnames(simulated_listwise) <- paste0("y", 1:12)
simulated_missing <- as.data.frame(lapply(
  as.data.frame(simulated_missing),
  ordered
))
simulated_listwise <- as.data.frame(lapply(
  as.data.frame(simulated_listwise),
  ordered
))

#### Models ####
two_factors_four_items <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4"
two_factors_five_items <- "A =~ A1 + A2 + A3 + A4 + A5\nC =~ C1 + C2 + C3 + C4 + C5"
three_factors_four_items <- "A =~ A1 + A2 + A3 + A4\nC =~ C1 + C2 + C3 + C4\nE =~ E1 + E2 + E3 + E4"
simulated_model <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6\nf2 =~ y7 + y8 + y9 + y10 + y11 + y12"
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
grid_with_tests <- c(
  "binary, one factor",
  "ordinal, std.lv",
  "ordinal, marker",
  "cross-loading",
  "equality constraints",
  "theta, std.lv",
  "theta, marker, residual covariance",
  "latent mean"
)
group_equal <- list(
  configural = character(0),
  loadings = "loadings",
  thresholds = "thresholds",
  thresholds_loadings = c("thresholds", "loadings"),
  thresholds_loadings_intercepts = c("thresholds", "loadings", "intercepts")
)

#### Cases ####
grid_case <- function(case) {
  parameterization <- if (is.null(case$parameterization)) {
    "delta"
  } else {
    case$parameterization
  }
  function(fit) {
    fit$cfa(
      case$model,
      case$data,
      std.lv = case$std.lv,
      parameterization = parameterization
    )
  }
}

grid_test_case <- function(case) {
  parameterization <- if (is.null(case$parameterization)) {
    "delta"
  } else {
    case$parameterization
  }
  function(fit) {
    fit$cfa(
      case$model,
      case$data,
      std.lv = case$std.lv,
      parameterization = parameterization,
      test = "mean.var.adjusted"
    )
  }
}

multigroup_cases <- function(parameterization) {
  setNames(
    lapply(group_equal, function(equal) {
      function(fit) {
        fit$cfa(
          "A =~ A1 + A2 + A3\nC =~ C1 + C2 + C3",
          bfi[1:1000, ],
          group = "gender",
          group.equal = equal,
          std.lv = parameterization == "delta",
          parameterization = parameterization,
          test = "mean.var.adjusted"
        )
      }
    }),
    paste0("multigroup/", parameterization, "/", names(group_equal))
  )
}

lavaan_pml_cases <- c(
  setNames(lapply(grid, grid_case), paste0("grid/", names(grid))),
  setNames(
    lapply(grid[grid_with_tests], grid_test_case),
    paste0("tests/", grid_with_tests)
  ),
  list(
    "nested/correlated" = function(fit) {
      fit$cfa(
        two_factors_five_items,
        bfi,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "nested/orthogonal" = function(fit) {
      fit$cfa(
        two_factors_five_items,
        bfi,
        std.lv = TRUE,
        orthogonal = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "nested/equal loadings" = function(fit) {
      fit$cfa(
        "A =~ A1 + a*A2 + a*A3 + A4 + A5\nC =~ C1 + b*C2 + b*C3 + C4 + C5",
        bfi,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "missing/listwise" = function(fit) {
      fit$cfa(
        simulated_model,
        simulated_listwise,
        std.lv = TRUE,
        missing = "listwise",
        test = "mean.var.adjusted"
      )
    },
    "missing/pairwise" = function(fit) {
      fit$cfa(
        simulated_model,
        simulated_missing,
        std.lv = TRUE,
        missing = "pairwise",
        test = "mean.var.adjusted"
      )
    },
    "missing/available.cases" = function(fit) {
      fit$cfa(
        simulated_model,
        simulated_missing,
        std.lv = TRUE,
        missing = "available.cases",
        test = "mean.var.adjusted"
      )
    },
    "missing/pairwise, orthogonal" = function(fit) {
      fit$cfa(
        simulated_model,
        simulated_missing,
        std.lv = TRUE,
        orthogonal = TRUE,
        missing = "pairwise",
        test = "mean.var.adjusted"
      )
    },
    "missing/available.cases, orthogonal" = function(fit) {
      fit$cfa(
        simulated_model,
        simulated_missing,
        std.lv = TRUE,
        orthogonal = TRUE,
        missing = "available.cases",
        test = "mean.var.adjusted"
      )
    },
    "weights/design" = function(fit) {
      fit$cfa(
        two_factors_four_items,
        bfi_weights,
        std.lv = TRUE,
        test = "mean.var.adjusted",
        sampling.weights = "w"
      )
    },
    "weights/frequency" = function(fit) {
      fit$cfa(
        two_factors_four_items,
        bfi_weights,
        std.lv = TRUE,
        test = "mean.var.adjusted",
        sampling.weights = "w",
        sampling.weights.type = "frequency"
      )
    },
    "weights/missing" = function(fit) {
      fit$cfa(
        two_factors_four_items,
        bfi_weights_missing,
        std.lv = TRUE,
        test = "mean.var.adjusted",
        sampling.weights = "w",
        missing = "pairwise"
      )
    },
    "clusters/clusters" = function(fit) {
      fit$cfa(
        two_factors_four_items,
        bfi_weights,
        std.lv = TRUE,
        test = "none",
        cluster = "cl"
      )
    },
    "clusters/weights and clusters" = function(fit) {
      fit$cfa(
        two_factors_four_items,
        bfi_weights,
        std.lv = TRUE,
        test = "none",
        sampling.weights = "w",
        cluster = "cl"
      )
    },
    "constraints/linear" = function(fit) {
      fit$cfa(
        "A =~ A1 + a*A2 + b*A3 + A4\nC =~ C1 + C2 + C3 + C4\na == 0.5*b + 0.1",
        bfi_600,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "constraints/fixed threshold" = function(fit) {
      fit$cfa(
        paste(two_factors_four_items, "\nA1 | -0.4*t1"),
        bfi_600,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "constraints/equal variances" = function(fit) {
      fit$cfa(
        paste(two_factors_four_items, "\nA ~~ v*A\nC ~~ v*C"),
        bfi_600,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "constraints/variance ratio" = function(fit) {
      fit$cfa(
        paste(two_factors_four_items, "\nA ~~ v1*A\nC ~~ v2*C\nv1 == 2*v2"),
        bfi_600,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/free" = function(fit) {
      fit$cfa(
        three_factors_four_items,
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/zero" = function(fit) {
      fit$cfa(
        paste(three_factors_four_items, "\nA ~~ 0*E"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/fixed" = function(fit) {
      fit$cfa(
        paste(three_factors_four_items, "\nA ~~ 0.3*C"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/equal" = function(fit) {
      fit$cfa(
        paste(three_factors_four_items, "\nA ~~ r*C\nA ~~ r*E"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/fixed, marker" = function(fit) {
      fit$cfa(
        paste(three_factors_four_items, "\nA ~~ 0.1*E"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "latent covariances/equal, marker" = function(fit) {
      fit$cfa(
        paste(three_factors_four_items, "\nA ~~ r*C\nA ~~ r*E"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "bifactor" = function(fit) {
      fit$cfa(
        paste(
          three_factors_four_items,
          "\nG =~ A1 + A2 + C1 + C2 + E1 + E2\nG ~~ 0*A\nG ~~ 0*C\nG ~~ 0*E"
        ),
        bfi_800,
        std.lv = TRUE
      )
    },
    "structural/residual covariances" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nA1 ~~ A2\nC1 ~~ E1"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "structural/residual covariances, marker" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nA1 ~~ A2\nC1 ~~ E1"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "structural/regressions" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nC ~ A\nE ~ A + C"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "structural/regressions, marker" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nC ~ A\nE ~ A + C"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "structural/two regressions, marker" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nE ~ A + C"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    },
    "structural/two regressions" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nE ~ A + C"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "structural/equal regressions" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nE ~ b*A + b*C"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "structural/regressions and residual covariance" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nC ~ A\nE ~ A + C\nA1 ~~ A2"),
        bfi_800,
        std.lv = TRUE,
        test = "mean.var.adjusted"
      )
    },
    "structural/higher order, marker" = function(fit) {
      fit$sem(
        paste(three_factors_four_items, "\nG =~ A + C + E"),
        bfi_800,
        std.lv = FALSE,
        test = "mean.var.adjusted"
      )
    }
  ),
  multigroup_cases("delta"),
  multigroup_cases("theta")
)

#### Summaries ####
pml_summary <- function(fit, with_vcov) {
  estimates <- parameterEstimates(fit)
  list(
    estimates = estimates[, intersect(
      c("lhs", "op", "rhs", "group", "est", "se"),
      names(estimates)
    )],
    vcov = if (with_vcov) unclass(vcov(fit)),
    fitted = fitted(fit),
    logl = fit@optim$logl,
    nobs = lavaan::lavInspect(fit, "nobs"),
    test = lapply(fit@test, function(x) {
      x[intersect(
        names(x),
        c("stat", "stat.group", "df", "pvalue", "scaling.factor")
      )]
    }),
    fit_measures = if (!identical(fit@Options$test, "none")) {
      unclass(fitMeasures(fit))
    }
  )
}

nested_tests <- function(fits, lrt) {
  columns <- c("Df", "PL_AIC", "PL_BIC", "Chisq")
  chain <- function(parameterization) {
    level <- function(name) {
      fits[[paste0("multigroup/", parameterization, "/", name)]]
    }
    lrt(
      level("configural"),
      level("thresholds"),
      level("thresholds_loadings"),
      level("thresholds_loadings_intercepts")
    )
  }
  lapply(
    list(
      orthogonal = lrt(
        fits[["nested/correlated"]],
        fits[["nested/orthogonal"]]
      ),
      equal_loadings = lrt(
        fits[["nested/correlated"]],
        fits[["nested/equal loadings"]]
      ),
      invariance_delta = chain("delta"),
      invariance_theta = chain("theta")
    ),
    function(x) unlist(x[, columns])
  )
}
