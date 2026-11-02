#### Data and models ####
three_factors <- build_constr_mat(9, 3, "simple")
latent_correlations <- matrix(.3, 3, 3)
diag(latent_correlations) <- 1
items <- sim_data(
  600,
  gen_loadings(three_factors, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
data <- set_latpl_data(items)
regressions <- matrix(0, 3, 3)
regressions[2, 1] <- regressions[3, 1] <- regressions[3, 2] <- NA
residual_covariances <- matrix(0, 9, 9)
residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
residual_covariances[5, 6] <- residual_covariances[6, 5] <- 0.1

#### theta2par() and par2theta() invert each other ####
for (latent_form in c("orthogonal", "lkj", "direct")) {
  for (stdlv in c(TRUE, FALSE)) {
    for (structural in c(FALSE, TRUE)) {
      local({
        model <- set_latpl_model(
          three_factors,
          STDLV = stdlv,
          CORRFLAG = latent_form != "orthogonal",
          LATCOV = if (latent_form == "direct") "direct" else "lkj",
          REGMAT = if (structural) regressions,
          RESCOV = if (structural) residual_covariances
        )
        set.seed(1)
        theta <- runif(length(theta_names(data, model)), -0.3, 0.3)
        par <- theta2par(model, theta)
        expect_equal(par2theta(model, par), theta)
        expect_equal(par$residual_covariances[6, 5], if (structural) 0.1 else 0)
      })
    }
  }
}

#### theta_names() names covariances, regressions and log standard deviations ####
local({
  expect_equal(
    tail(
      theta_names(
        data,
        set_latpl_model(
          three_factors,
          LATCOV = "direct",
          REGMAT = regressions,
          RESCOV = residual_covariances
        )
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
    tail(theta_names(data, set_latpl_model(three_factors, STDLV = FALSE)), 4),
    c("eta[2,3]", "logsd[1]", "logsd[2]", "logsd[3]")
  )
})

#### Two groups split the parameter vector group by group ####
local({
  model <- set_latpl_model(build_constr_mat(9, 3, "simple"))
  data_two_groups <- set_latpl_data(items, GROUP = rep(c("b", "a"), each = 300))
  model_two_groups <- set_latpl_groups(list(model, model))
  names_one_group <- theta_names(data_two_groups, model)
  expect_identical(
    theta_names(data_two_groups, model_two_groups),
    c(names_one_group, paste0(names_one_group, ".g2"))
  )
  set.seed(1)
  theta <- runif(2 * length(names_one_group), -0.3, 0.3)
  par <- theta2par(model_two_groups, theta)
  expect_length(par, 2)
  expect_equal(par[[2]], theta2par(model, theta[-seq_along(names_one_group)]))
  expect_equal(par2theta(model_two_groups, par), theta)
})

#### get_par() returns the raw, listed and transformed parameters ####
local({
  bfi_items <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1
  fit <- fit_latpl(
    set_latpl_data(bfi_items),
    set_latpl_model(build_constr_mat(6, 2, "simple"), STDLV = FALSE)
  )
  expect_identical(get_par(fit, OPTION = "raw"), fit$theta)
  par <- get_par(fit)
  expect_identical(
    names(par),
    c(
      "thresholds",
      "loadings",
      "latent_correlations",
      "regressions",
      "residual_covariances",
      "intercepts",
      "means",
      "scales"
    )
  )
  expect_equal(par, theta2par(fit$model, fit$theta))
  transformed <- get_par(fit, OPTION = "transformed")
  expect_equal(
    transformed[length(transformed) - 2],
    cov2cor(par$latent_correlations)[1, 2]
  )
  expect_equal(tail(transformed, 2), sqrt(diag(par$latent_correlations)))
  expect_stdout(print(fit), "Free parameters")
})

#### LKJ and direct latent covariances give the same fit ####
for (stdlv in c(TRUE, FALSE)) {
  local({
    fit_lkj <- fit_latpl(data, set_latpl_model(three_factors, STDLV = stdlv))
    fit_direct <- fit_latpl(
      data,
      set_latpl_model(three_factors, STDLV = stdlv, LATCOV = "direct")
    )
    vars_lkj <- compute_var(fit_lkj)
    vars_direct <- compute_var(fit_direct)
    test_columns <- c("stat", "df", "stat.scaled", "df.scaled")
    expect_equal(
      fit_direct$numFit$value,
      fit_lkj$numFit$value,
      tolerance = 1e-8
    )
    expect_equal(
      get_par(fit_direct, OPTION = "transformed"),
      get_par(fit_lkj, OPTION = "transformed"),
      tolerance = 1e-4
    )
    expect_equal(
      vars_direct$asymptotic_variance,
      vars_lkj$asymptotic_variance,
      tolerance = 1e-3
    )
    expect_equal(
      compute_plrt(fit_direct, vars_direct)$model[test_columns],
      compute_plrt(fit_lkj, vars_lkj)$model[test_columns],
      tolerance = 1e-3
    )
  })
}

#### A direct latent covariance can be fixed by a constraint ####
local({
  model <- set_latpl_model(three_factors, LATCOV = "direct")
  names_direct <- theta_names(data, model)
  fit <- fit_latpl(
    data,
    set_latpl_model(
      three_factors,
      LATCOV = "direct",
      CONSTRAINTS = list(
        A = matrix(
          replace(
            numeric(length(names_direct)),
            match("psi[1,3]", names_direct),
            1
          ),
          1
        ),
        b = 0
      )
    )
  )
  expect_equal(get_par(fit, OPTION = "list")$latent_correlations[1, 3], 0)
})

#### A saturated structural model is equivalent to correlated factors ####
local({
  fit_correlated <- fit_latpl(
    data,
    set_latpl_model(three_factors, STDLV = FALSE)
  )
  fit_saturated <- fit_latpl(
    data,
    set_latpl_model(
      three_factors,
      STDLV = FALSE,
      CORRFLAG = FALSE,
      REGMAT = regressions
    )
  )
  fit_restricted <- fit_latpl(
    data,
    set_latpl_model(
      three_factors,
      STDLV = FALSE,
      CORRFLAG = FALSE,
      REGMAT = replace(regressions, 3, 0)
    )
  )
  par_saturated <- get_par(fit_saturated, OPTION = "list")
  total_effects <- solve(diag(3) - par_saturated$regressions)
  expect_equal(
    fit_saturated$numFit$value,
    fit_correlated$numFit$value,
    tolerance = 1e-6
  )
  expect_equal(
    total_effects %*% par_saturated$latent_correlations %*% t(total_effects),
    get_par(fit_correlated, OPTION = "list")$latent_correlations,
    tolerance = 1e-3
  )
  expect_equal(
    compute_plrt(fit_saturated)$model$stat,
    compute_plrt(fit_correlated)$model$stat,
    tolerance = 1e-4
  )
  expect_equal(
    compute_plrt(fit_restricted, H1 = fit_correlated)$model[c("stat", "df")],
    compute_plrt(fit_restricted, H1 = fit_saturated)$model[c("stat", "df")],
    tolerance = 1e-4
  )
})

#### Fixed residual covariances keep their value ####
local({
  fixed_and_free <- matrix(0, 9, 9)
  fixed_and_free[1, 2] <- fixed_and_free[2, 1] <- 0.1
  fixed_and_free[4, 7] <- fixed_and_free[7, 4] <- NA
  fit <- fit_latpl(
    data,
    set_latpl_model(
      three_factors,
      REGMAT = regressions,
      RESCOV = fixed_and_free
    )
  )
  expect_equal(get_par(fit, OPTION = "list")$residual_covariances[2, 1], 0.1)
})
