source("fixtures/simulate.R")
#### Data and models ####
three_factors <- build_constr_mat(9, 3, "simple")
three_markers <- replace(three_factors, cbind(c(1, 4, 7), 1:3), 1)
latent_correlations <- matrix(.3, 3, 3)
diag(latent_correlations) <- 1
items <- sim_data(
  600,
  gen_loadings(three_factors, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
data <- lavaan.pl:::lavaanpl_data(items)
regressions <- matrix(0, 3, 3)
regressions[2, 1] <- regressions[3, 1] <- regressions[3, 2] <- NA
residual_covariances <- matrix(0, 9, 9)
residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
residual_covariances[5, 6] <- residual_covariances[6, 5] <- 0.1

#### lavaan.pl:::vector_to_matrices() and lavaan.pl:::matrices_to_vector() invert each other ####
for (latent_form in c("orthogonal", "lkj", "direct")) {
  for (stdlv in c(TRUE, FALSE)) {
    for (structural in c(FALSE, TRUE)) {
      local({
        model <- lavaan.pl:::lavaanpl_model(
          if (stdlv) three_factors else three_markers,
          PSI = replace(matrix(if (latent_form == "orthogonal") 0 else NA, 3, 3), cbind(1:3, 1:3), if (stdlv) 1 else NA),
          LATCOV = if (latent_form == "direct") "direct" else "lkj",
          BETA = if (structural) regressions,
          THETA = if (structural) residual_covariances
        )
        set.seed(1)
        theta <- runif(length(lavaan.pl:::parameter_names(data, model)), -0.3, 0.3)
        par <- lavaan.pl:::vector_to_matrices(model, theta)
        expect_equal(lavaan.pl:::matrices_to_vector(model, par), theta)
        expect_equal(par$theta[6, 5], if (structural) 0.1 else 0)
      })
    }
  }
}

#### lavaan.pl:::parameter_names() names covariances, regressions and log standard deviations ####
local({
  expect_equal(
    tail(
      lavaan.pl:::parameter_names(
        data,
        lavaan.pl:::lavaanpl_model(
          three_factors,
          LATCOV = "direct",
          BETA = regressions,
          THETA = residual_covariances
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
    tail(lavaan.pl:::parameter_names(data, lavaan.pl:::lavaanpl_model(three_markers, PSI = matrix(NA, 3, 3))), 4),
    c("eta[2,3]", "logsd[1]", "logsd[2]", "logsd[3]")
  )
})

#### Two groups split the parameter vector group by group ####
local({
  model <- lavaan.pl:::lavaanpl_model(build_constr_mat(9, 3, "simple"))
  data_two_groups <- lavaan.pl:::lavaanpl_data(items, GROUP = rep(c("b", "a"), each = 300))
  model_two_groups <- lavaan.pl:::lavaanpl_groups(list(model, model))
  names_one_group <- lavaan.pl:::parameter_names(data_two_groups, model)
  expect_identical(
    lavaan.pl:::parameter_names(data_two_groups, model_two_groups),
    c(names_one_group, paste0(names_one_group, ".g2"))
  )
  set.seed(1)
  theta <- runif(2 * length(names_one_group), -0.3, 0.3)
  par <- lavaan.pl:::vector_to_matrices(model_two_groups, theta)
  expect_length(par, 2)
  expect_equal(par[[2]], lavaan.pl:::vector_to_matrices(model, theta[-seq_along(names_one_group)]))
  expect_equal(lavaan.pl:::matrices_to_vector(model_two_groups, par), theta)
})

#### lavaan.pl:::vector_to_matrices() and lavaan.pl:::transform_parameters() return the model matrices and the interpretable parameters ####
local({
  bfi_items <- as.matrix(data.frame(lapply(bfi[, 1:6], as.numeric))) - 1
  fit <- lavaan.pl:::fit_lavaanpl(
    lavaan.pl:::lavaanpl_data(bfi_items),
    lavaan.pl:::lavaanpl_model(replace(build_constr_mat(6, 2, "simple"), cbind(c(1, 4), 1:2), 1), PSI = matrix(NA, 2, 2)),
    CONTROL = set_control("QN")
  )
  par <- lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]
  expect_identical(
    names(par),
    c("tau", "lambda", "psi", "beta", "theta", "nu", "alpha", "delta")
  )
  transformed <- lavaan.pl:::transform_parameters(fit$model, fit$parameter_vector)
  expect_equal(
    transformed[length(transformed) - 2],
    cov2cor(par$psi)[1, 2]
  )
  expect_equal(tail(transformed, 2), sqrt(diag(par$psi)))
  expect_stdout(print(fit), "Free parameters")
})

#### LKJ and direct latent covariances give the same fit ####
for (stdlv in c(TRUE, FALSE)) {
  local({
    loadings <- if (stdlv) three_factors else three_markers
    psi <- if (!stdlv) matrix(NA, 3, 3)
    fit_lkj <- lavaan.pl:::fit_lavaanpl(data, lavaan.pl:::lavaanpl_model(loadings, psi), CONTROL = set_control("QN"))
    fit_direct <- lavaan.pl:::fit_lavaanpl(
      data,
      lavaan.pl:::lavaanpl_model(loadings, psi, LATCOV = "direct"),
      CONTROL = set_control("QN")
    )
    vars_lkj <- lavaan.pl:::compute_var(fit_lkj)
    vars_direct <- lavaan.pl:::compute_var(fit_direct)
    test_columns <- c("stat", "df", "stat.scaled", "df.scaled")
    expect_equal(
      fit_direct$numFit$value,
      fit_lkj$numFit$value,
      tolerance = 1e-8
    )
    expect_equal(
      lavaan.pl:::transform_parameters(fit_direct$model, fit_direct$parameter_vector),
      lavaan.pl:::transform_parameters(fit_lkj$model, fit_lkj$parameter_vector),
      tolerance = 1e-4
    )
    expect_equal(
      vars_direct$asymptotic_variance,
      vars_lkj$asymptotic_variance,
      tolerance = 1e-3
    )
    expect_equal(
      lavaan.pl:::compute_plrt(fit_direct, vars_direct)$model[test_columns],
      lavaan.pl:::compute_plrt(fit_lkj, vars_lkj)$model[test_columns],
      tolerance = 1e-3
    )
  })
}

#### A direct latent covariance can be fixed by a constraint ####
local({
  model <- lavaan.pl:::lavaanpl_model(three_factors, LATCOV = "direct")
  names_direct <- lavaan.pl:::parameter_names(data, model)
  fit <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(
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
    ),
    CONTROL = set_control("QN")
  )
  expect_equal(lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]$psi[1, 3], 0)
})

#### A saturated structural model is equivalent to correlated factors ####
local({
  fit_correlated <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(three_markers, PSI = matrix(NA, 3, 3)),
    CONTROL = set_control("QN")
  )
  fit_saturated <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(
      three_markers,
      PSI = diag(NA, 3),
      BETA = regressions
    ),
    CONTROL = set_control("QN")
  )
  fit_restricted <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(
      three_markers,
      PSI = diag(NA, 3),
      BETA = replace(regressions, 3, 0)
    ),
    CONTROL = set_control("QN")
  )
  par_saturated <- lavaan.pl:::vector_to_matrices(fit_saturated$model, fit_saturated$parameter_vector)[[1]]
  total_effects <- solve(diag(3) - par_saturated$beta)
  expect_equal(
    fit_saturated$numFit$value,
    fit_correlated$numFit$value,
    tolerance = 1e-6
  )
  expect_equal(
    total_effects %*% par_saturated$psi %*% t(total_effects),
    lavaan.pl:::vector_to_matrices(fit_correlated$model, fit_correlated$parameter_vector)[[1]]$psi,
    tolerance = 1e-3
  )
  expect_equal(
    lavaan.pl:::compute_plrt(fit_saturated)$model$stat,
    lavaan.pl:::compute_plrt(fit_correlated)$model$stat,
    tolerance = 1e-4
  )
  expect_equal(
    lavaan.pl:::compute_plrt(fit_restricted, H1 = fit_correlated)$model[c("stat", "df")],
    lavaan.pl:::compute_plrt(fit_restricted, H1 = fit_saturated)$model[c("stat", "df")],
    tolerance = 1e-4
  )
})

#### Fixed residual covariances keep their value ####
local({
  fixed_and_free <- matrix(0, 9, 9)
  fixed_and_free[1, 2] <- fixed_and_free[2, 1] <- 0.1
  fixed_and_free[4, 7] <- fixed_and_free[7, 4] <- NA
  fit <- lavaan.pl:::fit_lavaanpl(
    data,
    lavaan.pl:::lavaanpl_model(
      three_factors,
      BETA = regressions,
      THETA = fixed_and_free
    ),
    CONTROL = set_control("QN")
  )
  expect_equal(lavaan.pl:::vector_to_matrices(fit$model, fit$parameter_vector)[[1]]$theta[2, 1], 0.1)
})
