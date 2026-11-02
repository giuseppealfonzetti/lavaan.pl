#### Data and models ####
two_factors <- build_constr_mat(6, 2, "simple")
latent_correlations <- matrix(.3, 2, 2)
diag(latent_correlations) <- 1
items <- sim_data(
  300,
  gen_loadings(two_factors, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
items_missing <- items
set.seed(4)
for (j in c(2:3, 5:6)) {
  items_missing[runif(300) < c(.1, .3, .3, .6)[items[, if (j < 4) 1 else 4] + 1], j] <- NA
}
set.seed(5)
datasets <- list(
  available_cases = set_latpl_data(items_missing, MISSING = "available.cases"),
  weights_and_clusters = set_latpl_data(items, WEIGHTS = runif(300, .2, 3), CLUSTER = rep(1:40, length.out = 300)),
  complete = set_latpl_data(items)
)
fits <- lapply(datasets, function(x) fit_latpl(x, set_latpl_model(two_factors)))
vars <- lapply(fits, compute_var)
tests <- Map(compute_plrt, fits, vars)

#### Helpers ####
saturated <- function(fit, tau, rho, order, scores) {
  engine_data <- list(C_VEC = fit$data$cat, FREQ = fit$data$freq, UFREQ = fit$data$ufreq, Y = fit$data$Y,
                      UW = fit$data$uweights)
  lavaan.pl:::cpp_saturated(tau, rho, engine_data, order, scores)
}

implied_correlations <- function(fit, theta) {
  par <- get_par(replace(fit, "theta", list(theta)), OPTION = "list")
  implied <- par$loadings %*% par$latent_correlations %*% t(par$loadings)
  implied[upper.tri(implied)]
}

#### The saturated fit is a maximum and the degrees of freedom count pairs and thresholds ####
for (case in names(fits)) {
  local({
    fit <- fits[[case]]
    test <- tests[[case]]
    gradient <- saturated(fit, test$saturated$tau, test$saturated$rho, 1, FALSE)$gradient
    expect_true(max(abs(gradient)) / fit$dims$n < 1e-6)
    expect_equal(test$model$df, fit$dims$nthr + fit$dims$pairs - fit$dims$npar)
    expect_equal(test$baseline$df, fit$dims$pairs)
  })
}

#### The mean and variance adjustment matches the dense computation ####
for (case in names(fits)) {
  local({
    fit <- fits[[case]]
    test <- tests[[case]]
    n <- fit$dims$n
    nthr <- fit$dims$nthr
    structural <- (nthr + 1):fit$dims$npar
    at_model <- saturated(fit, fit$theta[1:nthr], implied_correlations(fit, fit$theta), 2, TRUE)
    hessian_saturated <- rbind(cbind(at_model$htt, at_model$htr), cbind(t(at_model$htr), diag(at_model$hrr))) / n
    scores_saturated <- rbind(at_model$st, at_model$sr)
    if (length(fit$data$cluster)) {
      cluster_totals <- t(rowsum(t(scores_saturated) * fit$data$weights, fit$data$cluster))
      cluster_totals <- cluster_totals - rowMeans(cluster_totals)
      variability_saturated <- ncol(cluster_totals) / (ncol(cluster_totals) - 1) * tcrossprod(cluster_totals) / n
    } else {
      variability_saturated <- tcrossprod(scores_saturated) / n
    }
    inverse_saturated <- solve(hessian_saturated)
    correlations_information <- solve(inverse_saturated[-(1:nthr), -(1:nthr)])
    p_saturated <- correlations_information %*%
      (inverse_saturated %*% variability_saturated %*% inverse_saturated)[-(1:nthr), -(1:nthr)]
    inverse_model <- solve(vars[[case]]$H)
    sandwich_model <- (inverse_model %*% vars[[case]]$J %*% inverse_model)[structural, structural]
    p_model <- solve(inverse_model[structural, structural]) %*% sandwich_model
    jacobian <- numDeriv::jacobian(function(x) implied_correlations(fit, x), fit$theta)[, structural]
    mean_statistic <- sum(diag(p_saturated)) - sum(diag(p_model))
    variance_statistic <- 2 * sum(diag(p_saturated %*% p_saturated)) + 2 * sum(diag(p_model %*% p_model)) -
      4 * sum(diag(t(jacobian) %*% correlations_information %*% jacobian %*% sandwich_model %*% p_model))
    expect_equal(test$model$scaling.factor, (variance_statistic / 2) / mean_statistic)
    expect_equal(test$model$df.scaled, mean_statistic^2 / (variance_statistic / 2))
    expect_equal(test$model$stat, 2 * (test$saturated$logl - at_model$ll))
  })
}

#### Closed-form independence thresholds maximise the independence pairwise likelihood ####
for (case in names(fits)) {
  local({
    fit <- fits[[case]]
    test <- tests[[case]]
    gradient <- saturated(fit, test$independence$tau, rep(0, fit$dims$pairs), 1, FALSE)$gradient
    expect_true(max(abs(gradient[1:fit$dims$nthr])) / fit$dims$n < 1e-6)
    expect_equal(test$baseline$stat, 2 * (test$saturated$logl - test$independence$logl))
  })
}

#### Nested test of orthogonal factors against correlated factors ####
local({
  fit_correlated <- fits$complete
  fit_orthogonal <- fit_latpl(datasets$complete, set_latpl_model(two_factors, CORRFLAG = FALSE))
  nested <- compute_plrt(fit_orthogonal, H1 = fit_correlated)$model
  hj <- lavaan.pl:::cpp_HJ(c(fit_orthogonal$theta, 0), fit_correlated$engine$model, fit_correlated$engine$data)
  inverse_hessian <- solve(hj$H)
  correlation <- fit_correlated$dims$npar
  expect_equal(nested$df, 1)
  expect_equal(
    nested$scaling.factor,
    (inverse_hessian %*% hj$J %*% inverse_hessian)[correlation, correlation] / inverse_hessian[correlation, correlation]
  )
  expect_equal(nested$df.scaled, 1)
  expect_equal(nested$stat.scaled, nested$stat / nested$scaling.factor)
  expect_error(compute_plrt(fit_correlated, H1 = fit_orthogonal), "not nested")
})

#### anova() falls back to stats for other models ####
local({
  smaller <- lm(mpg ~ wt, mtcars)
  larger <- lm(mpg ~ wt + hp, mtcars)
  expect_equal(lavaan.pl::anova(smaller, larger), stats::anova(smaller, larger))
})

#### A fixed threshold uses the traces over the full parameter vectors ####
local({
  data <- datasets$complete
  names_free <- theta_names(data, set_latpl_model(two_factors))
  fixed_threshold <- matrix(replace(numeric(length(names_free)), match("tau[2,1]", names_free), 1), 1)
  fit <- fit_latpl(data, set_latpl_model(two_factors, CONSTRAINTS = list(A = fixed_threshold, b = -1)))
  model_vars <- compute_var(fit)
  test <- compute_plrt(fit, model_vars)
  n <- fit$dims$n
  nthr <- fit$dims$nthr
  constraints <- fit$lincon
  correlations_free <- function(x) implied_correlations(fit, as.vector(constraints$K %*% x + constraints$c))
  theta_free <- fit$theta[constraints$free_idx]
  at_model <- lavaan.pl:::cpp_saturated(
    fit$theta[1:nthr],
    correlations_free(theta_free),
    list(C_VEC = fit$data$cat, FREQ = fit$data$freq, UFREQ = numeric(0), Y = fit$data$Y, UW = numeric(0)),
    2,
    TRUE
  )
  hessian_saturated <- rbind(cbind(at_model$htt, at_model$htr), cbind(t(at_model$htr), diag(at_model$hrr))) / n
  p_saturated <- solve(hessian_saturated, tcrossprod(rbind(at_model$st, at_model$sr)) / n)
  p_model <- solve(model_vars$H, model_vars$J)
  jacobian <- rbind(constraints$K[1:nthr, ], numDeriv::jacobian(correlations_free, theta_free))
  mean_statistic <- sum(diag(p_saturated)) - sum(diag(p_model))
  variance_statistic <- 2 * sum(diag(p_saturated %*% p_saturated)) + 2 * sum(diag(p_model %*% p_model)) -
    4 * sum(diag(t(jacobian) %*% hessian_saturated %*% jacobian %*% p_model %*% p_model %*% solve(model_vars$H)))
  expect_equal(test$model$df, nthr + fit$dims$pairs - fit$dims$npar)
  expect_equal(test$model$scaling.factor, (variance_statistic / 2) / mean_statistic)
  expect_equal(test$model$df.scaled, mean_statistic^2 / (variance_statistic / 2))
})

#### Two unconstrained groups: statistic, degrees of freedom and traces add up over groups ####
local({
  data_two_groups <- set_latpl_data(items, GROUP = rep(1:2, each = 150))
  fit_two_groups <- fit_latpl(data_two_groups, set_latpl_groups(list(set_latpl_model(two_factors),
                                                                     set_latpl_model(two_factors))))
  test_separate <- lapply(1:2, function(g) {
    compute_plrt(fit_latpl(set_latpl_data(items[data_two_groups$group == g - 1, ]), set_latpl_model(two_factors)))
  })
  test_two_groups <- compute_plrt(fit_two_groups)
  statistic_and_traces <- function(x) {
    c(x$stat, x$df, x$df.scaled * x$scaling.factor, x$df.scaled * x$scaling.factor^2)
  }
  for (part in c("model", "baseline")) {
    expect_equal(
      statistic_and_traces(test_two_groups[[part]]),
      statistic_and_traces(test_separate[[1]][[part]]) + statistic_and_traces(test_separate[[2]][[part]]),
      tolerance = 1e-6
    )
  }
})

#### Nested test of an equality across groups ####
local({
  data_two_groups <- set_latpl_data(items, GROUP = rep(1:2, each = 150))
  model_two_groups <- set_latpl_groups(list(set_latpl_model(two_factors), set_latpl_model(two_factors)))
  fit_free <- fit_latpl(data_two_groups, model_two_groups)
  names_free <- theta_names(data_two_groups, model_two_groups)
  equal_loadings <- matrix(0, 1, length(names_free))
  equal_loadings[1, match(c("lambda[2,1]", "lambda[2,1].g2"), names_free)] <- c(1, -1)
  fit_equal <- fit_latpl(data_two_groups, set_latpl_groups(model_two_groups$MODELS, list(A = equal_loadings, b = 0)))
  nested <- compute_plrt(fit_equal, H1 = fit_free)$model
  pl <- function(fit) lavaan.pl:::cpp_pl(fit$theta, fit$engine$model, fit$engine$data)$pl
  expect_equal(nested$df, 1)
  expect_equal(nested$stat, 2 * (pl(fit_free) - pl(fit_equal)))
})
