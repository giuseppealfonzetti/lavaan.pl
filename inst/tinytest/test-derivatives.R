#### Data and models ####
loadings_structure <- build_constr_mat(9, 3, "simple")
loadings_structure[4, 1] <- NA
latent_correlations <- matrix(.3, 3, 3)
diag(latent_correlations) <- 1
items <- sim_data(
  300,
  gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
items[, 1:3] <- 1 * (items[, 1:3] >= 2)
items_missing <- items
set.seed(4)
items_missing[sample(length(items_missing), 150)] <- NA
data_missing <- set_latpl_data(items_missing, MISSING = "available.cases")
regressions <- matrix(0, 3, 3)
regressions[2, 1] <- regressions[3, 1] <- regressions[3, 2] <- NA
residual_covariances <- matrix(0, 9, 9)
residual_covariances[1, 2] <- residual_covariances[2, 1] <- NA
residual_covariances[5, 6] <- residual_covariances[6, 5] <- 0.1

#### Helpers ####
unit_scores <- function(fit, theta) {
  data <- fit$data
  t(vapply(
    seq_len(data$n),
    function(i) {
      response <- data$Y[i, , drop = FALSE]
      unit_data <- lapply(fit$engine$data, function(x) {
        x$FREQ[] <- 0
        x$UFREQ[] <- 0
        x
      })
      group <- data$group[i] + 1
      unit_data[[group]]$FREQ <- pairs_freq(response, data$cat)
      unit_data[[group]]$UFREQ <- unlist(lapply(seq_len(data$p), function(j) {
        data$uweights[i] *
          (seq_len(data$cat[j]) - 1 == response[j] & !is.na(response[j]))
      }))
      lavaan.pl:::cpp_pl(theta, fit$engine$model, unit_data, 1L)$gradient
    },
    numeric(length(theta))
  ))
}

expect_derivatives <- function(fit) {
  block <- lavaan.pl:::theta_block(fit$model, sum(fit$data$cat - 1))
  size <- function(name) sum(block == name)
  theta <- numeric(length(block))
  theta[block == "tau"] <- unlist(rep(
    lapply(fit$data$cat - 1, function(m) seq(-0.8, 0.8, length.out = m)),
    fit$data$ngroups
  ))
  theta[block == "lambda"] <- seq(.4, .6, length.out = size("lambda"))
  theta[block == "corr"] <- seq(.05, .15, length.out = size("corr"))
  theta[block == "logsd"] <- -.5
  theta[block == "beta"] <- seq(.2, .3, length.out = size("beta"))
  theta[block == "rescov"] <- .05
  theta[block == "nu"] <- seq(-.2, .2, length.out = size("nu"))
  theta[block == "alpha"] <- seq(.3, -.2, length.out = size("alpha"))
  theta[block == "scale"] <- seq(-.1, .1, length.out = size("scale"))

  pl <- function(x, order) {
    lavaan.pl:::cpp_pl(x, fit$engine$model, fit$engine$data, order)
  }
  rho <- function(x) lavaan.pl:::cpp_rho(x, fit$engine$model, fit$engine$data)
  gradient <- pl(theta, 1L)$gradient
  hj <- lavaan.pl:::cpp_HJ(theta, fit$engine$model, fit$engine$data)
  hessian_numeric <- -numDeriv::jacobian(
    function(x) pl(x, 1L)$gradient,
    theta
  ) /
    fit$dims$n
  scores <- unit_scores(fit, theta)
  expect_equal(
    gradient,
    numDeriv::grad(function(x) pl(x, 0L)$pl, theta),
    tolerance = 1e-6
  )
  expect_equal(
    hj$H,
    (hessian_numeric + t(hessian_numeric)) / 2,
    tolerance = 1e-6
  )
  expect_equal(hj$gradient, gradient)
  expect_equal(colSums(scores), gradient)
  expect_equal(hj$J, crossprod(scores) / fit$dims$n)
  expect_equal(
    rho(theta)$jacobian,
    numDeriv::jacobian(function(x) rho(x)$rho, theta),
    tolerance = 1e-6
  )
  expect_equal(
    rho(theta)$tau_jacobian,
    numDeriv::jacobian(function(x) rho(x)$tau, theta),
    tolerance = 1e-6
  )
}

engine_only <- function(data, model) {
  suppressWarnings(fit_latpl(
    data,
    model,
    CONTROL = set_ucminf_control(MAXEVAL = 1)
  ))
}

#### Gradient, H, J and Jacobians of the latent covariance forms and structural models ####
for (latent_form in c("orthogonal", "lkj", "direct")) {
  for (stdlv in c(TRUE, FALSE)) {
    for (structural in c(FALSE, TRUE)) {
      expect_derivatives(engine_only(
        data_missing,
        set_latpl_model(
          loadings_structure,
          STDLV = stdlv,
          CORRFLAG = latent_form != "orthogonal",
          LATCOV = if (latent_form == "direct") "direct" else "lkj",
          REGMAT = if (structural) regressions,
          RESCOV = if (structural) residual_covariances
        )
      ))
    }
  }
}

#### Gradient, H, J and Jacobians with intercepts, latent means, scales and the theta parameterization ####
for (parameterization in c("delta", "theta")) {
  for (latent_form in c("lkj", "direct")) {
    for (structural in c(FALSE, TRUE)) {
      local({
        rescov <- if (structural) residual_covariances else diag(0, 9)
        if (parameterization == "theta") {
          diag(rescov) <- c(NA, 1.2, NA, 1, NA, 0.8, NA, 1, NA)
        }
        expect_derivatives(engine_only(
          data_missing,
          set_latpl_model(
            loadings_structure,
            STDLV = FALSE,
            LATCOV = latent_form,
            REGMAT = if (structural) regressions,
            RESCOV = if (structural || parameterization == "theta") rescov,
            PARAMETERIZATION = parameterization,
            INTERCEPTS = c(NA, 0.1, NA, NA, 0, NA, -0.1, NA, NA),
            MEANS = c(NA, 0.2, NA),
            SCALES = if (parameterization == "delta") {
              c(NA, 1.1, NA, NA, 1, NA, 0.9, NA, NA)
            }
          )
        ))
      })
    }
  }
}

#### Gradient, H, J and Jacobians with two groups, the second with a mean structure ####
local({
  data_two_groups <- set_latpl_data(
    items_missing,
    MISSING = "available.cases",
    GROUP = rep(c("a", "b"), length.out = nrow(items_missing))
  )
  model_two_groups <- set_latpl_groups(list(
    set_latpl_model(loadings_structure, LATCOV = "direct"),
    set_latpl_model(
      loadings_structure,
      STDLV = FALSE,
      PARAMETERIZATION = "theta",
      RESCOV = diag(c(NA, 1, NA, 1, 1, 1, 1, 1, 1)),
      INTERCEPTS = c(NA, 0.1, rep(0, 7)),
      MEANS = c(NA, 0, 0.2)
    )
  ))
  expect_derivatives(engine_only(data_two_groups, model_two_groups))
})

#### Frequency weights give the H and J of replicated rows ####
local({
  items_complete <- items[, 4:9]
  model <- set_latpl_model(build_constr_mat(6, 2, "simple"))
  set.seed(2)
  weights <- sample(1:3, nrow(items_complete), replace = TRUE)
  fit_weighted <- fit_latpl(
    set_latpl_data(
      items_complete,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency"
    ),
    model
  )
  fit_replicated <- fit_latpl(
    set_latpl_data(items_complete[
      rep(seq_len(nrow(items_complete)), weights),
    ]),
    model
  )
  theta <- fit_weighted$theta + 0.05
  hj_weighted <- lavaan.pl:::cpp_HJ(
    theta,
    fit_weighted$engine$model,
    fit_weighted$engine$data
  )
  hj_replicated <- lavaan.pl:::cpp_HJ(
    theta,
    fit_replicated$engine$model,
    fit_replicated$engine$data
  )
  expect_equal(hj_weighted$H, hj_replicated$H)
  expect_equal(hj_weighted$J, hj_replicated$J)
})

#### J with design weights and clusters sums the weighted unit scores within clusters ####
local({
  items_complete <- items[, 4:9]
  model <- set_latpl_model(build_constr_mat(6, 2, "simple"))
  set.seed(2)
  weights <- sample(1:3, nrow(items_complete), replace = TRUE)
  clusters <- rep(1:30, length.out = nrow(items_complete))
  fit_unweighted <- fit_latpl(set_latpl_data(items_complete), model)
  fit_clustered <- fit_latpl(
    set_latpl_data(items_complete, WEIGHTS = weights, CLUSTER = clusters),
    model
  )
  theta <- fit_unweighted$theta + 0.05
  hj <- lavaan.pl:::cpp_HJ(
    theta,
    fit_clustered$engine$model,
    fit_clustered$engine$data
  )
  weighted_scores <- unit_scores(fit_unweighted, theta) *
    fit_clustered$data$weights
  cluster_totals <- rowsum(weighted_scores, clusters)
  cluster_totals <- sweep(cluster_totals, 2, colMeans(cluster_totals))
  expect_equal(hj$J, 30 / 29 * crossprod(cluster_totals) / nrow(items_complete))
  expect_equal(hj$gradient, colSums(weighted_scores))
})

#### Gradient and blocks of the Hessian of the saturated model ####
local({
  data <- set_latpl_data(items_missing[, 4:9], MISSING = "available.cases")
  engine_data <- list(
    N = data$n,
    C_VEC = data$cat,
    FREQ = data$freq,
    UFREQ = data$ufreq,
    Y = data$Y,
    UW = data$uweights
  )
  nthr <- sum(data$cat - 1)
  thresholds_and_correlations <- c(
    unlist(lapply(data$cat - 1, function(m) seq(-0.8, 0.8, length.out = m))),
    seq(.1, .3, length.out = 15)
  )
  saturated <- function(x, order) {
    lavaan.pl:::cpp_saturated(
      x[1:nthr],
      x[-(1:nthr)],
      engine_data,
      order,
      FALSE
    )
  }
  fit <- saturated(thresholds_and_correlations, 2L)
  hessian_numeric <- -numDeriv::jacobian(
    function(x) saturated(x, 1L)$gradient,
    thresholds_and_correlations
  )
  correlations_block <- hessian_numeric[-(1:nthr), -(1:nthr)]
  expect_equal(
    fit$gradient,
    numDeriv::grad(
      function(x) saturated(x, 0L)$ll,
      thresholds_and_correlations
    ),
    tolerance = 1e-6
  )
  expect_equal(
    cbind(fit$htt, fit$htr),
    hessian_numeric[1:nthr, ],
    tolerance = 1e-6
  )
  expect_equal(fit$hrr, diag(correlations_block), tolerance = 1e-6)
  expect_true(
    max(abs(correlations_block - diag(diag(correlations_block)))) < 1e-6
  )
  expect_equal(fit$hschur, fit$htt - fit$htr %*% (t(fit$htr) / fit$hrr))
})
