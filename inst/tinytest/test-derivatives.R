#### Data and models ####
loadings_structure <- lavaan.pl:::build_constr_mat(9, 3, "simple")
loadings_structure[4, 1] <- NA
latent_correlations <- matrix(.3, 3, 3)
diag(latent_correlations) <- 1
items <- lavaan.pl:::sim_data(
  300,
  lavaan.pl:::gen_loadings(loadings_structure, SEED = 1, LB = .4, UB = .8),
  c(-1, 0.3, 1),
  latent_correlations,
  SEED = 3
)
items[, 1:3] <- 1 * (items[, 1:3] >= 2)
items_missing <- items
set.seed(4)
items_missing[sample(length(items_missing), 150)] <- NA
data_missing <- lavaan.pl:::set_latpl_data(items_missing, MISSING = "available.cases")
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
        x$W[] <- 0
        x
      })
      group <- data$group[i] + 1
      if (any(data$numeric)) {
        unit_data[[group]]$W[sum(data$group[seq_len(i)] == data$group[i])] <- 1
      }
      stratum <- data$stratum[i]
      freq <- lavaan.pl:::pairs_freq(response, data$cat)
      ufreq <- unlist(lapply(seq_len(data$p), function(j) {
        data$uweights[i] *
          (seq_len(data$cat[j]) - 1 == response[j] & !is.na(response[j]))
      }))
      unit_data[[group]]$FREQ[stratum * length(freq) + seq_along(freq)] <- freq
      if (length(unit_data[[group]]$UFREQ)) {
        unit_data[[group]]$UFREQ[stratum * length(ufreq) + seq_along(ufreq)] <- ufreq
      }
      lavaan.pl:::cpp_pl(theta, fit$engine$model, unit_data, 1L)$gradient
    },
    numeric(length(theta))
  ))
}

pair_gradients <- function(fit, theta) {
  data <- fit$data
  first_item <- unlist(lapply(2:data$p, function(k) rep(k, k - 1)))
  second_item <- unlist(lapply(2:data$p, function(k) seq_len(k - 1)))
  pair_offsets <- cumsum(c(0, data$cat[first_item] * data$cat[second_item]))
  item_offsets <- cumsum(c(0, data$cat))
  stacked <- function(cells, size, strata) {
    unlist(lapply(seq_len(strata) - 1, function(t) t * size + cells))
  }
  Reduce(
    cbind,
    lapply(seq_along(fit$engine$data), function(group) {
      group_data <- fit$engine$data[[group]]
      strata <- ncol(group_data$X)
      gradient <- function(freq, ufreq) {
        engine_data <- lapply(fit$engine$data, function(x) {
          x$FREQ[] <- 0
          x$UFREQ[] <- 0
          x
        })
        engine_data[[group]]$FREQ <- freq
        engine_data[[group]]$UFREQ <- ufreq
        lavaan.pl:::cpp_pl(theta, fit$engine$model, engine_data, 1L)$gradient
      }
      univariate <- sapply(seq_len(data$p), function(j) {
        cells <- stacked(item_offsets[j] + seq_len(data$cat[j]), sum(data$cat), strata)
        gradient(
          0 * group_data$FREQ,
          replace(0 * group_data$UFREQ, cells, group_data$UFREQ[cells])
        )
      })
      sapply(seq_along(first_item), function(r) {
        cells <- stacked((pair_offsets[r] + 1):pair_offsets[r + 1], max(pair_offsets), strata)
        gradient(
          replace(0 * group_data$FREQ, cells, group_data$FREQ[cells]),
          0 * group_data$UFREQ
        ) +
          (univariate[, first_item[r]] + univariate[, second_item[r]]) /
            (data$p - 1)
      })
    })
  )
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
  theta[block == "gamma"] <- seq(.4, -.3, length.out = size("gamma"))
  theta[block == "kappa"] <- seq(-.2, .3, length.out = size("kappa"))

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
  if (!any(fit$data$numeric)) {
    per_pair <- pair_gradients(fit, theta)
    expect_equal(rowSums(per_pair), gradient)
    expect_equal(
      lavaan.pl:::cpp_HJ(theta, fit$engine$model, fit$engine$data, TRUE)$pairvar,
      tcrossprod(per_pair)
    )
  }
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
  suppressWarnings(lavaan.pl:::fit_latpl(
    data,
    model,
    CONTROL = set_qn_control(MAXEVAL = 1)
  ))
}

#### Gradient, H, J and Jacobians of the latent covariance forms and structural models ####
for (latent_form in c("orthogonal", "lkj", "direct")) {
  for (stdlv in c(TRUE, FALSE)) {
    for (structural in c(FALSE, TRUE)) {
      expect_derivatives(engine_only(
        data_missing,
        lavaan.pl:::set_latpl_model(
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
          lavaan.pl:::set_latpl_model(
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
  data_two_groups <- lavaan.pl:::set_latpl_data(
    items_missing,
    MISSING = "available.cases",
    GROUP = rep(c("a", "b"), length.out = nrow(items_missing))
  )
  model_two_groups <- lavaan.pl:::set_latpl_groups(list(
    lavaan.pl:::set_latpl_model(loadings_structure, LATCOV = "direct"),
    lavaan.pl:::set_latpl_model(
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

#### Gradient, H, J and Jacobians with covariates ####
if (at_home()) local({
  set.seed(5)
  covariates <- cbind(
    rbinom(nrow(items_missing), 1, .5),
    sample(c(-1, 0, 2), nrow(items_missing), TRUE)
  )
  data_covariates <- lavaan.pl:::set_latpl_data(
    items_missing,
    MISSING = "available.cases",
    X = covariates
  )
  gamma <- matrix(c(NA, NA, 0.2, NA, 0, NA), 3, 2)
  kappa <- matrix(0, 9, 2)
  kappa[c(1, 4), 1] <- NA
  kappa[6, 2] <- 0.1
  kappa[8, 2] <- NA
  for (latent_form in c("lkj", "direct")) {
    for (stdlv in c(TRUE, FALSE)) {
      for (structural in c(FALSE, TRUE)) {
        expect_derivatives(engine_only(
          data_covariates,
          lavaan.pl:::set_latpl_model(
            loadings_structure,
            STDLV = stdlv,
            LATCOV = latent_form,
            REGMAT = if (structural) regressions,
            RESCOV = if (structural) residual_covariances,
            GAMMA = gamma,
            KAPPA = kappa
          )
        ))
      }
    }
  }
  expect_derivatives(engine_only(
    data_covariates,
    lavaan.pl:::set_latpl_model(
      loadings_structure,
      STDLV = FALSE,
      PARAMETERIZATION = "theta",
      RESCOV = diag(c(NA, 1.2, NA, 1, NA, 0.8, NA, 1, NA)),
      INTERCEPTS = c(NA, 0.1, NA, NA, 0, NA, -0.1, NA, NA),
      MEANS = c(NA, 0.2, NA),
      GAMMA = gamma,
      KAPPA = kappa
    )
  ))
  expect_derivatives(engine_only(
    data_covariates,
    lavaan.pl:::set_latpl_model(
      loadings_structure,
      STDLV = FALSE,
      REGMAT = regressions,
      INTERCEPTS = c(NA, 0.1, NA, NA, 0, NA, -0.1, NA, NA),
      MEANS = c(NA, 0.2, NA),
      SCALES = c(NA, 1.1, NA, NA, 1, NA, 0.9, NA, NA),
      GAMMA = gamma,
      KAPPA = kappa
    )
  ))
  data_covariate_groups <- lavaan.pl:::set_latpl_data(
    items_missing,
    MISSING = "available.cases",
    GROUP = rep(c("a", "b"), length.out = nrow(items_missing)),
    X = covariates
  )
  expect_derivatives(engine_only(
    data_covariate_groups,
    lavaan.pl:::set_latpl_groups(list(
      lavaan.pl:::set_latpl_model(loadings_structure, GAMMA = gamma),
      lavaan.pl:::set_latpl_model(loadings_structure, KAPPA = kappa)
    ))
  ))
})

#### Gradient, H, J and Jacobians with numeric items ####
if (at_home()) local({
  numeric_items <- c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)
  set.seed(6)
  mixed <- items
  mixed[, numeric_items] <- 0.5 * mixed[, numeric_items] + rnorm(3 * nrow(items))
  mixed[sample(length(mixed), 100)] <- NA
  covariates <- cbind(rnorm(nrow(items)), rbinom(nrow(items), 1, .5))
  data_mixed <- lavaan.pl:::set_latpl_data(mixed, MISSING = "pairwise", NUMERIC = numeric_items)
  data_covariates <- lavaan.pl:::set_latpl_data(
    mixed,
    MISSING = "pairwise",
    NUMERIC = numeric_items,
    X = covariates
  )
  data_groups <- lavaan.pl:::set_latpl_data(
    mixed,
    MISSING = "pairwise",
    NUMERIC = numeric_items,
    X = covariates,
    GROUP = rep(c("a", "b"), length.out = nrow(items))
  )
  gamma <- matrix(c(NA, NA, 0.2, NA, 0, NA), 3, 2)
  kappa <- matrix(0, 9, 2)
  kappa[c(1, 2), 1] <- NA
  kappa[5, 2] <- NA
  for (latent_form in c("lkj", "direct")) {
    for (stdlv in c(TRUE, FALSE)) {
      for (structural in c(FALSE, TRUE)) {
        model_arguments <- list(
          loadings_structure,
          STDLV = stdlv,
          LATCOV = latent_form,
          REGMAT = if (structural) regressions,
          RESCOV = if (structural) residual_covariances,
          NUMERIC = numeric_items
        )
        expect_derivatives(engine_only(
          data_mixed,
          do.call(lavaan.pl:::set_latpl_model, model_arguments)
        ))
        expect_derivatives(engine_only(
          data_covariates,
          do.call(lavaan.pl:::set_latpl_model, c(model_arguments, list(GAMMA = gamma, KAPPA = kappa)))
        ))
      }
    }
  }
  rescov_theta <- diag(c(NA, 1.2, NA, 1, NA, 0.8, NA, 1, NA))
  rescov_theta[1, 2] <- rescov_theta[2, 1] <- NA
  expect_derivatives(engine_only(
    data_covariates,
    lavaan.pl:::set_latpl_model(
      loadings_structure,
      STDLV = FALSE,
      PARAMETERIZATION = "theta",
      RESCOV = rescov_theta,
      INTERCEPTS = c(NA, 0.1, NA, NA, NA, NA, -0.1, NA, NA),
      MEANS = c(NA, 0.2, NA),
      GAMMA = gamma,
      KAPPA = kappa,
      NUMERIC = numeric_items
    )
  ))
  expect_derivatives(engine_only(
    data_covariates,
    lavaan.pl:::set_latpl_model(
      loadings_structure,
      STDLV = FALSE,
      REGMAT = regressions,
      INTERCEPTS = c(NA, 0.1, NA, NA, 0.2, NA, -0.1, NA, NA),
      MEANS = c(NA, 0.2, NA),
      SCALES = c(NA, 1.1, NA, NA, 1.3, NA, 0.9, NA, NA),
      GAMMA = gamma,
      KAPPA = kappa,
      NUMERIC = numeric_items
    )
  ))
  expect_derivatives(engine_only(
    data_groups,
    lavaan.pl:::set_latpl_groups(list(
      lavaan.pl:::set_latpl_model(loadings_structure, GAMMA = gamma, NUMERIC = numeric_items),
      lavaan.pl:::set_latpl_model(loadings_structure, KAPPA = kappa, NUMERIC = numeric_items)
    ))
  ))
})

#### Gradient, blocks of the Hessian and scores of the saturated model with covariates and numeric items ####
local({
  set.seed(8)
  numeric_items <- c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE)
  mixed <- items[, 4:9]
  mixed[, numeric_items] <- mixed[, numeric_items] + rnorm(2 * nrow(mixed))
  mixed[sample(length(mixed), 40)] <- NA
  covariates <- cbind(rnorm(nrow(mixed)), rbinom(nrow(mixed), 1, .5))
  cases <- list(
    lavaan.pl:::set_latpl_data(items[, 4:9], X = covariates[, 2, drop = FALSE]),
    lavaan.pl:::set_latpl_data(items_missing[, 4:9], MISSING = "available.cases", X = covariates[, 2, drop = FALSE]),
    lavaan.pl:::set_latpl_data(items[, 4:9], X = covariates),
    lavaan.pl:::set_latpl_data(mixed, MISSING = "pairwise", NUMERIC = numeric_items),
    lavaan.pl:::set_latpl_data(mixed, MISSING = "pairwise", NUMERIC = numeric_items, X = covariates)
  )
  for (data in cases) {
    numeric_item <- data$numeric
    nx <- ncol(data$X)
    engine_data <- lavaan.pl:::engine(
      data,
      lavaan.pl:::set_latpl_model(
        lavaan.pl:::build_constr_mat(6, 2, "simple"),
        GAMMA = if (nx) matrix(0, 2, nx),
        NUMERIC = numeric_item
      )
    )$data[[1]]
    items_block <- c(
      unlist(lapply((data$cat - 1)[!numeric_item], function(m) seq(-0.8, 0.8, length.out = m))),
      seq(-.1, .2, length.out = 6 * nx),
      rep(-.3, sum(numeric_item)),
      rep(.9, sum(numeric_item))
    )
    ni <- length(items_block)
    point <- c(items_block, seq(.1, .3, length.out = 15))
    saturated <- function(x, order, scores = FALSE) {
      lavaan.pl:::cpp_saturated(x[1:ni], x[-(1:ni)], engine_data, order, scores)
    }
    at_point <- saturated(point, 2L)
    expect_equal(
      at_point$gradient,
      numDeriv::grad(function(x) saturated(x, 0L)$ll, point),
      tolerance = 1e-6
    )
    hessian_numeric <- -numDeriv::jacobian(function(x) saturated(x, 1L)$gradient, point)
    hessian_numeric <- (hessian_numeric + t(hessian_numeric)) / 2
    expect_equal(at_point$htt, hessian_numeric[1:ni, 1:ni], tolerance = 1e-5)
    expect_equal(at_point$htr, hessian_numeric[1:ni, -(1:ni)], tolerance = 1e-5)
    expect_equal(at_point$hrr, diag(hessian_numeric[-(1:ni), -(1:ni)]), tolerance = 1e-5)
    with_scores <- saturated(point, 1L, TRUE)
    weights_rows <- engine_data$W
    expect_equal(
      as.vector(with_scores$st %*% weights_rows),
      at_point$gradient[1:ni],
      tolerance = 1e-8
    )
    expect_equal(
      as.vector(with_scores$sr %*% weights_rows),
      at_point$gradient[-(1:ni)],
      tolerance = 1e-8
    )
  }
})

#### Frequency weights give the H and J of replicated rows ####
local({
  items_complete <- items[, 4:9]
  model <- lavaan.pl:::set_latpl_model(lavaan.pl:::build_constr_mat(6, 2, "simple"))
  set.seed(2)
  weights <- sample(1:3, nrow(items_complete), replace = TRUE)
  fit_weighted <- lavaan.pl:::fit_latpl(
    lavaan.pl:::set_latpl_data(
      items_complete,
      WEIGHTS = weights,
      WEIGHTS_TYPE = "frequency"
    ),
    model,
    CONTROL = set_qn_control()
  )
  fit_replicated <- lavaan.pl:::fit_latpl(
    lavaan.pl:::set_latpl_data(items_complete[
      rep(seq_len(nrow(items_complete)), weights),
    ]),
    model,
    CONTROL = set_qn_control()
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
  model <- lavaan.pl:::set_latpl_model(lavaan.pl:::build_constr_mat(6, 2, "simple"))
  set.seed(2)
  weights <- sample(1:3, nrow(items_complete), replace = TRUE)
  clusters <- rep(1:30, length.out = nrow(items_complete))
  fit_unweighted <- lavaan.pl:::fit_latpl(lavaan.pl:::set_latpl_data(items_complete), model, CONTROL = set_qn_control())
  fit_clustered <- lavaan.pl:::fit_latpl(
    lavaan.pl:::set_latpl_data(items_complete, WEIGHTS = weights, CLUSTER = clusters),
    model,
    CONTROL = set_qn_control()
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

#### J with clusters in two groups sums the unit scores within the clusters of each group ####
local({
  items_complete <- items[, 4:9]
  model <- lavaan.pl:::set_latpl_model(lavaan.pl:::build_constr_mat(6, 2, "simple"))
  model_two_groups <- lavaan.pl:::set_latpl_groups(list(model, model))
  group <- rep(1:2, each = 150)
  clusters <- rep(1:15, length.out = 300)
  fit_unclustered <- lavaan.pl:::fit_latpl(
    lavaan.pl:::set_latpl_data(items_complete, GROUP = group),
    model_two_groups,
    CONTROL = set_qn_control()
  )
  fit_clustered <- lavaan.pl:::fit_latpl(
    lavaan.pl:::set_latpl_data(items_complete, CLUSTER = clusters, GROUP = group),
    model_two_groups,
    CONTROL = set_qn_control()
  )
  theta <- fit_unclustered$theta + 0.05
  hj <- lavaan.pl:::cpp_HJ(
    theta,
    fit_clustered$engine$model,
    fit_clustered$engine$data
  )
  scores <- unit_scores(fit_unclustered, theta)
  within_group <- lapply(1:2, function(g) {
    cluster_totals <- rowsum(scores[group == g, ], clusters[group == g])
    cluster_totals <- sweep(cluster_totals, 2, colMeans(cluster_totals))
    15 / 14 * crossprod(cluster_totals)
  })
  expect_equal(hj$J, (within_group[[1]] + within_group[[2]]) / 300)
})

#### Gradient and blocks of the Hessian of the saturated model ####
local({
  data <- lavaan.pl:::set_latpl_data(items_missing[, 4:9], MISSING = "available.cases")
  engine_data <- list(
    N = data$n,
    C_VEC = data$cat,
    TABLES = TRUE,
    FREQ = data$freq,
    UFREQ = data$ufreq,
    X = matrix(0, 0, 1),
    S = numeric(data$n),
    W = rep(1, data$n),
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
