#' Fit factor models for ordinal and numeric data by pairwise likelihood
#'
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_model()], or from [lavaanpl_groups()] for several groups.
#' @param CONTROL Control object from [set_control()], which selects the method.
#' @param START Starting parameter vector, in the order of [parameter_names()]. If `NULL`, default values.
#' @param NCORES Number of threads.
#' @param VERBOSE If `TRUE`, print progress.
#'
#' @return An object of class `lavaanpl_internal`: a list with the inputs `data`, `model` and `control`, the
#'   dimensions `dims`, the `method`, the constraints `lincon`, the parameter blocks `blocks`, the [cpp_inputs()]
#'   `cpp_inputs`, the starting values `init`, the estimates `parameter_vector`, the optimiser output `stoFit` (`SA`,
#'   and the burn-in of `SA-QN`) and `numFit` (`QN`, `SA-QN`), the variance `vars` (`SA`), and the timing `RTime`.
#' @seealso [vector_to_matrices()], [compute_var()]
#' @keywords internal
fit_lavaanpl <- function(
  DATA,
  MODEL,
  CONTROL = set_control(),
  START = NULL,
  NCORES = 1,
  VERBOSE = FALSE
) {
  start_time <- Sys.time()
  stopifnot(
    inherits(DATA, "lavaanpl_data"),
    inherits(MODEL, c("lavaanpl_model", "lavaanpl_groups")),
    inherits(CONTROL, c("lavaanpl_qn_control", "lavaanpl_sa_control")),
    DATA$p == MODEL$p
  )
  if (inherits(MODEL, "lavaanpl_model")) {
    MODEL <- lavaanpl_groups(list(MODEL))
  }
  models <- MODEL$MODELS
  if (
    !all(vapply(
      models,
      function(x) identical(x$NUMERIC, unname(DATA$numeric)),
      logical(1)
    ))
  ) {
    stop("DATA and MODEL must declare the same numeric items.")
  }
  if (DATA$ngroups != length(models)) {
    stop("DATA has ", DATA$ngroups, " groups, MODEL has ", length(models), ".")
  }
  direct <- vapply(models, function(x) x$LATCOV == "direct", logical(1))
  nthr <- sum(DATA$cat - 1)
  block <- parameter_blocks(MODEL, nthr)
  METHOD <- control_method(CONTROL)
  if (METHOD == "SA" && CONTROL$average_epochs == 0) {
    stop(
      "The SA method needs average_epochs of at least 1: use set_control(\"SA-QN\") to skip the averaged phase."
    )
  }
  dims <- model_dims(DATA, MODEL)
  lc <- lincon(MODEL, nthr)
  dims$npar <- length(lc$free_idx)
  RcppParallel::setThreadOptions(numThreads = NCORES)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  inputs <- cpp_inputs(DATA, MODEL)
  fit <- structure(
    list(
      data = DATA,
      model = MODEL,
      control = CONTROL,
      dims = dims,
      method = METHOD,
      lincon = lc,
      blocks = block,
      cpp_inputs = inputs,
      version = as.character(utils::packageVersion("lavaan.pl")),
      cores = NCORES
    ),
    class = "lavaanpl_internal"
  )
  full <- function(x) as.vector(lc$K %*% x + lc$c)
  fit$init <- start_values(START, DATA, MODEL, VERBOSE)
  if (is.null(START)) {
    x <- fit$init[lc$free_idx]
    free_block <- block[lc$free_idx]
    for (i in 1:20) {
      fit$init <- full(x)
      if (implied_variance(MODEL, fit$init) < .95) {
        break
      }
      x[free_block == "lambda"] <- .8 * x[free_block == "lambda"]
      x[free_block == "logsd"] <- x[free_block == "logsd"] + log(.8)
    }
  }
  std <- standardise_covariates(DATA, MODEL, lc, nthr)
  inputs_std <- if (is.null(std$model)) inputs else cpp_inputs(std$data, std$model)
  nll <- function(vector) -cpp_pl(vector, inputs$model, inputs$data)$pl / DATA$n

  # With the direct parametrisation, points with a non positive definite latent
  # covariance matrix are outside the parameter space
  pd <- function(vector) {
    if (!any(direct)) {
      return(TRUE)
    }
    par <- vector_to_matrices(MODEL, vector)
    all(
      !direct |
        vapply(
          par,
          function(x) {
            !inherits(
              try(chol(x$psi), silent = TRUE),
              "try-error"
            )
          },
          logical(1)
        )
    )
  }

  if (METHOD == "QN" && any(direct) && pd(fit$init)) {
    uncov <- function(C, B) {
      keep <- rowSums(C$A[, B == "corr", drop = FALSE] != 0) == 0
      if (any(keep)) list(A = C$A[keep, , drop = FALSE], b = C$b[keep])
    }
    warm_models <- lapply(models, function(x) {
      if (!is.null(x$CONSTRAINTS)) {
        x$CONSTRAINTS <- uncov(x$CONSTRAINTS, parameter_blocks(x, nthr))
      }
      x$LATCOV <- "lkj"
      x
    })
    warm_model <- lavaanpl_groups(
      warm_models,
      if (!is.null(MODEL$CONSTRAINTS)) uncov(MODEL$CONSTRAINTS, block)
    )
    warm <- fit_lavaanpl(
      DATA,
      warm_model,
      CONTROL,
      matrices_to_vector(warm_model, vector_to_matrices(MODEL, fit$init)),
      NCORES
    )
    RcppParallel::setThreadOptions(numThreads = NCORES)
    vector <- matrices_to_vector(
      MODEL,
      vector_to_matrices(warm$model, warm$parameter_vector)
    )
    x <- as.vector(solve(crossprod(lc$K), crossprod(lc$K, vector - lc$c)))
    cov_free <- block[lc$free_idx] == "corr"
    for (i in 1:20) {
      vector <- full(x)
      if (all(is.finite(vector)) && pd(vector)) {
        break
      }
      x[cov_free] <- 0.8 * x[cov_free]
    }
    start <- full(fit$init[lc$free_idx])
    if (
      all(is.finite(vector)) &&
        pd(vector) &&
        isTRUE(nll(vector) < if (pd(start)) nll(start) else Inf)
    ) {
      fit$init <- vector
    }
  }
  fit$init <- full(fit$init[lc$free_idx])

  if (METHOD != "SA") {
    x0 <- std$to(fit$init)
    qn_control <- list(
      grtol = CONTROL$grtol,
      xtol = CONTROL$xtol,
      stepmax = CONTROL$stepmax,
      maxeval = CONTROL$maxeval
    )
    fn <- function(x) {
      vector <- full(x)
      val <- if (pd(vector)) {
        -cpp_pl(vector, inputs_std$model, inputs_std$data)$pl / DATA$n
      } else {
        Inf
      }
      val <- if (is.finite(val)) val else Inf
      evals <<- evals + 1L
      if (val < best) {
        best <<- val
        trail[[length(trail) + 1]] <<- list(
          epoch = evals,
          vector = vector,
          nll = val
        )
      }
      val
    }
    gr <- function(x) {
      vector <- full(x)
      g <- if (pd(vector)) {
        -as.vector(crossprod(
          lc$K,
          cpp_pl(vector, inputs_std$model, inputs_std$data, 1L)$gradient
        )) /
          DATA$n
      } else {
        0 * x
      }
      if (all(is.finite(g))) g else 0 * x
    }
    nstarts <- if (METHOD == "SA-QN") CONTROL$nstarts else 1L
    values <- numeric(nstarts)
    for (s in seq_len(nstarts)) {
      start <- x0
      start_control <- qn_control
      if (METHOD == "SA-QN") {
        if (VERBOSE) {
          message("- Burn-in via pairwise stochastic approximation...")
        }
        burn_control <- CONTROL
        burn_control$seed <- CONTROL$seed + s - 1L
        burn <- run_sa(x0, std$data, inputs_std, lc, burn_control, VERBOSE)
        start <- burn$parameter_vector
        if (length(burn$invh) > 0) {
          start_control$invhessian.lt <- burn$invh[lower.tri(
            burn$invh,
            diag = TRUE
          )]
        }
      }
      if (VERBOSE) {
        message("- Optimising with ucminf...")
      }
      evals <- 0L
      best <- Inf
      trail <- list()
      opt <- ucminf::ucminf(
        start[lc$free_idx],
        fn,
        gr,
        hessian = 0,
        control = start_control
      )
      opt$path <- list(
        epoch = vapply(trail, function(x) x$epoch, integer(1)),
        parameter_vector = lapply(trail, function(x) std$from(x$vector)),
        nll = vapply(trail, function(x) x$nll, numeric(1))
      )
      values[s] <- opt$value
      if (s == 1 || opt$value < fit$numFit$value) {
        fit$numFit <- opt
        if (METHOD == "SA-QN") {
          fit$stoFit <- burn[setdiff(
            names(burn),
            c("parameter_vector", "last_iterate", "burn_in_iterate", "invh")
          )]
          fit$stoFit$path$parameter_vector <- lapply(
            burn$path$parameter_vector,
            std$from
          )
        }
      }
    }
    if (METHOD == "SA-QN") {
      fit$stoFit$start_logl <- -values * DATA$n
    }
    if (DATA$n * diff(range(values)) > 1e-3) {
      warning(
        "The ",
        nstarts,
        " starts ended at pairwise log-likelihoods up to ",
        signif(DATA$n * diff(range(values)), 3),
        " apart, so the pairwise likelihood may have several optima: the highest is returned.",
        call. = FALSE
      )
    }
    fit$parameter_vector <- std$from(full(fit$numFit$par))
    fit$numFit$par <- fit$parameter_vector[lc$free_idx]
    if (
      !fit$numFit$convergence %in% 1:2 &&
        fit$numFit$info[["maxgradient"]] > 1e-3
    ) {
      warning(
        "ucminf stopped with a gradient of ",
        signif(fit$numFit$info[["maxgradient"]], 3),
        " per observation: the estimates may not be at an optimum.",
        call. = FALSE
      )
    }
  } else {
    if (VERBOSE) {
      message("- Optimising via pairwise stochastic approximation...")
    }
    opt <- run_sa(std$to(fit$init), std$data, inputs_std, lc, CONTROL, VERBOSE)
    opt[c("parameter_vector", "last_iterate", "burn_in_iterate")] <- lapply(
      opt[c("parameter_vector", "last_iterate", "burn_in_iterate")],
      std$from
    )
    opt$path$parameter_vector <- lapply(opt$path$parameter_vector, std$from)
    if (opt$skipped_averaging > 0) {
      warning(
        opt$skipped_averaging,
        " SA iterations of the averaged phase were skipped because no feasible step was found.",
        call. = FALSE
      )
    }
    opt$converged <- opt$skipped_averaging == 0 &&
      all(is.finite(opt$parameter_vector))
    fit$stoFit <- opt[setdiff(names(opt), c("parameter_vector", "invh"))]
    fit$parameter_vector <- opt$parameter_vector
    fit <- check_sa(fit)
  }
  fit$RTime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (VERBOSE) {
    message("Done! (", round(fit$RTime, 2), " secs)")
  }
  fit
}


#' End-of-run check of the stochastic approximation
#'
#' @param FIT Object of class `lavaanpl_internal`.
#' @return `FIT`, with the variance `vars` and, in `stoFit`, the `check` and the updated `converged`.
#' @seealso [fit_lavaanpl()], [compute_var()]
#' @keywords internal
check_sa <- function(FIT) {
  FIT$stoFit$check <- list(done = FALSE, pd = NA, distance = NA, noise = NA)
  if (FIT$dims$npar > 2000) {
    if (!FIT$stoFit$plateau) {
      FIT$stoFit$converged <- FALSE
      warning(
        "The Adam phase did not reach a plateau of the pairwise likelihood: increase adam_epochs.",
        call. = FALSE
      )
    }
    return(FIT)
  }
  vars <- tryCatch(compute_var(FIT), lavaanpl_not_pd = function(e) NULL)
  FIT$stoFit$check$done <- TRUE
  FIT$stoFit$check$pd <- !is.null(vars)
  if (is.null(vars)) {
    FIT$stoFit$converged <- FALSE
    warning(
      "The observed information is not positive definite at the SA estimate, which is not an optimum of the pairwise likelihood: use set_control(\"SA-QN\") or other starting values.",
      call. = FALSE
    )
    return(FIT)
  }
  se <- sqrt(vars$asymptotic_variance / FIT$data$n)
  distance <- FIT$stoFit$check$distance <- max(vars$distance)
  FIT$stoFit$check$noise <- max((sqrt(pmax(vars$optimisation_noise, 0)) / se)[
    se > 0
  ])
  if (distance > max(1, 4 * FIT$stoFit$check$noise)) {
    FIT$stoFit$converged <- FALSE
    warning(
      "The SA estimate is ",
      signif(distance, 3),
      " standard errors from the optimum of the pairwise likelihood on its worst parameter, more than the optimisation noise explains: use other starting values or set_control(\"SA-QN\").",
      call. = FALSE
    )
  } else if (distance > 1) {
    warning(
      "The SA estimate is ",
      signif(distance, 3),
      " standard errors from the optimum of the pairwise likelihood on its worst parameter, because of the optimisation noise, which the standard errors account for. For more precision, increase average_epochs or use set_control(\"SA-QN\"); see diagnostics().",
      call. = FALSE
    )
  }
  FIT$vars <- vars
  FIT
}

#' Standardised covariates for the optimisers
#'
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param LINCON Output of [lincon()].
#' @param NTHR Number of thresholds in each group.
#' @return A list with the data `data` and model `model` for the optimisers, and the maps `to` and `from` their
#'   parameter vector. Without covariates, or with constraints on the thresholds, intercepts or effects of covariates,
#'   `DATA`, no `model` and identities.
#' @seealso [fit_lavaanpl()]
#' @keywords internal
standardise_covariates <- function(DATA, MODEL, LINCON, NTHR) {
  raw <- list(data = DATA, model = NULL, to = identity, from = identity)
  if (ncol(DATA$X) == 0) {
    return(raw)
  }
  models <- MODEL$MODELS
  if (!any(vapply(models, function(x) ncol(x$GAMMA) > 0, logical(1)))) {
    return(raw)
  }
  centre <- colMeans(DATA$X)
  spread <- apply(DATA$X, 2, stats::sd)
  spread[spread == 0] <- 1
  block <- parameter_blocks(MODEL, NTHR)
  touched <- block %in%
    c("tau", "gamma", "kappa") |
    (block == "nu" & any(DATA$numeric))
  fixed_nu <- any(vapply(
    models,
    function(x) any(!is.na(x$NU[x$NUMERIC])),
    logical(1)
  ))
  if (fixed_nu || !all(which(touched) %in% LINCON$free_idx)) {
    if (any(abs(centre) > 3 * spread)) {
      message(
        "Constraints on thresholds, intercepts or effects of covariates prevent standardising the covariates: consider centring them."
      )
    }
    return(raw)
  }
  scale_x <- function(M) sweep(M, 2, spread, "*")
  rescale <- function(M) {
    M$GAMMA <- scale_x(M$GAMMA)
    M$KAPPA <- scale_x(M$KAPPA)
    M
  }
  std_models <- lapply(models, rescale)
  std_model <- replace(MODEL, "MODELS", list(std_models))
  standard <- function(X) sweep(sweep(X, 2, centre, "-"), 2, spread, "/")
  std_data <- DATA
  std_data$X <- standard(DATA$X)
  std_data$strata <- lapply(DATA$strata, standard)
  if (!is.null(DATA$val)) {
    std_data$val$strata <- lapply(DATA$val$strata, standard)
  }
  map <- function(vector, from_model, to_model, to_standard) {
    pars <- vector_to_matrices(from_model, vector)
    pars <- Map(
      function(par, m) {
        if (is.null(par$gamma)) {
          return(par)
        }
        gamma_raw <- if (to_standard) {
          par$gamma
        } else {
          sweep(par$gamma, 2, spread, "/")
        }
        kappa_raw <- if (to_standard) {
          par$kappa
        } else {
          sweep(par$kappa, 2, spread, "/")
        }
        loadings <- par$lambda %*% solve(diag(ncol(par$lambda)) - par$beta)
        offset <- as.vector((loadings %*% gamma_raw + kappa_raw) %*% centre) *
          (if (to_standard) 1 else -1)
        par$tau <- par$tau - rep(offset, DATA$cat - 1)
        par$nu[m$NUMERIC] <- par$nu[m$NUMERIC] + offset[m$NUMERIC]
        par$gamma <- if (to_standard) scale_x(par$gamma) else gamma_raw
        par$kappa <- if (to_standard) scale_x(par$kappa) else kappa_raw
        par
      },
      pars,
      models
    )
    matrices_to_vector(to_model, pars)
  }
  list(
    data = std_data,
    model = std_model,
    to = function(vector) map(vector, MODEL, std_model, TRUE),
    from = function(vector) map(vector, std_model, MODEL, FALSE)
  )
}

#' Inputs of the C++ routines
#'
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param TABLES Use the frequency tables (`TRUE`) or the individual observations (`FALSE`). If `NA`, as set by
#'   [lavaanpl_data()].
#' @return A list with `model` and `data`, each with one element per group: the model matrices, and the frequency
#'   tables, covariate strata, distinct response patterns and weights.
#' @seealso [fit_lavaanpl()], [run_sa()]
#' @keywords internal
cpp_inputs <- function(DATA, MODEL, TABLES = NA) {
  models <- MODEL$MODELS
  gw <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
  jw <- if (length(DATA$jweights)) DATA$jweights else rep(1, DATA$n)
  ns <- vapply(DATA$strata, nrow, numeric(1))
  uc <- if (length(DATA$ufreq)) sum(DATA$cat) else 0
  ou <- cumsum(c(0, ns[-length(ns)])) * uc
  ot <- cumsum(c(0, (ns * DATA$tables)[-length(ns)]))
  nc <- if (any(DATA$tables)) length(DATA$freq) / sum(ns * DATA$tables) else 0
  nx <- ncol(DATA$X)
  list(
    model = lapply(models, function(x) {
      list(
        LAMBDA = x$LAMBDA,
        LOGSD = log(sqrt(diag(x$PSI))),
        LATCOV = x$LATCOV,
        BETA = if (all(x$BETA %in% 0)) matrix(0, 0, 0) else x$BETA,
        THETA = if (all(x$THETA[upper.tri(x$THETA)] %in% 0)) {
          matrix(0, 0, 0)
        } else {
          x$THETA
        },
        NU = x$NU,
        ALPHA = x$ALPHA,
        LOGSCALE = log(item_scales(x)),
        PARAMETERIZATION = x$PARAMETERIZATION,
        GAMMA = if (ncol(x$GAMMA) == nx) x$GAMMA else matrix(0, x$q, nx),
        KAPPA = if (ncol(x$KAPPA) == nx) x$KAPPA else matrix(0, x$p, nx)
      )
    }),
    data = lapply(seq_len(DATA$ngroups) - 1, function(g) {
      rows <- which(DATA$group == g)
      key <- if (any(DATA$numeric)) {
        seq_along(rows)
      } else {
        do.call(
          paste,
          as.data.frame(cbind(DATA$stratum[rows], DATA$Y[rows, , drop = FALSE]))
        )
      }
      first <- !duplicated(key)
      id <- match(key, key[first])
      cluster_weights <- matrix(0, 0, 0)
      if (length(DATA$cluster)) {
        clusters <- factor(DATA$cluster[rows])
        C <- nlevels(clusters)
        A <- tapply(
          gw[rows],
          list(
            factor(id, seq_len(sum(first))),
            clusters
          ),
          sum,
          default = 0
        )
        cluster_weights <- (A - rowMeans(A)) * sqrt(C / (C - 1))
      }
      tables <- if (is.na(TABLES)) DATA$tables[g + 1] else TABLES
      freq <- if (tables == DATA$tables[g + 1]) {
        if (tables) {
          DATA$freq[ot[g + 1] * nc + seq_len(ns[g + 1] * nc)]
        } else {
          numeric(0)
        }
      } else if (tables) {
        cpp_pairs_freq(
          DATA$Y[rows, , drop = FALSE],
          DATA$cat,
          if (length(DATA$weights)) DATA$weights[rows] else numeric(0),
          DATA$stratum[rows]
        )
      } else {
        numeric(0)
      }
      list(
        N = length(rows),
        CAT = DATA$cat,
        TABLES = tables,
        FREQ = freq,
        UFREQ = DATA$ufreq[ou[g + 1] + seq_len(ns[g + 1] * uc)],
        STRATA = t(DATA$strata[[g + 1]]),
        STRATUM = DATA$stratum[rows[first]],
        WEIGHTS = as.numeric(rowsum(gw[rows], id)),
        Y = DATA$Y[rows[first], , drop = FALSE],
        UWEIGHTS = DATA$uweights[rows[first]],
        JWEIGHTS = sqrt(as.numeric(rowsum(jw[rows]^2, id))),
        CLUSTER_WEIGHTS = cluster_weights
      )
    })
  )
}

#' Starting values
#'
#' @param START Starting parameter vector, or `NULL`.
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param VERBOSE If `TRUE`, report user-supplied starting values.
#' @return `START` when supplied, otherwise thresholds from the margins of each group, loadings from
#'   [signed_loadings()], means and standard deviations of the numeric items, and fixed values for the rest.
#' @seealso [fit_lavaanpl()], [spec_from_lavaan()]
#' @keywords internal
start_values <- function(START, DATA, MODEL, VERBOSE) {
  models <- MODEL$MODELS
  nthr <- sum(DATA$cat - 1)
  vector <- START
  if (is.null(vector)) {
    w <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
    vector <- unlist(lapply(seq_along(models), function(g) {
      m <- models[[g]]
      rows <- DATA$group == g - 1
      block <- parameter_blocks(m, nthr)
      th <- rep(0, length(block))
      th[block == "tau"] <- unlist(lapply(seq_len(DATA$p), function(j) {
        pj <- vapply(
          seq_len(DATA$cat[j]) - 1,
          function(a) sum(w[rows & DATA$Y[, j] %in% a]),
          numeric(1)
        )
        stats::qnorm(pmin(pmax(cumsum(pj) / sum(pj), 1e-4), 1 - 1e-4))[
          -DATA$cat[j]
        ]
      }))
      th[block == "beta"] <- 0.25
      spread <- rep(1, DATA$p)
      centre <- rep(0, DATA$p)
      for (j in which(DATA$numeric)) {
        seen <- rows & !is.na(DATA$Y[, j])
        centre[j] <- stats::weighted.mean(DATA$Y[seen, j], w[seen])
        spread[j] <- sqrt(stats::weighted.mean(
          (DATA$Y[seen, j] - centre[j])^2,
          w[seen]
        ))
      }
      start <- signed_loadings(DATA$Y[rows, , drop = FALSE], m, spread)
      free_scales <- is.na(item_scales(m))
      th[block == "lambda"] <- start$loadings[is.na(m$LAMBDA)]
      th[block == "logsd"] <- log(start$sd[is.na(diag(m$PSI))])
      th[block == "nu"] <- centre[is.na(m$NU)]
      th[block == "scale"][m$NUMERIC[free_scales]] <- log(spread[
        m$NUMERIC & free_scales
      ])
      if (m$PARAMETERIZATION == "theta") {
        par <- vector_to_matrices(m, th)
        L <- par$lambda %*% solve(diag(m$q) - par$beta)
        th[block == "tau"] <- th[block == "tau"] *
          rep(
            sqrt(
              rowSums((L %*% par$psi) * L) +
                diag(par$theta)
            ),
            DATA$cat - 1
          )
      }
      th
    }))
  } else if (VERBOSE) {
    message("- Initialised at user-supplied values.")
  }
  stopifnot(
    is.numeric(vector),
    length(vector) == length(parameter_blocks(MODEL, nthr)),
    all(is.finite(vector))
  )
  if (!is.null(START) && implied_variance(MODEL, vector) > 1) {
    warning("Negative error variances implied by init parameters")
  }
  vector
}

#' Signed starting loadings
#'
#' @param Y Data matrix of the group.
#' @param MODEL Model object from [lavaanpl_model()].
#' @param SPREAD Standard deviation of the numeric items, one for the others.
#' @return A list with the starting loadings `loadings`, from the first eigenvector of the item correlations of each
#'   factor, and the latent standard deviations `sd`.
#' @seealso [start_values()]
#' @keywords internal
signed_loadings <- function(Y, MODEL, SPREAD) {
  rows <- unique(round(seq(1, nrow(Y), length.out = min(nrow(Y), 2000))))
  L <- MODEL$LAMBDA
  variance <- diag(MODEL$PSI)
  sd <- ifelse(is.na(variance), .25, sqrt(variance))
  L[is.na(L)] <- .25
  for (h in seq_len(MODEL$q)) {
    S <- which(is.na(MODEL$LAMBDA[, h]) | MODEL$LAMBDA[, h] != 0)
    if (length(S) < 2) {
      next
    }
    R <- suppressWarnings(stats::cor(
      Y[rows, S, drop = FALSE],
      use = "pairwise.complete.obs"
    ))
    R[is.na(R)] <- 0
    e <- eigen(R, symmetric = TRUE)
    v <- pmin(pmax(e$vectors[, 1] * sqrt(max(e$values[1] - 1, .1)), -.9), .9)
    v <- if (v[1] < 0) -v else v
    fix <- which(!is.na(MODEL$LAMBDA[S, h]))[1]
    if (is.na(variance[h]) && !is.na(fix)) {
      sd[h] <- abs(v[fix] * SPREAD[S[fix]] / MODEL$LAMBDA[S[fix], h])
      v <- v * sign(v[fix] * MODEL$LAMBDA[S[fix], h]) * SPREAD[S] / sd[h]
    } else {
      if (is.na(variance[h])) {
        sd[h] <- 1
      }
      if (!is.na(fix)) {
        v <- v * sign(v[fix] * MODEL$LAMBDA[S[fix], h])
      }
      v <- v * SPREAD[S] / sd[h]
    }
    free <- is.na(MODEL$LAMBDA[S, h])
    L[S[free], h] <- v[free]
  }
  list(loadings = L, sd = sd)
}

#' Largest variance of the latent responses implied by the factors
#'
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param PARAMETER_VECTOR Parameter vector.
#' @return The largest implied variance of the ordinal items under `PARAMETERIZATION = "delta"`, zero if there are
#'   none. It must be below one.
#' @seealso [fit_lavaanpl()], [start_values()]
#' @keywords internal
implied_variance <- function(MODEL, PARAMETER_VECTOR) {
  max(
    0,
    unlist(Map(
      function(m, par) {
        if (m$PARAMETERIZATION != "delta") {
          return(0)
        }
        L <- par$lambda %*% solve(diag(m$q) - par$beta)
        (rowSums((L %*% par$psi) * L) * par$delta^2)[!m$NUMERIC]
      },
      MODEL$MODELS,
      vector_to_matrices(MODEL, PARAMETER_VECTOR)
    ))
  )
}

#' Stochastic approximation
#'
#' @param START Starting parameter vector.
#' @param DATA Data object from [lavaanpl_data()].
#' @param CPP_INPUTS Output of [cpp_inputs()].
#' @param LINCON Output of [lincon()].
#' @param CONTROL Control object from [set_control()].
#' @param VERBOSE If `TRUE`, print progress.
#' @return A list with the estimate `parameter_vector`, the iterates `last_iterate` and `burn_in_iterate`, the epoch
#'   length and phase counts, the rejected and skipped steps, the trajectory `path` and the timings `time`.
#' @seealso [fit_lavaanpl()], [set_control()]
#' @keywords internal
run_sa <- function(START, DATA, CPP_INPUTS, LINCON, CONTROL, VERBOSE) {
  npar <- length(LINCON$free_idx)
  free <- npar < length(START)
  adam_epochs <- CONTROL$adam_epochs
  if (is.null(adam_epochs)) {
    adam_epochs <- max(
      5L,
      ceiling(
        150 / ceiling(DATA$pairs / min(CONTROL$pairs_per_iteration, DATA$pairs))
      )
    )
  }
  hessian_gain <- switch(
    CONTROL$gain,
    auto = npar <= 2000,
    hessian = TRUE,
    diagonal = FALSE
  )
  out <- cpp_sa(
    START = START[LINCON$free_idx],
    K = if (free) LINCON$K else matrix(0, 0, 0),
    C = LINCON$c,
    MODEL = CPP_INPUTS$model,
    DATA = CPP_INPUTS$data,
    PAIRS_PER_ITERATION = CONTROL$pairs_per_iteration,
    ALPHA = CONTROL$alpha,
    ADAM_EPOCHS = as.integer(adam_epochs),
    PRECONDITION_EPOCHS = CONTROL$precondition_epochs,
    MAX_PRECONDITION_EPOCHS = CONTROL$max_precondition_epochs,
    AVERAGE_EPOCHS = CONTROL$average_epochs,
    TOL = CONTROL$tol,
    DECAY = CONTROL$decay,
    AVG_FROM = CONTROL$avg_from,
    HESSIAN_GAIN = hessian_gain,
    SEED = CONTROL$seed,
    VERBOSE = VERBOSE
  )
  out$path$nll <- c(out$path$nll, rep(NA, sum(out$path$phase == 3)))
  out
}
