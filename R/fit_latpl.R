#'Fit factor models for ordinal and numeric data with pairwise likelihood methods
#'
#'@description
#'
#'Fit latent variable models with one of three methods: the pairwise likelihood estimate by
#'quasi-Newton BFGS started from the burn-in of stochastic approximation and its curvature
#'(`SA-QN`, the default), pairwise stochastic approximation with Polyak-Ruppert averaging, whose
#'standard errors account for the optimisation noise (`SA`), or quasi-Newton BFGS from the
#'starting values alone (`QN`).
#'
#'@param DATA Data object from [set_latpl_data()].
#'@param MODEL Model object from [set_latpl_model()], or from [set_latpl_groups()] for several groups.
#'@param CONTROL Control object, which selects the method: [set_saqn_control()] for
#'  `SA-QN`, [set_sa_control()] for `SA`, [set_qn_control()] for `QN`.
#'@param INIT Starting values object from [set_init()].
#'@param NCORES Integer value setting the number of threads to carry out the
#'  estimation.
#'@param VERBOSE `TRUE` for verbose output
#'
#'@return Object of class `latpl_fit`: a list with the inputs `data`, `model` and
#'  `control`, the starting values `init`, the estimates `theta` (use
#'  [get_par()]), the dimensions `dims`, the `method`, the optimiser output `stoFit` (`SA`,
#'  the burn-in of `SA-QN`, with the pairwise log-likelihood reached from each start in
#'  `start_logl`) or `numFit` (`QN`, `SA-QN`), the variance `vars` at the estimate (`SA`), and
#'  the timing `RTime`.
#' @keywords internal
fit_latpl <- function(
  DATA,
  MODEL,
  CONTROL = set_saqn_control(),
  INIT = set_init(),
  NCORES = 1,
  VERBOSE = FALSE
) {
  start_time <- Sys.time()
  stopifnot(
    inherits(DATA, "latpl_data"),
    inherits(MODEL, c("latpl_model", "latpl_groups")),
    inherits(CONTROL, c("latpl_qn_control", "latpl_sa_control")),
    inherits(INIT, "latpl_init"),
    DATA$p == MODEL$p
  )
  grouped <- inherits(MODEL, "latpl_groups")
  models <- if (grouped) MODEL$MODELS else list(MODEL)
  if (!all(vapply(models, function(x) identical(x$NUMERIC, unname(DATA$numeric)), logical(1)))) {
    stop("DATA and MODEL must declare the same numeric items.")
  }
  if (DATA$ngroups != length(models)) {
    stop("DATA has ", DATA$ngroups, " groups, MODEL has ", length(models), ".")
  }
  direct <- vapply(models, function(x) corrflag_code(x) == 2, logical(1))
  nthr <- sum(DATA$cat - 1)
  METHOD <- control_method(CONTROL)
  if (METHOD == "SA" && CONTROL$AVERAGE_EPOCHS == 0) {
    stop("The SA method needs AVERAGE_EPOCHS of at least 1: use set_saqn_control() to skip the averaged phase.")
  }
  dims <- latpl_dims(DATA, MODEL)
  lc <- lincon(MODEL, nthr)
  dims$npar <- length(lc$free_idx)
  RcppParallel::setThreadOptions(numThreads = NCORES)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  eng <- engine(DATA, MODEL)
  fit <- structure(
    list(
      data = DATA,
      model = MODEL,
      control = CONTROL,
      dims = dims,
      method = METHOD,
      lincon = lc,
      engine = eng,
      version = as.character(utils::packageVersion("lavaan.pl")),
      cores = NCORES
    ),
    class = "latpl_fit"
  )
  full <- function(x) as.vector(lc$K %*% x + lc$c)
  fit$init <- init_theta(INIT, DATA, dims, MODEL, VERBOSE)
  if (is.null(INIT$THETA)) {
    x <- fit$init[lc$free_idx]
    free_block <- theta_block(MODEL, nthr)[lc$free_idx]
    for (i in 1:20) {
      fit$init <- full(x)
      if (implied_variance(MODEL, fit$init) < .95) break
      x[free_block == "lambda"] <- .8 * x[free_block == "lambda"]
      x[free_block == "logsd"] <- x[free_block == "logsd"] + log(.8)
    }
  }
  std <- standardise_covariates(DATA, MODEL, lc, nthr)
  eng_std <- if (is.null(std$model)) eng else engine(std$data, std$model)
  if (METHOD == "QN" && INIT$METHOD == "SA") {
    if (VERBOSE) {
      message("- Initialising with the SA Adam phase...")
    }
    fit$init <- std$from(run_sa(
      std$to(fit$init),
      std$data,
      eng_std,
      lc,
      INIT$CONTROL,
      TRUE,
      FALSE
    )$theta)
  }
  nll <- function(theta) -cpp_pl(theta, eng$model, eng$data)$pl / dims$n

  # With the direct parametrisation, points with a non positive definite latent
  # covariance matrix are outside the parameter space
  pd <- function(theta) {
    if (!any(direct)) {
      return(TRUE)
    }
    par <- if (grouped) {
      theta2par(MODEL, theta)
    } else {
      list(theta2par(MODEL, theta))
    }
    all(
      !direct |
        vapply(
          par,
          function(x) {
            !inherits(
              try(chol(x$latent_correlations), silent = TRUE),
              "try-error"
            )
          },
          logical(1)
        )
    )
  }

  if (METHOD == "QN" && any(direct) && pd(fit$init)) {
    block <- theta_block(MODEL, nthr)
    uncov <- function(C, B) {
      keep <- rowSums(C$A[, B == "corr", drop = FALSE] != 0) == 0
      if (any(keep)) list(A = C$A[keep, , drop = FALSE], b = C$b[keep])
    }
    warm_models <- lapply(models, function(x) {
      if (!is.null(x$CONSTRAINTS)) {
        x$CONSTRAINTS <- uncov(x$CONSTRAINTS, theta_block(x, nthr))
      }
      x$LATCOV <- "lkj"
      x
    })
    warm_model <- if (!grouped) {
      warm_models[[1]]
    } else {
      set_latpl_groups(
        warm_models,
        if (!is.null(MODEL$CONSTRAINTS)) uncov(MODEL$CONSTRAINTS, block)
      )
    }
    warm <- fit_latpl(
      DATA,
      warm_model,
      CONTROL,
      set_init(THETA = par2theta(warm_model, theta2par(MODEL, fit$init))),
      NCORES
    )
    RcppParallel::setThreadOptions(numThreads = NCORES)
    theta <- par2theta(MODEL, get_par(warm, OPTION = "list"))
    x <- as.vector(solve(crossprod(lc$K), crossprod(lc$K, theta - lc$c)))
    cov_free <- block[lc$free_idx] == "corr"
    for (i in 1:20) {
      theta <- full(x)
      if (all(is.finite(theta)) && pd(theta)) {
        break
      }
      x[cov_free] <- 0.8 * x[cov_free]
    }
    start <- full(fit$init[lc$free_idx])
    if (
      all(is.finite(theta)) &&
        pd(theta) &&
        isTRUE(nll(theta) < if (pd(start)) nll(start) else Inf)
    ) {
      fit$init <- theta
    }
  }
  fit$init <- full(fit$init[lc$free_idx])

  if (METHOD %in% c("QN", "SA-QN")) {
    x0 <- std$to(fit$init)
    qn_control <- list(
      grtol = CONTROL$GRTOL,
      xtol = CONTROL$XTOL,
      stepmax = CONTROL$STEPMAX,
      maxeval = CONTROL$MAXEVAL
    )
    fn <- function(x) {
      theta <- full(x)
      val <- if (pd(theta)) -cpp_pl(theta, eng_std$model, eng_std$data)$pl / dims$n else Inf
      if (is.finite(val)) val else Inf
    }
    gr <- function(x) {
      theta <- full(x)
      g <- if (pd(theta)) {
        -as.vector(crossprod(
          lc$K,
          cpp_pl(theta, eng_std$model, eng_std$data, 1L)$gradient
        )) /
          dims$n
      } else {
        0 * x
      }
      if (all(is.finite(g))) g else 0 * x
    }
    nstarts <- if (METHOD == "SA-QN") CONTROL$NSTARTS else 1L
    values <- numeric(nstarts)
    for (s in seq_len(nstarts)) {
      start <- x0
      start_control <- qn_control
      if (METHOD == "SA-QN") {
        if (VERBOSE) {
          message("- Burn-in via pairwise stochastic approximation...")
        }
        burn_control <- CONTROL
        burn_control$SEED <- CONTROL$SEED + s - 1L
        burn <- run_sa(x0, std$data, eng_std, lc, burn_control, FALSE, VERBOSE)
        start <- burn$theta
        if (length(burn$invh) > 0) start_control$invhessian.lt <- burn$invh[lower.tri(burn$invh, diag = TRUE)]
      }
      if (VERBOSE) {
        message("- Optimising with ucminf...")
      }
      opt <- ucminf::ucminf(
        start[lc$free_idx],
        fn,
        gr,
        hessian = 0,
        control = start_control
      )
      values[s] <- opt$value
      if (s == 1 || opt$value < fit$numFit$value) {
        fit$numFit <- opt
        if (METHOD == "SA-QN") fit$stoFit <- burn[setdiff(names(burn), c("theta", "theta_last", "theta_burn", "invh", "path"))]
      }
    }
    if (METHOD == "SA-QN") fit$stoFit$start_logl <- -values * dims$n
    if (dims$n * diff(range(values)) > 1e-3) {
      warning(
        "The ",
        nstarts,
        " starts ended at pairwise log-likelihoods up to ",
        signif(dims$n * diff(range(values)), 3),
        " apart, so the pairwise likelihood may have several optima: the highest is returned.",
        call. = FALSE
      )
    }
    fit$theta <- std$from(full(fit$numFit$par))
    fit$numFit$par <- fit$theta[lc$free_idx]
    if (!fit$numFit$convergence %in% 1:2 && fit$numFit$info[["maxgradient"]] > 1e-3) {
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
    opt <- run_sa(std$to(fit$init), std$data, eng_std, lc, CONTROL, FALSE, VERBOSE)
    opt[c("theta", "theta_last", "theta_burn")] <- lapply(opt[c("theta", "theta_last", "theta_burn")], std$from)
    opt$path$theta <- lapply(opt$path$theta, std$from)
    if (opt$skipped_averaging > 0) {
      warning(
        opt$skipped_averaging,
        " SA iterations of the averaged phase were skipped because no feasible step was found.",
        call. = FALSE
      )
    }
    opt$converged <- opt$skipped_averaging == 0 && all(is.finite(opt$theta))
    fit$stoFit <- opt[setdiff(names(opt), c("theta", "invh"))]
    fit$theta <- opt$theta
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
#' Stores the variance at the estimate, which the check needs and [compute_var()] reuses. The fit has failed if the
#' observed information is not positive definite, or if the estimate is further from the optimum of the pairwise
#' likelihood than one standard error and four times the predicted noise on the worst parameter, because the noise cannot
#' explain it. A larger distance that the noise explains only warns: the standard errors account for it. Above 2000 free
#' parameters there is no check, and the fit has failed if the Adam phase did not reach a plateau.
#'
#' @param FIT Object of class latpl_fit.
#' @return `FIT`, with `vars` and, in `stoFit`, the check `check` and the updated `converged`.
#' @seealso [fit_latpl()], [compute_var()]
#' @keywords internal
check_sa <- function(FIT) {
  FIT$stoFit$check <- list(done = FALSE, pd = NA, distance = NA, noise = NA)
  if (FIT$dims$npar > 2000) {
    if (!FIT$stoFit$plateau) {
      FIT$stoFit$converged <- FALSE
      warning("The Adam phase did not reach a plateau of the pairwise likelihood: increase ADAM_EPOCHS.", call. = FALSE)
    }
    return(FIT)
  }
  vars <- tryCatch(compute_var(FIT), latpl_not_pd = function(e) NULL)
  FIT$stoFit$check$done <- TRUE
  FIT$stoFit$check$pd <- !is.null(vars)
  if (is.null(vars)) {
    FIT$stoFit$converged <- FALSE
    warning(
      "The observed information is not positive definite at the SA estimate, which is not an optimum of the pairwise likelihood: use set_saqn_control() or other starting values.",
      call. = FALSE
    )
    return(FIT)
  }
  relative_noise <- function(v) {
    se <- sqrt(v$asymptotic_variance / FIT$dims$n)
    max((sqrt(pmax(v$optimisation_noise, 0)) / se)[se > 0])
  }
  FIT$stoFit$check$distance <- max(vars$distance)
  FIT$stoFit$check$noise <- relative_noise(vars)
  distance <- max(vars$distance)
  if (distance > max(1, 4 * relative_noise(vars))) {
    FIT$stoFit$converged <- FALSE
    warning(
      "The SA estimate is ",
      signif(distance, 3),
      " standard errors from the optimum of the pairwise likelihood on its worst parameter, more than the optimisation noise explains: use other starting values or set_saqn_control().",
      call. = FALSE
    )
  } else if (distance > 1) {
    warning(
      "The SA estimate is ",
      signif(distance, 3),
      " standard errors from the optimum of the pairwise likelihood on its worst parameter, because of the optimisation noise, which the standard errors account for. For more precision, increase AVERAGE_EPOCHS or use set_saqn_control(); see sa_diagnostics().",
      call. = FALSE
    )
  }
  FIT$vars <- vars
  FIT
}

#' Standardised covariates for the optimisers
#'
#' Covariates far from zero make the thresholds and their effects almost collinear, and the optimisers unstable. The
#' optimisers then work on standardised covariates, with the matching thresholds and intercepts, which is the same
#' model. This needs thresholds, intercepts of numeric items and effects of covariates free of constraints.
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param LINCON Output of [lincon()].
#' @param NTHR Number of thresholds in each group.
#' @return A list with the data `data` and model `model` for the optimisers, and the maps `to` and `from` the
#'   parameter vector of the optimisers. Without standardisation, `DATA`, no `model` and identities.
#' @seealso [fit_latpl()]
#' @keywords internal
standardise_covariates <- function(DATA, MODEL, LINCON, NTHR) {
  raw <- list(data = DATA, model = NULL, to = identity, from = identity)
  if (ncol(DATA$X) == 0) return(raw)
  grouped <- inherits(MODEL, "latpl_groups")
  models <- if (grouped) MODEL$MODELS else list(MODEL)
  if (!any(vapply(models, function(x) ncol(x$GAMMA) > 0, logical(1)))) return(raw)
  centre <- colMeans(DATA$X)
  spread <- apply(DATA$X, 2, stats::sd)
  spread[spread == 0] <- 1
  block <- theta_block(MODEL, NTHR)
  touched <- block %in% c("tau", "gamma", "kappa") | (block == "nu" & any(DATA$numeric))
  fixed_nu <- any(vapply(models, function(x) any(!is.na(x$CONSTRINT[x$NUMERIC])), logical(1)))
  if (fixed_nu || !all(which(touched) %in% LINCON$free_idx)) {
    if (any(abs(centre) > 3 * spread)) {
      message("Constraints on thresholds, intercepts or effects of covariates prevent standardising the covariates: consider centring them.")
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
  std_model <- if (grouped) replace(MODEL, "MODELS", list(std_models)) else std_models[[1]]
  standard <- function(X) sweep(sweep(X, 2, centre, "-"), 2, spread, "/")
  std_data <- DATA
  std_data$X <- standard(DATA$X)
  std_data$strata <- lapply(DATA$strata, standard)
  if (!is.null(DATA$val)) std_data$val$strata <- lapply(DATA$val$strata, standard)
  map <- function(theta, from_model, to_model, to_standard) {
    pars <- theta2par(from_model, theta)
    if (!grouped) pars <- list(pars)
    pars <- Map(function(par, m) {
      if (is.null(par$gamma)) return(par)
      gamma_raw <- if (to_standard) par$gamma else sweep(par$gamma, 2, spread, "/")
      kappa_raw <- if (to_standard) par$kappa else sweep(par$kappa, 2, spread, "/")
      loadings <- par$loadings %*% solve(diag(ncol(par$loadings)) - par$regressions)
      offset <- as.vector((loadings %*% gamma_raw + kappa_raw) %*% centre) * (if (to_standard) 1 else -1)
      par$thresholds <- par$thresholds - rep(offset, DATA$cat - 1)
      par$intercepts[m$NUMERIC] <- par$intercepts[m$NUMERIC] + offset[m$NUMERIC]
      par$gamma <- if (to_standard) scale_x(par$gamma) else gamma_raw
      par$kappa <- if (to_standard) scale_x(par$kappa) else kappa_raw
      par
    }, pars, models)
    par2theta(to_model, if (grouped) pars else pars[[1]])
  }
  list(
    data = std_data,
    model = std_model,
    to = function(theta) map(theta, MODEL, std_model, TRUE),
    from = function(theta) map(theta, std_model, MODEL, FALSE)
  )
}

#' Inputs of the C++ routines
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param TABLES Use the frequency tables (`TRUE`) or the individual observations (`FALSE`). If `NA`, as set by
#'   [set_latpl_data()].
#' @return A list with `model` and `data`, each a list with one element per group: the model specification, and
#'   the frequencies, distinct covariate rows, distinct response patterns with their strata, unit weights and cluster
#'   weights of the group.
#' @seealso [fit_latpl()], [run_sa()]
#' @keywords internal
engine <- function(DATA, MODEL, TABLES = NA) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
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
        CONSTRMAT = x$CONSTRMAT,
        CONSTRLOGSD = x$CONSTRLOGSD,
        CORRFLAG = corrflag_code(x),
        REGMAT = x$REGMAT,
        RESCOV = x$RESCOV,
        CONSTRINT = x$CONSTRINT,
        CONSTRMEAN = x$CONSTRMEAN,
        CONSTRLOGSCALE = x$CONSTRLOGSCALE,
        THETAPAR = x$PARAMETERIZATION == "theta",
        GAMMA = if (ncol(x$GAMMA) == nx) x$GAMMA else matrix(0, x$q, nx),
        KAPPA = if (ncol(x$KAPPA) == nx) x$KAPPA else matrix(0, x$p, nx)
      )
    }),
    data = lapply(seq_len(DATA$ngroups) - 1, function(g) {
      rows <- which(DATA$group == g)
      key <- if (any(DATA$numeric)) seq_along(rows) else do.call(paste, as.data.frame(cbind(DATA$stratum[rows], DATA$Y[rows, , drop = FALSE])))
      first <- !duplicated(key)
      id <- match(key, key[first])
      WM <- matrix(0, 0, 0)
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
        WM <- (A - rowMeans(A)) * sqrt(C / (C - 1))
      }
      tables <- if (is.na(TABLES)) DATA$tables[g + 1] else TABLES
      freq <- if (tables == DATA$tables[g + 1]) {
        if (tables) DATA$freq[ot[g + 1] * nc + seq_len(ns[g + 1] * nc)] else numeric(0)
      } else if (tables) {
        pairs_freq(DATA$Y[rows, , drop = FALSE], DATA$cat, if (length(DATA$weights)) DATA$weights[rows] else numeric(0), DATA$stratum[rows])
      } else {
        numeric(0)
      }
      list(
        N = length(rows),
        C_VEC = DATA$cat,
        TABLES = tables,
        FREQ = freq,
        UFREQ = DATA$ufreq[ou[g + 1] + seq_len(ns[g + 1] * uc)],
        X = t(DATA$strata[[g + 1]]),
        S = DATA$stratum[rows[first]],
        W = as.numeric(rowsum(gw[rows], id)),
        Y = DATA$Y[rows[first], , drop = FALSE],
        UW = DATA$uweights[rows[first]],
        WV = sqrt(as.numeric(rowsum(jw[rows]^2, id))),
        WM = WM
      )
    })
  )
}

#' Starting values
#'
#' @param INIT Initialisation object from [set_init()].
#' @param DATA Data object from [set_latpl_data()].
#' @param DIMS Dimensions from [latpl_dims()].
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param VERBOSE If `TRUE`, report user-supplied starting values.
#' @return The starting parameter vector: `INIT$THETA` when supplied, otherwise thresholds from the margins of each
#'   group and fixed values for the other parameters.
#' @seealso [fit_latpl()], [spec_from_lavaan()]
#' @keywords internal
init_theta <- function(INIT, DATA, DIMS, MODEL, VERBOSE) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
  nthr <- sum(DATA$cat - 1)
  theta <- INIT$THETA
  if (is.null(theta)) {
    w <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
    theta <- unlist(lapply(seq_along(models), function(g) {
      m <- models[[g]]
      rows <- DATA$group == g - 1
      block <- theta_block(m, nthr)
      th <- rep(0, length(block))
      th[block == "tau"] <- unlist(lapply(seq_len(DIMS$p), function(j) {
        pj <- vapply(
          seq_len(DIMS$cat[j]) - 1,
          function(a) sum(w[rows & DATA$Y[, j] %in% a]),
          numeric(1)
        )
        stats::qnorm(pmin(pmax(cumsum(pj) / sum(pj), 1e-4), 1 - 1e-4))[
          -DIMS$cat[j]
        ]
      }))
      th[block == "beta"] <- 0.25
      spread <- rep(1, DIMS$p)
      centre <- rep(0, DIMS$p)
      for (j in which(DATA$numeric)) {
        seen <- rows & !is.na(DATA$Y[, j])
        centre[j] <- stats::weighted.mean(DATA$Y[seen, j], w[seen])
        spread[j] <- sqrt(stats::weighted.mean((DATA$Y[seen, j] - centre[j])^2, w[seen]))
      }
      start <- signed_loadings(DATA$Y[rows, , drop = FALSE], m, spread)
      th[block == "lambda"] <- start$loadings[is.na(m$CONSTRMAT)]
      th[block == "logsd"] <- log(start$sd[is.na(m$CONSTRLOGSD)])
      th[block == "nu"] <- centre[is.na(m$CONSTRINT)]
      th[block == "scale"][m$NUMERIC[is.na(m$CONSTRLOGSCALE)]] <- log(spread[m$NUMERIC & is.na(m$CONSTRLOGSCALE)])
      if (m$PARAMETERIZATION == "theta") {
        par <- theta2par(m, th)
        L <- par$loadings %*% solve(diag(m$q) - par$regressions)
        th[block == "tau"] <- th[block == "tau"] *
          rep(
            sqrt(
              rowSums((L %*% par$latent_correlations) * L) +
                diag(par$residual_covariances)
            ),
            DIMS$cat - 1
          )
      }
      th
    }))
  } else if (VERBOSE) {
    message("- Initialised at user-supplied values.")
  }
  stopifnot(is.numeric(theta), length(theta) == DIMS$d, all(is.finite(theta)))
  if (!is.null(INIT$THETA) && implied_variance(MODEL, theta) > 1) {
    warning("Negative error variances implied by init parameters")
  }
  theta
}

#' Signed starting loadings
#'
#' The first eigenvector of the correlations among the items of each factor, so that reversed items start with a
#' negative loading. Factors scaled by a fixed loading take the variance that matches it, and factors without items keep
#' a small variance.
#'
#' @param Y Data matrix of the group.
#' @param MODEL Model object from [set_latpl_model()].
#' @param SPREAD Standard deviation of numeric items, one for the others.
#' @return A list with the matrix `loadings` and the latent standard deviations `sd`.
#' @seealso [init_theta()]
#' @keywords internal
signed_loadings <- function(Y, MODEL, SPREAD) {
  rows <- unique(round(seq(1, nrow(Y), length.out = min(nrow(Y), 2000))))
  L <- MODEL$CONSTRMAT
  sd <- ifelse(is.na(MODEL$CONSTRLOGSD), .25, exp(MODEL$CONSTRLOGSD))
  L[is.na(L)] <- .25
  for (h in seq_len(MODEL$q)) {
    S <- which(is.na(MODEL$CONSTRMAT[, h]) | MODEL$CONSTRMAT[, h] != 0)
    if (length(S) < 2) next
    R <- suppressWarnings(stats::cor(Y[rows, S, drop = FALSE], use = "pairwise.complete.obs"))
    R[is.na(R)] <- 0
    e <- eigen(R, symmetric = TRUE)
    v <- pmin(pmax(e$vectors[, 1] * sqrt(max(e$values[1] - 1, .1)), -.9), .9)
    v <- if (v[1] < 0) -v else v
    fix <- which(!is.na(MODEL$CONSTRMAT[S, h]))[1]
    if (is.na(MODEL$CONSTRLOGSD[h]) && !is.na(fix)) {
      sd[h] <- abs(v[fix] * SPREAD[S[fix]] / MODEL$CONSTRMAT[S[fix], h])
      v <- v * sign(v[fix] * MODEL$CONSTRMAT[S[fix], h]) * SPREAD[S] / sd[h]
    } else {
      if (is.na(MODEL$CONSTRLOGSD[h])) sd[h] <- 1
      if (!is.na(fix)) v <- v * sign(v[fix] * MODEL$CONSTRMAT[S[fix], h])
      v <- v * SPREAD[S] / sd[h]
    }
    free <- is.na(MODEL$CONSTRMAT[S, h])
    L[S[free], h] <- v[free]
  }
  list(loadings = L, sd = sd)
}

#' Largest variance of the latent responses implied by the factors
#'
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param THETA Parameter vector.
#' @return The largest implied variance among the ordinal items under `PARAMETERIZATION = "delta"`, which must be below
#'   one, and zero if there are none.
#' @seealso [fit_latpl()], [init_theta()]
#' @keywords internal
implied_variance <- function(MODEL, THETA) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
  pars <- if (inherits(MODEL, "latpl_groups")) theta2par(MODEL, THETA) else list(theta2par(MODEL, THETA))
  max(0, unlist(Map(function(m, par) {
    if (m$PARAMETERIZATION != "delta") return(0)
    L <- par$loadings %*% solve(diag(m$q) - par$regressions)
    (rowSums((L %*% par$latent_correlations) * L) * par$scales^2)[!m$NUMERIC]
  }, models, pars)))
}

#' Stochastic approximation
#'
#' @param START Starting parameter vector.
#' @param DATA Data object from [set_latpl_data()].
#' @param ENGINE Output of [engine()].
#' @param LINCON Output of [lincon()].
#' @param CONTROL Control object from [set_sa_control()] or [set_sa_init_control()].
#' @param ADAM_ONLY If `TRUE`, stop after the Adam phase.
#' @param VERBOSE If `TRUE`, print progress.
#' @return A list with the estimate `theta`, the last and end-of-burn-in iterates `theta_last` and `theta_burn`,
#'   the epoch length and phase counts, the rejected and skipped steps, the trajectory `path` and the timings
#'   `time`.
#' @seealso [fit_latpl()], [set_sa_control()]
#' @keywords internal
run_sa <- function(START, DATA, ENGINE, LINCON, CONTROL, ADAM_ONLY, VERBOSE) {
  valdata <- ENGINE$data
  if (!is.null(DATA$val)) {
    valdata <- lapply(seq_len(DATA$ngroups) - 1, function(g) {
      rows <- DATA$val$group == g
      list(
        TABLES = DATA$val$tables[g + 1],
        FREQ = if (DATA$val$tables[g + 1]) DATA$val$freq[[g + 1]] else numeric(0),
        UFREQ = numeric(0),
        X = t(DATA$val$strata[[g + 1]]),
        S = DATA$val$stratum[rows],
        W = rep(1, sum(rows)),
        Y = DATA$val$Y[rows, , drop = FALSE]
      )
    })
  }
  npar <- length(LINCON$free_idx)
  free <- npar < length(START)
  pairs <- DATA$ngroups * DATA$p * (DATA$p - 1) / 2
  adam_epochs <- CONTROL$ADAM_EPOCHS
  if (is.null(adam_epochs)) adam_epochs <- max(5L, ceiling(150 / ceiling(pairs / min(CONTROL$PAIRS_PER_ITERATION, pairs))))
  gain <- switch(
    CONTROL$GAIN,
    auto = as.integer(npar <= 2000),
    hessian = 1L,
    diagonal = 0L
  )
  out <- cpp_sa(
    START = START[LINCON$free_idx],
    K = if (free) LINCON$K else matrix(0, 0, 0),
    C = LINCON$c,
    MODEL = ENGINE$model,
    DATA = ENGINE$data,
    VALDATA = valdata,
    PAIRS_PER_ITERATION = CONTROL$PAIRS_PER_ITERATION,
    ALPHA = CONTROL$ALPHA,
    ADAM_EPOCHS = as.integer(adam_epochs),
    PRECONDITION_EPOCHS = CONTROL$PRECONDITION_EPOCHS,
    MAX_PRECONDITION_EPOCHS = CONTROL$MAX_PRECONDITION_EPOCHS,
    AVERAGE_EPOCHS = CONTROL$AVERAGE_EPOCHS,
    TOL = CONTROL$TOL,
    DECAY = CONTROL$DECAY,
    AVG_FROM = CONTROL$AVG_FROM,
    GAIN = gain,
    SEED = CONTROL$SEED,
    ADAM_ONLY = ADAM_ONLY,
    VERBOSE = VERBOSE
  )
  out$path$nll <- c(out$path$nll, rep(NA, sum(out$path$phase == 3)))
  out
}
