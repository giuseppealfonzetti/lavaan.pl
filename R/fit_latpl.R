#'Fit factor models for ordinal data with pairwise likelihood methods
#'
#'@description
#'
#'Fit latent variable models for ordinal variable using pairwise likelihood
#'methods running via quasi-Newton BFGS or stochastic approximations.
#'
#'@param DATA Data object from [set_latpl_data()].
#'@param MODEL Model object from [set_latpl_model()], or from [set_latpl_groups()] for several groups.
#'@param CONTROL Control object of the main optimiser: [set_ucminf_control()]
#'  for the quasi-Newton BFGS optimiser from the \code{ucminf} package,
#'  [set_sa_control()] for Stochastic Approximations.
#'@param INIT Starting values object from [set_init()].
#'@param NCORES Integer value setting the number of threads to carry out the
#'  estimation.
#'@param VERBOSE `TRUE` for verbose output
#'
#'@return Object of class `latpl_fit`: a list with the inputs `data`, `model` and
#'  `control`, the starting values `init`, the estimates `theta` (use
#'  [get_par()]), the dimensions `dims`, the optimiser output `numFit` (ucminf)
#'  or `stoFit` (SA), and the timing `RTime`.
#'@export
fit_latpl <- function(
  DATA,
  MODEL,
  CONTROL = set_ucminf_control(),
  INIT = set_init(),
  NCORES = 1,
  VERBOSE = FALSE
) {
  start_time <- Sys.time()
  stopifnot(
    inherits(DATA, "latpl_data"),
    inherits(MODEL, c("latpl_model", "latpl_groups")),
    inherits(CONTROL, c("latpl_ucminf_control", "latpl_sa_control")),
    inherits(INIT, "latpl_init"),
    DATA$p == MODEL$p
  )
  grouped <- inherits(MODEL, "latpl_groups")
  models <- if (grouped) MODEL$MODELS else list(MODEL)
  if (DATA$ngroups != length(models)) {
    stop("DATA has ", DATA$ngroups, " groups, MODEL has ", length(models), ".")
  }
  direct <- vapply(models, function(x) corrflag_code(x) == 2, logical(1))
  nthr <- sum(DATA$cat - 1)
  METHOD <- if (inherits(CONTROL, "latpl_sa_control")) "SA" else "ucminf"
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
  fit$init <- init_theta(INIT, DATA, dims, MODEL, VERBOSE)
  if (METHOD == "ucminf" && INIT$METHOD == "SA") {
    if (VERBOSE) {
      message("- Initialising with the SA Adam phase...")
    }
    fit$init <- run_sa(
      list(fit$init),
      DATA,
      eng,
      lc,
      INIT$CONTROL,
      TRUE,
      FALSE
    )$theta
  }
  full <- function(x) as.vector(lc$K %*% x + lc$c)
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

  if (METHOD == "ucminf" && any(direct) && pd(fit$init)) {
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

  if (METHOD == "ucminf") {
    if (VERBOSE) {
      message("- Optimising with ucminf...")
    }
    fn <- function(x) {
      theta <- full(x)
      val <- if (pd(theta)) nll(theta) else Inf
      if (is.finite(val)) val else Inf
    }
    gr <- function(x) {
      theta <- full(x)
      g <- if (pd(theta)) {
        -as.vector(crossprod(
          lc$K,
          cpp_pl(theta, eng$model, eng$data, 1L)$gradient
        )) /
          dims$n
      } else {
        0 * x
      }
      if (all(is.finite(g))) g else 0 * x
    }
    fit$numFit <- ucminf::ucminf(
      fit$init[lc$free_idx],
      fn,
      gr,
      hessian = 0,
      control = list(
        grtol = CONTROL$GRTOL,
        xtol = CONTROL$XTOL,
        stepmax = CONTROL$STEPMAX,
        maxeval = CONTROL$MAXEVAL
      )
    )
    fit$theta <- full(fit$numFit$par)
  } else {
    if (VERBOSE) {
      message("- Optimising via pairwise stochastic approximation...")
    }
    opt <- run_sa(list(fit$init), DATA, eng, lc, CONTROL, FALSE, VERBOSE)
    if (!opt$plateau) {
      warning(
        "The SA burn-in did not reach a plateau of the pairwise likelihood: increase BURN_EPOCHS, or check the estimates against set_ucminf_control().",
        call. = FALSE
      )
    }
    if (opt$skipped_averaging > 0) {
      warning(
        opt$skipped_averaging,
        " SA iterations of the averaged phase were skipped because no feasible step was found.",
        call. = FALSE
      )
    }
    opt$converged <- opt$plateau &&
      opt$skipped_averaging == 0 &&
      all(is.finite(opt$theta))
    fit$stoFit <- opt[setdiff(names(opt), "theta")]
    fit$theta <- opt$theta
  }
  fit$RTime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (VERBOSE) {
    message("Done! (", round(fit$RTime, 2), " secs)")
  }
  fit
}


#' Inputs of the C++ routines
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @return A list with `model` and `data`, each a list with one element per group: the model specification, and
#'   the frequencies, distinct response patterns, unit weights and cluster weights of the group.
#' @seealso [fit_latpl()], [run_sa()]
#' @keywords internal
engine <- function(DATA, MODEL) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
  gw <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
  jw <- if (length(DATA$jweights)) DATA$jweights else rep(1, DATA$n)
  nc <- length(DATA$freq) / DATA$ngroups
  nu <- length(DATA$ufreq) / DATA$ngroups
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
        THETAPAR = x$PARAMETERIZATION == "theta"
      )
    }),
    data = lapply(seq_len(DATA$ngroups) - 1, function(g) {
      rows <- which(DATA$group == g)
      key <- do.call(paste, as.data.frame(DATA$Y[rows, , drop = FALSE]))
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
      list(
        N = length(rows),
        C_VEC = DATA$cat,
        FREQ = DATA$freq[g * nc + seq_len(nc)],
        UFREQ = DATA$ufreq[g * nu + seq_len(nu)],
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
      th[block %in% c("lambda", "beta")] <- 0.25
      fixed <- if (length(m$REGMAT)) {
        rbind(m$CONSTRMAT, m$REGMAT)
      } else {
        m$CONSTRMAT
      }
      marker <- is.na(m$CONSTRLOGSD) & colSums(!is.na(fixed) & fixed != 0) > 0
      th[block == "lambda"][
        col(m$CONSTRMAT)[is.na(m$CONSTRMAT)] %in% which(marker)
      ] <- 1
      th[block == "logsd"] <- ifelse(marker[is.na(m$CONSTRLOGSD)], log(0.25), 0)
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
  par <- if (inherits(MODEL, "latpl_groups")) {
    theta2par(MODEL, theta)
  } else {
    list(theta2par(MODEL, theta))
  }
  for (g in seq_along(models)) {
    L <- par[[g]]$loadings %*% solve(diag(models[[g]]$q) - par[[g]]$regressions)
    if (
      models[[g]]$PARAMETERIZATION == "delta" &&
        any(
          rowSums((L %*% par[[g]]$latent_correlations) * L) *
            par[[g]]$scales^2 >
            1
        )
    ) {
      warning("Negative error variances implied by init parameters")
    }
  }
  theta
}

#' Stochastic approximation
#'
#' @param STARTS List of starting parameter vectors.
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
run_sa <- function(STARTS, DATA, ENGINE, LINCON, CONTROL, ADAM_ONLY, VERBOSE) {
  val <- is.null(DATA$valfreq)
  npar <- length(LINCON$free_idx)
  free <- npar < length(STARTS[[1]])
  gain <- switch(
    CONTROL$GAIN,
    auto = as.integer(npar <= 2000),
    hessian = 1L,
    diagonal = 0L
  )
  out <- cpp_sa(
    STARTS = matrix(
      vapply(STARTS, function(x) x[LINCON$free_idx], numeric(npar)),
      npar
    ),
    K = if (free) LINCON$K else matrix(0, 0, 0),
    C = LINCON$c,
    MODEL = ENGINE$model,
    DATA = ENGINE$data,
    VALFREQ = if (val) DATA$freq else DATA$valfreq,
    VALUFREQ = if (val) DATA$ufreq else numeric(0),
    PAIRS_PER_ITERATION = CONTROL$PAIRS_PER_ITERATION,
    ALPHA = CONTROL$ALPHA,
    ADAM_EPOCHS = CONTROL$ADAM_EPOCHS,
    BURN_EPOCHS = CONTROL$BURN_EPOCHS,
    EPOCHS = CONTROL$EPOCHS,
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
