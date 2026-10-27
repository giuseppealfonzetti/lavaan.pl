#'Fit factor models for ordinal data with pairwise likelihood methods
#'
#'@description
#'
#'Fit latent variable models for ordinal variable using pairwise likelihood
#'methods running via quasi-Newton BFGS or stochastic approximations.
#'
#'@param DATA Data object from [set_latpl_data()].
#'@param MODEL Model object from [set_latpl_model()].
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
){
  start_time <- Sys.time()
  stopifnot(
    inherits(DATA, "latpl_data"),
    inherits(MODEL, "latpl_model"),
    inherits(CONTROL, c("latpl_ucminf_control", "latpl_sa_control")),
    inherits(INIT, "latpl_init"),
    DATA$p == MODEL$p
  )
  METHOD <- if (inherits(CONTROL, "latpl_sa_control")) "SA" else "ucminf"
  dims <- latpl_dims(DATA, MODEL)
  lc <- lincon(MODEL, dims$nthr)
  if (METHOD == "SA") {
    if (dims$p < 40) warning("SA defaults are tuned for models with at least 40 items: check the estimates against set_ucminf_control().", call. = FALSE)
    if (corrflag_code(MODEL) == 2) stop("LATCOV = \"direct\" is not currently supported with SA; use set_ucminf_control().", call. = FALSE)
    if (!sa_model(MODEL)) stop("REGMAT and RESCOV are not currently supported with SA; use set_ucminf_control().", call. = FALSE)
    if (length(lc$free_idx) < dims$d) stop("CONSTRAINTS are not currently supported with SA; use set_ucminf_control().", call. = FALSE)
  }
  dims$npar <- length(lc$free_idx)
  RcppParallel::setThreadOptions(numThreads = NCORES)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  eng <- engine(DATA, MODEL)
  fit <- structure(list(data = DATA, model = MODEL, control = CONTROL, dims = dims, method = METHOD, lincon = lc,
                        engine = eng, version = as.character(utils::packageVersion("lavaan.pl")), cores = NCORES),
                   class = "latpl_fit")
  fit$init <- init_theta(INIT, dims, MODEL, eng, VERBOSE)
  full <- function(x) as.vector(lc$K %*% x + lc$c)
  nll <- function(theta) -cpp_pl(theta, eng$model, eng$data)$pl / dims$n

  # With the direct parametrisation, points with a non positive definite latent
  # covariance matrix are outside the parameter space
  pd <- function(theta) {
    corrflag_code(MODEL) != 2 || !inherits(try(chol(theta2par(MODEL, theta)$latent_correlations), silent = TRUE), "try-error")
  }

  # Direct parametrisation: warm start from the LKJ fit without the constraints on
  # the latent covariances, which is better conditioned
  if (METHOD == "ucminf" && corrflag_code(MODEL) == 2 && pd(fit$init)) {
    block <- theta_block(MODEL, dims$nthr)
    warm_model <- MODEL
    warm_model$LATCOV <- "lkj"
    if (!is.null(MODEL$CONSTRAINTS)) {
      keep <- rowSums(MODEL$CONSTRAINTS$A[, block == "corr", drop = FALSE] != 0) == 0
      warm_model$CONSTRAINTS <- if (any(keep)) list(A = MODEL$CONSTRAINTS$A[keep, , drop = FALSE], b = MODEL$CONSTRAINTS$b[keep])
    }
    warm <- fit_latpl(DATA, warm_model, CONTROL, set_init(THETA = par2theta(warm_model, theta2par(MODEL, fit$init))), NCORES)
    RcppParallel::setThreadOptions(numThreads = NCORES)
    theta <- par2theta(MODEL, get_par(warm, OPTION = "list"))
    # Least-squares projection on the constraints, shrinking the free latent
    # covariances until positive definite
    x <- as.vector(solve(crossprod(lc$K), crossprod(lc$K, theta - lc$c)))
    cov_free <- block[lc$free_idx] == "corr"
    for (i in 1:20) {
      theta <- full(x)
      if (all(is.finite(theta)) && pd(theta)) break
      x[cov_free] <- 0.8 * x[cov_free]
    }
    start <- full(fit$init[lc$free_idx])
    if (all(is.finite(theta)) && pd(theta) && isTRUE(nll(theta) < if (pd(start)) nll(start) else Inf)) fit$init <- theta
  }
  fit$init <- full(fit$init[lc$free_idx])

  if (METHOD == "ucminf") {
    if (VERBOSE) message("- Optimising with ucminf...")
    fn <- function(x) {
      theta <- full(x)
      val <- if (pd(theta)) nll(theta) else Inf
      if (is.finite(val)) val else Inf
    }
    gr <- function(x) {
      theta <- full(x)
      g <- if (pd(theta)) -as.vector(crossprod(lc$K, cpp_pl(theta, eng$model, eng$data, 1L)$gradient)) / dims$n else 0 * x
      if (all(is.finite(g))) g else 0 * x
    }
    fit$numFit <- ucminf::ucminf(fit$init[lc$free_idx], fn, gr, hessian = 0,
                                 control = list(grtol = CONTROL$GRTOL, xtol = CONTROL$XTOL, stepmax = CONTROL$STEPMAX, maxeval = CONTROL$MAXEVAL))
    fit$theta <- full(fit$numFit$par)
  } else {
    if (VERBOSE) message("- Optimising via pairwise stochastic approximation...")
    if (is.null(CONTROL$STEP0)) CONTROL$STEP0 <- sqrt(1 / dims$d)
    fit$control <- CONTROL
    val <- is.null(DATA$valfreq)
    opt <- run_sa(fit$init, eng, CONTROL, if (val) DATA$freq else DATA$valfreq, if (val) DATA$ufreq else numeric(0), VERBOSE)
    if (VERBOSE) {
      if (!opt$convergence_burn) message("Burn-in period did not reach tolerance level ", CONTROL$TOL_BURN, ". Stopped at ", CONTROL$BURNE, " epochs Try increasing BURNE or STEP0.")
      if (opt$proj_after_burn) message("Projections performed after the burn-in. Trajectories might be unstable. Check trajectories visually and evaluate decreasing STEP0.")
      if (opt$convergence_full == 0) message("Algorithm ended before reaching tolerance level ", CONTROL$TOL_END, ". Stopped at ", CONTROL$MAXE, ". epochs Try increasing MAXE or STEP0.")
      if (opt$convergence_full == -1) message("Stopped because of divergent trajectories. Decrease STEP0.")
      if (!opt$finite_gr) message("The algorithm automatically dealt with divergent trajectories. Check trajectories visually and evaluate decreasing STEP0.")
      if (opt$neg_pdiff) message("Possible divergent trajectories detected monitoring the full negative pairwise likelihood. Check trajectories visually and evaluate decreasing STEP0.")
    }
    fit$stoFit <- list(convergence = opt[c("convergence_full", "convergence_burn", "proj_after_burn", "finite_gr", "neg_pdiff")],
                       path = list(iters = opt$path_iters, theta = opt$path_theta, avtheta = opt$path_avtheta, nll = opt$path_nll),
                       nll = opt$nll, last_iter = opt$last_iter, burnt = opt$burnt)
    fit$theta <- opt$avtheta
  }
  fit$RTime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (VERBOSE) message("Done! (", round(fit$RTime, 2), " secs)")
  fit
}

# Lists passed to the C++ engine. The units are collapsed to their distinct response
# patterns, with the square roots of the summed squared score weights WV or, with
# clusters, the centred and scaled cluster aggregation matrix WM
engine <- function(DATA, MODEL) {
  key <- do.call(paste, as.data.frame(DATA$Y))
  first <- !duplicated(key)
  id <- match(key, key[first])
  gw <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
  jw <- if (length(DATA$jweights)) DATA$jweights else rep(1, DATA$n)
  WM <- matrix(0, 0, 0)
  if (length(DATA$cluster)) {
    G <- length(unique(DATA$cluster))
    A <- tapply(gw, list(factor(id, seq_len(sum(first))), factor(DATA$cluster, seq_len(G) - 1)), sum, default = 0)
    WM <- (A - rowMeans(A)) * sqrt(G / (G - 1))
  }
  list(model = list(CONSTRMAT = MODEL$CONSTRMAT, CONSTRLOGSD = MODEL$CONSTRLOGSD, CORRFLAG = corrflag_code(MODEL),
                    REGMAT = MODEL$REGMAT, RESCOV = MODEL$RESCOV),
       data = list(N = DATA$n, C_VEC = DATA$cat, FREQ = DATA$freq, UFREQ = DATA$ufreq, Y = DATA$Y[first, , drop = FALSE],
                   UW = DATA$uweights[first], WV = sqrt(as.numeric(rowsum(jw^2, id))), WM = WM))
}

init_theta <- function(INIT, DIMS, MODEL, ENGINE, VERBOSE) {
  proj <- function(theta) if (sa_model(MODEL)) cpp_sa_proj(theta, ENGINE$model, ENGINE$data) else theta
  theta <- INIT$THETA
  if (is.null(theta)) {
    thr <- unlist(lapply(DIMS$cat - 1, function(m) if (m == 1) 0 else seq(-1, 1, length.out = m)))
    theta <- proj(c(thr, rep(0.25, DIMS$nload), rep(0, DIMS$d - DIMS$nthr - DIMS$nload)))
  } else if (VERBOSE) {
    message("- Initialised at user-supplied values.")
  }
  stopifnot(is.numeric(theta), length(theta) == DIMS$d, all(is.finite(theta)))
  if (INIT$METHOD == "SA" && DIMS$p >= 15 && sa_model(MODEL)) {
    if (VERBOSE) message("- Initialising with SA...")
    CONTROL <- INIT$CONTROL
    if (is.null(CONTROL$STEP0)) CONTROL$STEP0 <- sqrt(1 / DIMS$d)
    theta <- run_sa(theta, ENGINE, CONTROL, ENGINE$data$FREQ, ENGINE$data$UFREQ, FALSE)$avtheta
  }
  par <- theta2par(MODEL, theta)
  L <- par$loadings %*% solve(diag(DIMS$q) - par$regressions)
  if (any(rowSums((L %*% par$latent_correlations) * L) > 1)) warning("Negative error variances implied by init parameters")
  proj(theta)
}

run_sa <- function(THETA, ENGINE, CONTROL, VALFREQ, VALUFREQ, VERBOSE) {
  do.call(cpp_plSA2, c(list(THETA_INIT = THETA, MODEL = ENGINE$model, DATA = ENGINE$data, VALFREQ = VALFREQ,
                            VALUFREQ = VALUFREQ), unclass(CONTROL), VERBOSE = VERBOSE))
}
