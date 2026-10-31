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
  dims$npar <- length(lc$free_idx)
  RcppParallel::setThreadOptions(numThreads = NCORES)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  eng <- engine(DATA, MODEL)
  fit <- structure(list(data = DATA, model = MODEL, control = CONTROL, dims = dims, method = METHOD, lincon = lc,
                        engine = eng, version = as.character(utils::packageVersion("lavaan.pl")), cores = NCORES),
                   class = "latpl_fit")
  fit$init <- init_theta(INIT, DATA, dims, MODEL, VERBOSE)
  if (METHOD == "ucminf" && INIT$METHOD == "SA") {
    if (VERBOSE) message("- Initialising with the SA Adam phase...")
    fit$init <- run_sa(list(fit$init), DATA, eng, lc, INIT$CONTROL, TRUE, FALSE)$theta
  }
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
    opt <- run_sa(list(fit$init), DATA, eng, lc, CONTROL, FALSE, VERBOSE)
    if (!opt$plateau) warning("The SA burn-in did not reach a plateau of the pairwise likelihood: increase BURN_EPOCHS, or check the estimates against set_ucminf_control().", call. = FALSE)
    if (opt$skipped_averaging > 0) warning(opt$skipped_averaging, " SA iterations of the averaged phase were skipped because no feasible step was found.", call. = FALSE)
    opt$converged <- opt$plateau && opt$skipped_averaging == 0 && all(is.finite(opt$theta))
    fit$stoFit <- opt[setdiff(names(opt), "theta")]
    fit$theta <- opt$theta
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
                    REGMAT = MODEL$REGMAT, RESCOV = MODEL$RESCOV, CONSTRINT = MODEL$CONSTRINT, CONSTRMEAN = MODEL$CONSTRMEAN,
                    CONSTRLOGSCALE = MODEL$CONSTRLOGSCALE, THETAPAR = MODEL$PARAMETERIZATION == "theta"),
       data = list(N = DATA$n, C_VEC = DATA$cat, FREQ = DATA$freq, UFREQ = DATA$ufreq, Y = DATA$Y[first, , drop = FALSE],
                   UW = DATA$uweights[first], WV = sqrt(as.numeric(rowsum(jw^2, id))), WM = WM))
}

init_theta <- function(INIT, DATA, DIMS, MODEL, VERBOSE) {
  theta <- INIT$THETA
  if (is.null(theta)) {
    w <- if (length(DATA$weights)) DATA$weights else rep(1, DATA$n)
    thr <- unlist(lapply(seq_len(DIMS$p), function(j) {
      pj <- vapply(seq_len(DIMS$cat[j]) - 1, function(a) sum(w[DATA$Y[, j] %in% a]), numeric(1))
      stats::qnorm(pmin(pmax(cumsum(pj) / sum(pj), 1e-4), 1 - 1e-4))[-DIMS$cat[j]]
    }))
    theta <- c(thr, rep(0.25, DIMS$nload), rep(0, DIMS$d - DIMS$nthr - DIMS$nload))
    block <- theta_block(MODEL, DIMS$nthr)
    theta[block == "beta"] <- 0.25
    # Factors identified by a marker: the image of the standardised start (all loadings 0.25)
    marker <- is.na(MODEL$CONSTRLOGSD)
    theta[block == "lambda"][col(MODEL$CONSTRMAT)[is.na(MODEL$CONSTRMAT)] %in% which(marker)] <- 1
    theta[block == "logsd"] <- log(0.25)
    if (MODEL$PARAMETERIZATION == "theta") {
      par <- theta2par(MODEL, theta)
      L <- par$loadings %*% solve(diag(DIMS$q) - par$regressions)
      theta[block == "tau"] <- theta[block == "tau"] *
        rep(sqrt(rowSums((L %*% par$latent_correlations) * L) + diag(par$residual_covariances)), DIMS$cat - 1)
    }
  } else if (VERBOSE) {
    message("- Initialised at user-supplied values.")
  }
  stopifnot(is.numeric(theta), length(theta) == DIMS$d, all(is.finite(theta)))
  par <- theta2par(MODEL, theta)
  L <- par$loadings %*% solve(diag(DIMS$q) - par$regressions)
  if (MODEL$PARAMETERIZATION == "delta" && any(rowSums((L %*% par$latent_correlations) * L) * par$scales^2 > 1)) {
    warning("Negative error variances implied by init parameters")
  }
  theta
}

run_sa <- function(STARTS, DATA, ENGINE, LINCON, CONTROL, ADAM_ONLY, VERBOSE) {
  val <- is.null(DATA$valfreq)
  npar <- length(LINCON$free_idx)
  free <- npar < length(STARTS[[1]])
  gain <- switch(CONTROL$GAIN, auto = as.integer(npar <= 2000), hessian = 1L, diagonal = 0L)
  out <- cpp_sa(STARTS = matrix(vapply(STARTS, function(x) x[LINCON$free_idx], numeric(npar)), npar),
                K = if (free) LINCON$K else matrix(0, 0, 0), C = LINCON$c,
                MODEL = ENGINE$model, DATA = ENGINE$data, VALFREQ = if (val) DATA$freq else DATA$valfreq,
                VALUFREQ = if (val) DATA$ufreq else numeric(0), PAIRS_PER_ITERATION = CONTROL$PAIRS_PER_ITERATION,
                ALPHA = CONTROL$ALPHA, ADAM_EPOCHS = CONTROL$ADAM_EPOCHS, BURN_EPOCHS = CONTROL$BURN_EPOCHS,
                EPOCHS = CONTROL$EPOCHS, TOL = CONTROL$TOL, DECAY = CONTROL$DECAY, AVG_FROM = CONTROL$AVG_FROM,
                GAIN = gain, SEED = CONTROL$SEED, ADAM_ONLY = ADAM_ONLY, VERBOSE = VERBOSE)
  out$path$nll <- c(out$path$nll, rep(NA, sum(out$path$phase == 3)))
  out
}
