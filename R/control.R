#' Ordinal data for pairwise likelihood estimation
#'
#'@param DATA Integer data matrix of dimension \eqn{n*p}. Categories must be
#'  coded starting from zero.
#'@param VALDATA Validation data used to monitor convergence of stochastic
#'  approximations. If `NULL`, data passed via `DATA` is used for monitoring
#'  purposes.
#'@param MISSING Treatment of missing values (`NA`): `"listwise"` drops incomplete observations,
#'  `"pairwise"` uses the complete pairs and `"available.cases"` the available cases of
#'  Katsikatsou et al. (2022).
#'@param WEIGHTS Sampling weights, one per observation, normalised to sum to the number of
#'  observations. If `NULL`, unit weights.
#'@param WEIGHTS_TYPE Type of sampling weights, as in lavaan: `"design"` or `"frequency"`.
#'@param CLUSTER Cluster membership, one per observation, for cluster-robust standard errors. If
#'  `NULL`, observations are independent.
#'
#' @return An object of class `latpl_data` with the data matrix `Y`, its dimensions `n` and `p`,
#'   the number of categories of each item `cat`, the pairwise frequency tables `freq` and
#'   `valfreq`, the treatment of missing values `missing`, the unit weights `uweights` and
#'   univariate frequencies `ufreq` of the available-cases term, the normalised sampling weights
#'   `weights` and their factors on the scores `jweights`, the cluster indices `cluster`, and the
#'   time spent computing the frequencies `freqTime`.
#'
#' @export
set_latpl_data <- function(DATA, VALDATA = NULL, MISSING = c("listwise", "pairwise", "available.cases"),
                          WEIGHTS = NULL, WEIGHTS_TYPE = c("design", "frequency"), CLUSTER = NULL) {
  MISSING <- match.arg(MISSING)
  WEIGHTS_TYPE <- match.arg(WEIGHTS_TYPE)
  stopifnot(is.matrix(DATA), is.numeric(DATA), all(is.finite(DATA) | is.na(DATA)), min(DATA, na.rm = TRUE) == 0)
  Y <- round(DATA)
  stopifnot(
    is.null(WEIGHTS) || (is.numeric(WEIGHTS) && length(WEIGHTS) == nrow(Y) && all(is.finite(WEIGHTS)) && all(WEIGHTS >= 0)),
    is.null(CLUSTER) || (length(CLUSTER) == nrow(Y) && !anyNA(CLUSTER))
  )
  if (MISSING == "listwise" && anyNA(Y)) {
    keep <- stats::complete.cases(Y)
    message(sum(!keep), " incomplete observations removed (MISSING = \"listwise\").")
    Y <- Y[keep, , drop = FALSE]
    if (!is.null(WEIGHTS)) WEIGHTS <- WEIGHTS[keep]
    if (!is.null(CLUSTER)) CLUSTER <- CLUSTER[keep]
  }
  stopifnot(nrow(Y) > 0, colSums(!is.na(Y)) > 0)
  cat <- apply(Y, 2, max, na.rm = TRUE) + 1
  stopifnot(cat > 1)
  weights <- jweights <- cluster <- numeric(0)
  if (!is.null(WEIGHTS)) {
    weights <- WEIGHTS / sum(WEIGHTS) * nrow(Y)
    jweights <- if (WEIGHTS_TYPE == "design") weights else sqrt(weights)
  }
  if (!is.null(CLUSTER)) {
    cluster <- as.numeric(match(CLUSTER, unique(CLUSTER)) - 1)
    stopifnot(length(unique(cluster)) > 1)
  }
  freq_start_time <- Sys.time()
  freq <- pairs_freq(Y, cat, weights)
  uweights <- if (MISSING == "available.cases") rowSums(is.na(Y)) else rep(0, nrow(Y))
  uw <- if (length(weights)) uweights * weights else uweights
  ufreq <- unlist(lapply(seq_len(ncol(Y)), function(j) {
    vapply(seq_len(cat[j]) - 1, function(a) sum(uw[Y[, j] %in% a]), numeric(1))
  }))
  valfreq <- NULL
  if (!is.null(VALDATA)) {
    stopifnot(
      is.matrix(VALDATA),
      is.numeric(VALDATA),
      all(is.finite(VALDATA)),
      ncol(VALDATA) == ncol(Y),
      all(VALDATA >= 0),
      all(t(VALDATA) <= cat - 1)
    )
    valfreq <- pairs_freq(round(VALDATA), cat)
  }
  out <- list(
    Y = Y,
    n = nrow(Y),
    p = ncol(Y),
    cat = cat,
    freq = freq,
    valfreq = valfreq,
    missing = MISSING,
    uweights = uweights,
    ufreq = ufreq,
    weights = weights,
    jweights = jweights,
    cluster = cluster,
    freqTime = as.numeric(difftime(Sys.time(), freq_start_time, units = "secs"))
  )
  class(out) <- "latpl_data"
  out
}

#' Factor model for pairwise likelihood estimation
#'
#' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA`
#'  refers to free loading parameters. Elements set to numerical values denote
#'  fixed values constraints.
#' @param CONSTRVAR \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent
#'  variance parameters. Elements set to numerical values denote fixed values
#'  constraints.
#' @param CORRFLAG Logical indicator. Set it to
#'  `FALSE` if the latent variables are assumed to be independent. Set it `TRUE`
#'  otherwise.
#' @param STDLV Logical indicator. Set it to
#'  `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.
#' @param CONSTRAINTS Linear equality constraints `A %*% theta = b` on the parameter vector, as
#'  `list(A = , b = )`, with one column of `A` per parameter in the order of [theta_names()]. Each
#'  constraint involves either thresholds, loadings, latent covariances (with `LATCOV = "direct"`),
#'  latent log standard deviations, regression coefficients or residual covariances.
#' @param REGMAT \eqn{q*q}-dimensional matrix of structural regressions among the latent variables, with
#'  `REGMAT[i, j]` the effect of latent variable `j` on latent variable `i`. Elements set to `NA` are
#'  free, numerical values are fixed. The latent covariances then refer to the disturbances.
#' @param RESCOV \eqn{p*p}-dimensional symmetric matrix of residual covariances of the items.
#'  Off-diagonal elements set to `NA` are free, numerical values are fixed. The diagonal is ignored.
#' @param LATCOV Parametrisation of the latent correlations: `"lkj"` (Cholesky factor with Fisher-z
#'  angles, always positive definite) or `"direct"` (the latent covariances themselves, which allows
#'  linear constraints on them).
#'
#' @return An object of class `latpl_model` with the identified loading constraints `CONSTRMAT`,
#'   the log standard deviations of the latent variables `CONSTRLOGSD` (`NA` when free),
#'   `CORRFLAG`, `STDLV`, `CONSTRAINTS`, `LATCOV`, `REGMAT`, `RESCOV`, and the numbers of items `p` and
#'   latent variables `q`.
#'
#' @export
set_latpl_model <- function(CONSTRMAT, CONSTRVAR = NULL, CORRFLAG = TRUE, STDLV = TRUE, CONSTRAINTS = NULL,
                           LATCOV = c("lkj", "direct"), REGMAT = NULL, RESCOV = NULL) {
  LATCOV <- match.arg(LATCOV)
  stopifnot(is.matrix(CONSTRMAT), !is.infinite(CONSTRMAT), CORRFLAG %in% c(0, 1), is.logical(STDLV))
  storage.mode(CONSTRMAT) <- "double"
  p <- nrow(CONSTRMAT)
  q <- ncol(CONSTRMAT)
  if (!STDLV) {
    for (h in seq_len(q)) {
      f <- which(is.na(CONSTRMAT[, h]))[1]
      if (!is.na(f) && all(CONSTRMAT[seq_len(f - 1), h] == 0)) CONSTRMAT[f, h] <- 1
    }
  }
  if (is.null(CONSTRVAR)) CONSTRVAR <- rep(if (STDLV) 1 else NA, q)
  stopifnot(length(CONSTRVAR) == q, CONSTRVAR[!is.na(CONSTRVAR)] > 0)
  if (STDLV) CONSTRVAR[is.na(CONSTRVAR)] <- 1
  if (!is.null(CONSTRAINTS)) {
    stopifnot(
      is.list(CONSTRAINTS), is.matrix(CONSTRAINTS$A), is.numeric(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$A)),
      is.numeric(CONSTRAINTS$b), length(CONSTRAINTS$b) == nrow(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$b))
    )
  }
  if (LATCOV == "direct" && !CORRFLAG) stop("LATCOV = \"direct\" requires CORRFLAG = TRUE.")
  if (is.null(REGMAT)) REGMAT <- matrix(0, 0, 0)
  if (is.null(RESCOV)) RESCOV <- matrix(0, 0, 0)
  storage.mode(REGMAT) <- storage.mode(RESCOV) <- "double"
  if (length(REGMAT)) {
    stopifnot(is.matrix(REGMAT), dim(REGMAT) == q, !is.infinite(REGMAT))
    if (any(is.na(diag(REGMAT)) | diag(REGMAT) != 0)) stop("The diagonal of REGMAT must be zero.")
    if (!anyNA(REGMAT) && abs(det(diag(q) - REGMAT)) < 1e-8) stop("diag(q) - REGMAT must be invertible.")
  }
  if (length(RESCOV)) {
    stopifnot(is.matrix(RESCOV), dim(RESCOV) == p, !is.infinite(RESCOV))
    diag(RESCOV) <- 0
    if (!identical(is.na(RESCOV), t(is.na(RESCOV))) || !isSymmetric(unname(RESCOV))) stop("RESCOV must be symmetric.")
  }
  out <- list(CONSTRMAT = CONSTRMAT, CONSTRLOGSD = log(sqrt(CONSTRVAR)), CORRFLAG = as.logical(CORRFLAG), STDLV = STDLV,
              CONSTRAINTS = CONSTRAINTS, LATCOV = LATCOV, p = p, q = q, REGMAT = REGMAT, RESCOV = RESCOV)
  class(out) <- "latpl_model"
  out
}

#' Starting values for pairwise likelihood estimation
#'
#' @param METHOD \tabular{ll}{ \code{'standard'} \tab Uses cold
#'  initialisation where loadings are initialised at `0.25` and the latent
#'  covariance matrix as an Identity matrix. \cr \tab \cr \code{'SA'} \tab
#'  Computes the starting point using a short chain of stochastic updates.
#'  Updates start from the `standard` initialisation point. \cr }
#'  With fewer than 15 items, `'SA'` falls back to `'standard'`.
#' @param THETA Initialising vector. If `NULL`, the `standard` initialisation point is used.
#' @param CONTROL Control object from [set_sa_init_control()], used when `METHOD = "SA"`.
#'
#' @return An object of class `latpl_init` to be passed to [fit_latpl()] via `INIT`, or to [cfa()]
#'   via `estimator.args$init`.
#'
#' @export
set_init <- function(METHOD = c("standard", "SA"), THETA = NULL, CONTROL = set_sa_init_control()) {
  METHOD <- match.arg(METHOD)
  stopifnot(
    is.null(THETA) || (is.numeric(THETA) && is.vector(THETA) && all(is.finite(THETA))),
    inherits(CONTROL, "latpl_sa_control")
  )
  out <- list(METHOD = METHOD, THETA = THETA, CONTROL = CONTROL)
  class(out) <- "latpl_init"
  out
}

#' Control options for the ucminf optimiser
#'
#' @param GRTOL Tolerance on the largest absolute gradient element.
#' @param XTOL Tolerance on the relative step length.
#' @param STEPMAX Initial maximal step length of the trust region.
#' @param MAXEVAL Maximum number of function evaluations.
#'
#' @return An object of class `latpl_ucminf_control` to be passed to [fit_latpl()] via `CONTROL`
#'   when `METHOD = "ucminf"`, or to [cfa()] via `control`.
#'
#' @export
set_ucminf_control <- function(
  GRTOL = 1e-6,
  XTOL = 1e-12,
  STEPMAX = 1,
  MAXEVAL = 500
) {
  stopifnot(
    is.numeric(GRTOL), length(GRTOL) == 1, GRTOL > 0,
    is.numeric(XTOL), length(XTOL) == 1, XTOL > 0,
    is.numeric(STEPMAX), length(STEPMAX) == 1, STEPMAX > 0,
    is.numeric(MAXEVAL), length(MAXEVAL) == 1, MAXEVAL >= 1
  )
  out <- list(GRTOL = GRTOL, XTOL = XTOL, STEPMAX = STEPMAX, MAXEVAL = as.integer(MAXEVAL))
  class(out) <- "latpl_ucminf_control"
  out
}

#' Control options for the pairwise stochastic approximation
#'
#' @param PAIRS_PER_ITERATION Number of item pairs drawn at each iteration.
#' @param STEP0 Initial step length.
#' @param STEP1 Scaling of the step length.
#' @param STEP2 Rate of the step length decay.
#' @param STEP3 Exponent of the step length decay.
#' @param BURNE Maximum number of burn-in epochs before trajectory averaging.
#' @param MAXE Maximum number of epochs.
#' @param EPOCHS_PER_CHECK Number of epochs between evaluations of the full negative pairwise
#'   log-likelihood.
#' @param CHECK_TOL If `TRUE`, stop the burn-in and the optimisation when the relative change in
#'   the full negative pairwise log-likelihood falls below `TOL_BURN` and `TOL_END`.
#' @param MAX_TOL_COUNTER Number of consecutive checks below tolerance required to stop.
#' @param TOL_BURN Tolerance ending the burn-in.
#' @param TOL_END Tolerance ending the optimisation.
#' @param SEED Seed of the pair sampling.
#'
#' @return An object of class `latpl_sa_control` to be passed to [fit_latpl()] via `CONTROL` when
#'   `METHOD = "SA"`, or to [cfa()] via `control`.
#'
#' @export
set_sa_control <- function(
  PAIRS_PER_ITERATION = 16,
  STEP0 = 1,
  STEP1 = 1,
  STEP2 = 1e-8,
  STEP3 = 0.75,
  BURNE = 20,
  MAXE = 30,
  EPOCHS_PER_CHECK = 5,
  CHECK_TOL = TRUE,
  MAX_TOL_COUNTER = 1,
  TOL_BURN = 5e-7,
  TOL_END = 1e-7,
  SEED = 123
) {
  stopifnot(
    is.numeric(PAIRS_PER_ITERATION), length(PAIRS_PER_ITERATION) == 1, PAIRS_PER_ITERATION >= 1,
    is.null(STEP0) || (is.numeric(STEP0) && length(STEP0) == 1 && STEP0 > 0),
    is.numeric(STEP1), length(STEP1) == 1, STEP1 > 0,
    is.numeric(STEP2), length(STEP2) == 1, STEP2 > 0,
    is.numeric(STEP3), length(STEP3) == 1, STEP3 > 0,
    is.numeric(BURNE), length(BURNE) == 1, BURNE >= 0,
    is.numeric(MAXE), length(MAXE) == 1, MAXE >= 0,
    BURNE < MAXE,
    is.numeric(EPOCHS_PER_CHECK), length(EPOCHS_PER_CHECK) == 1, EPOCHS_PER_CHECK >= 1,
    is.logical(CHECK_TOL), length(CHECK_TOL) == 1,
    is.numeric(MAX_TOL_COUNTER), length(MAX_TOL_COUNTER) == 1, MAX_TOL_COUNTER >= 1,
    is.numeric(TOL_BURN), length(TOL_BURN) == 1, TOL_BURN > 0,
    is.numeric(TOL_END), length(TOL_END) == 1, TOL_END > 0,
    is.numeric(SEED), length(SEED) == 1, SEED > 0
  )
  out <- list(
    PAIRS_PER_ITERATION = as.integer(PAIRS_PER_ITERATION),
    STEP0 = STEP0,
    STEP1 = STEP1,
    STEP2 = STEP2,
    STEP3 = STEP3,
    BURNE = as.integer(BURNE),
    MAXE = as.integer(MAXE),
    EPOCHS_PER_CHECK = as.integer(EPOCHS_PER_CHECK),
    CHECK_TOL = CHECK_TOL,
    MAX_TOL_COUNTER = as.integer(MAX_TOL_COUNTER),
    TOL_BURN = TOL_BURN,
    TOL_END = TOL_END,
    SEED = as.integer(SEED)
  )
  class(out) <- "latpl_sa_control"
  out
}

#' Control options for the stochastic approximation initialisation
#'
#' @inheritParams set_sa_control
#' @param STEP0 Initial step length. If `NULL`, it is set to \eqn{\sqrt{1/d}} at fit time, where
#'   \eqn{d} is the number of free parameters.
#'
#' @return An object of class `latpl_sa_control` to be passed to [fit_latpl()] via `INIT_CONTROL`,
#'   or to [cfa()] via `estimator.args$init_control`, when the initialisation method is `"SA"`.
#'
#' @export
set_sa_init_control <- function(
  PAIRS_PER_ITERATION = 16,
  STEP0 = NULL,
  STEP1 = 1,
  STEP2 = 1e-8,
  STEP3 = 0.75,
  BURNE = 3,
  MAXE = 5,
  EPOCHS_PER_CHECK = 1,
  CHECK_TOL = FALSE,
  MAX_TOL_COUNTER = 1,
  TOL_BURN = 5e-7,
  TOL_END = 1e-7,
  SEED = 123
) {
  set_sa_control(
    PAIRS_PER_ITERATION = PAIRS_PER_ITERATION,
    STEP0 = STEP0,
    STEP1 = STEP1,
    STEP2 = STEP2,
    STEP3 = STEP3,
    BURNE = BURNE,
    MAXE = MAXE,
    EPOCHS_PER_CHECK = EPOCHS_PER_CHECK,
    CHECK_TOL = CHECK_TOL,
    MAX_TOL_COUNTER = MAX_TOL_COUNTER,
    TOL_BURN = TOL_BURN,
    TOL_END = TOL_END,
    SEED = SEED
  )
}
