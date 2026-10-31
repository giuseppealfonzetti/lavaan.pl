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
#'  Elements set to `NA` are free, numerical values are fixed. The diagonal holds the residual
#'  variances with `PARAMETERIZATION = "theta"` (default 1), and is ignored otherwise.
#' @param LATCOV Parametrisation of the latent correlations: `"lkj"` (Cholesky factor with Fisher-z
#'  angles, always positive definite) or `"direct"` (the latent covariances themselves, which allows
#'  linear constraints on them).
#' @param PARAMETERIZATION `"delta"` (latent responses with unit variances in the reference metric,
#'  scaled by `SCALES`) or `"theta"` (residual variances in the diagonal of `RESCOV`), as in lavaan.
#' @param INTERCEPTS \eqn{p}-dimensional vector of item intercepts. Elements set to `NA` are free,
#'  numerical values are fixed. If `NULL`, all zero.
#' @param MEANS \eqn{q}-dimensional vector of latent means (intercepts, with `REGMAT`). Elements set
#'  to `NA` are free, numerical values are fixed. If `NULL`, all zero.
#' @param SCALES \eqn{p}-dimensional vector of positive scale factors of the latent responses, with
#'  `PARAMETERIZATION = "delta"`. Elements set to `NA` are free, numerical values are fixed. If
#'  `NULL`, all one.
#'
#' @return An object of class `latpl_model` with the identified loading constraints `CONSTRMAT`,
#'   the log standard deviations of the latent variables `CONSTRLOGSD` (`NA` when free),
#'   `CORRFLAG`, `STDLV`, `CONSTRAINTS`, `LATCOV`, `REGMAT`, `RESCOV`, `PARAMETERIZATION`, the
#'   intercepts `CONSTRINT` and latent means `CONSTRMEAN`, the log scale factors or log residual
#'   standard deviations `CONSTRLOGSCALE` (`NA` when free), and the numbers of items `p` and latent
#'   variables `q`.
#'
#' @export
set_latpl_model <- function(CONSTRMAT, CONSTRVAR = NULL, CORRFLAG = TRUE, STDLV = TRUE, CONSTRAINTS = NULL,
                           LATCOV = c("lkj", "direct"), REGMAT = NULL, RESCOV = NULL,
                           PARAMETERIZATION = c("delta", "theta"), INTERCEPTS = NULL, MEANS = NULL, SCALES = NULL) {
  LATCOV <- match.arg(LATCOV)
  PARAMETERIZATION <- match.arg(PARAMETERIZATION)
  stopifnot(is.matrix(CONSTRMAT), !is.infinite(CONSTRMAT), CORRFLAG %in% c(0, 1), is.logical(STDLV))
  storage.mode(CONSTRMAT) <- "double"
  p <- nrow(CONSTRMAT)
  q <- ncol(CONSTRMAT)
  if (!STDLV && is.null(CONSTRVAR)) {
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
  if (is.null(INTERCEPTS)) INTERCEPTS <- rep(0, p)
  if (is.null(MEANS)) MEANS <- rep(0, q)
  stopifnot(length(INTERCEPTS) == p, !is.infinite(INTERCEPTS), length(MEANS) == q, !is.infinite(MEANS))
  if (PARAMETERIZATION == "theta") {
    if (!is.null(SCALES)) stop("SCALES requires PARAMETERIZATION = \"delta\"; set the residual variances in the diagonal of RESCOV.")
    SCALES <- if (length(RESCOV)) sqrt(diag(RESCOV)) else rep(1, p)
  }
  if (is.null(SCALES)) SCALES <- rep(1, p)
  stopifnot(length(SCALES) == p, !is.infinite(SCALES))
  if (any(SCALES <= 0, na.rm = TRUE)) stop(if (PARAMETERIZATION == "theta") "Residual variances must be positive." else "SCALES must be positive.")
  if (length(RESCOV)) {
    stopifnot(is.matrix(RESCOV), dim(RESCOV) == p, !is.infinite(RESCOV))
    diag(RESCOV) <- 0
    if (!identical(is.na(RESCOV), t(is.na(RESCOV))) || !isSymmetric(unname(RESCOV))) stop("RESCOV must be symmetric.")
  }
  out <- list(CONSTRMAT = CONSTRMAT, CONSTRLOGSD = log(sqrt(CONSTRVAR)), CORRFLAG = as.logical(CORRFLAG), STDLV = STDLV,
              CONSTRAINTS = CONSTRAINTS, LATCOV = LATCOV, p = p, q = q, REGMAT = REGMAT, RESCOV = RESCOV,
              PARAMETERIZATION = PARAMETERIZATION, CONSTRINT = as.numeric(INTERCEPTS), CONSTRMEAN = as.numeric(MEANS),
              CONSTRLOGSCALE = log(as.numeric(SCALES)))
  class(out) <- "latpl_model"
  out
}

#' Starting values for pairwise likelihood estimation
#'
#' @param METHOD \tabular{ll}{ \code{'standard'} \tab Thresholds from the univariate margins,
#'  loadings and regression coefficients at `0.25`, uncorrelated latent variables with unit
#'  variances, no residual covariances. \cr \tab \cr \code{'SA'} \tab The end point of the Adam phase of the
#'  pairwise stochastic approximation ([set_sa_init_control()]), started from `THETA` or from the
#'  `standard` values. With [set_sa_control()], the Adam phase is already part of the optimiser and
#'  `'SA'` adds nothing. \cr }
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
#' The algorithm draws `PAIRS_PER_ITERATION` item pairs at each iteration, uniformly and
#' independently across iterations, and runs in three phases: Adam, a stochastic Newton burn-in
#' and stochastic gradient updates with a frozen gain and Polyak-Ruppert averaging. An epoch is
#' the number of iterations needed to draw as many pairs as there are in the model.
#'
#' @param PAIRS_PER_ITERATION Number of item pairs drawn at each iteration.
#' @param ALPHA Step length of the Adam phase, divided by the square root of the number of
#'   iterations per epoch.
#' @param ADAM_EPOCHS Maximum number of epochs of the Adam phase.
#' @param BURN_EPOCHS Maximum number of epochs of the stochastic Newton burn-in.
#' @param EPOCHS Number of epochs of the averaged phase.
#' @param TOL Improvement of the pairwise log-likelihood at the epoch average, per free parameter,
#'   below which, at two consecutive epochs, a burn-in phase ends.
#' @param DECAY Exponent of the step length decay in the averaged phase, in \eqn{(1/2, 1)}.
#' @param AVG_FROM Fraction of the averaged phase discarded before averaging.
#' @param GAIN Gain of the burn-in and averaged phases: `"hessian"` (inverse observed information),
#'   `"diagonal"` (Adam's diagonal preconditioner, no Newton burn-in), or `"auto"` (`"hessian"` up
#'   to 2000 free parameters).
#' @param SEED Seed of the pair sampling.
#'
#' @return An object of class `latpl_sa_control` to be passed to [fit_latpl()] via `CONTROL`, or to
#'   [cfa()] via `control`.
#'
#' @export
set_sa_control <- function(
  PAIRS_PER_ITERATION = 16,
  ALPHA = 0.2,
  ADAM_EPOCHS = 5,
  BURN_EPOCHS = 30,
  EPOCHS = 20,
  TOL = 0.05,
  DECAY = 0.75,
  AVG_FROM = 0.25,
  GAIN = c("auto", "hessian", "diagonal"),
  SEED = 123
) {
  GAIN <- match.arg(GAIN)
  num <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  stopifnot(
    num(PAIRS_PER_ITERATION), PAIRS_PER_ITERATION >= 1,
    num(ALPHA), ALPHA > 0,
    num(ADAM_EPOCHS), ADAM_EPOCHS >= 1,
    num(BURN_EPOCHS), BURN_EPOCHS >= 0,
    num(EPOCHS), EPOCHS >= 0,
    num(TOL), TOL > 0,
    num(DECAY), DECAY > 0.5, DECAY < 1,
    num(AVG_FROM), AVG_FROM >= 0, AVG_FROM < 1,
    num(SEED), SEED > 0
  )
  out <- list(PAIRS_PER_ITERATION = as.integer(PAIRS_PER_ITERATION), ALPHA = ALPHA, ADAM_EPOCHS = as.integer(ADAM_EPOCHS),
              BURN_EPOCHS = as.integer(BURN_EPOCHS), EPOCHS = as.integer(EPOCHS), TOL = TOL, DECAY = DECAY,
              AVG_FROM = AVG_FROM, GAIN = GAIN, SEED = as.integer(SEED))
  class(out) <- "latpl_sa_control"
  out
}

#' Control options for the stochastic approximation initialisation
#'
#' Only the Adam phase of [set_sa_control()] runs, and its end point is the starting value of
#' [set_ucminf_control()].
#'
#' @inheritParams set_sa_control
#'
#' @return An object of class `latpl_sa_control` to be passed to [set_init()] via `CONTROL`.
#'
#' @export
set_sa_init_control <- function(PAIRS_PER_ITERATION = 16, ALPHA = 0.2, ADAM_EPOCHS = 5, TOL = 0.05, SEED = 123) {
  set_sa_control(PAIRS_PER_ITERATION = PAIRS_PER_ITERATION, ALPHA = ALPHA, ADAM_EPOCHS = ADAM_EPOCHS,
                 EPOCHS = 0, TOL = TOL, SEED = SEED)
}
