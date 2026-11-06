#' Ordinal and numeric data for pairwise likelihood estimation
#'
#'@param DATA Data matrix of dimension \eqn{n*p}. Categories of ordinal items must be
#'  integers coded starting from zero.
#'@param MISSING Treatment of missing values (`NA`): `"listwise"` drops incomplete observations,
#'  `"pairwise"` uses the complete pairs and `"available.cases"` the available cases of
#'  Katsikatsou et al. (2022).
#'@param WEIGHTS Sampling weights, one per observation, normalised to sum to the number of
#'  observations. If `NULL`, unit weights.
#'@param WEIGHTS_TYPE Type of sampling weights, as in lavaan: `"design"` or `"frequency"`.
#'@param CLUSTER Cluster membership, one per observation, for cluster-robust standard errors. If
#'  `NULL`, observations are independent.
#'@param GROUP Group membership, one per observation, for multigroup models (see
#'  [set_latpl_groups()]). Groups are numbered in order of appearance. If `NULL`, a single group.
#'@param WEIGHTS_NORMALIZATION Normalisation of the sampling weights, as in lavaan: to the total
#'  number of observations (`"total"`) or to the number of observations of each group (`"group"`).
#'@param X Numeric matrix of covariates, one row per observation, on which the responses are
#'  conditioned. Observations with missing covariates are removed. If `NULL`, no covariates.
#'@param NUMERIC Logical vector marking the numeric items. If `NULL`, all items are ordinal.
#'
#' @return An object of class `latpl_data` with the data matrix `Y`, its dimensions `n` and `p`,
#'   the number of categories of each item `cat` (one for numeric items), the indicator of numeric
#'   items `numeric`, the pairwise frequency tables `freq`, the treatment of missing values `missing`, the unit weights `uweights` and
#'   univariate frequencies `ufreq` of the available-cases term, the normalised sampling weights
#'   `weights` and their factors on the scores `jweights`, the cluster indices `cluster`, the group
#'   indices `group` (from zero) and their number `ngroups`, the covariates `X`, the stratum of each
#'   observation within its group `stratum` (from zero) and the distinct covariate rows of each group
#'   `strata`, whether each group uses frequency tables `tables` or its individual observations, and
#'   the time spent computing the frequencies `freqTime`. With several groups or strata, `freq` and
#'   `ufreq` stack the tables of each group and, within it, of each stratum.
#'
#' @keywords internal
set_latpl_data <- function(DATA, MISSING = c("listwise", "pairwise", "available.cases"),
                          WEIGHTS = NULL, WEIGHTS_TYPE = c("design", "frequency"), CLUSTER = NULL, GROUP = NULL,
                          WEIGHTS_NORMALIZATION = c("total", "group"), X = NULL, NUMERIC = NULL) {
  MISSING <- match.arg(MISSING)
  WEIGHTS_TYPE <- match.arg(WEIGHTS_TYPE)
  WEIGHTS_NORMALIZATION <- match.arg(WEIGHTS_NORMALIZATION)
  stopifnot(is.matrix(DATA), is.numeric(DATA), all(is.finite(DATA) | is.na(DATA)))
  numeric <- if (is.null(NUMERIC)) rep(FALSE, ncol(DATA)) else NUMERIC
  stopifnot(is.logical(numeric), length(numeric) == ncol(DATA), !anyNA(numeric), all(numeric) || min(DATA[, !numeric], na.rm = TRUE) == 0)
  if (MISSING == "available.cases" && any(numeric)) stop("MISSING = \"available.cases\" is not available with numeric items.")
  Y <- DATA
  storage.mode(Y) <- "double"
  Y[, !numeric] <- round(Y[, !numeric])
  stopifnot(
    is.null(WEIGHTS) || (is.numeric(WEIGHTS) && length(WEIGHTS) == nrow(Y) && all(is.finite(WEIGHTS)) && all(WEIGHTS >= 0)),
    is.null(CLUSTER) || (length(CLUSTER) == nrow(Y) && !anyNA(CLUSTER)),
    is.null(GROUP) || (length(GROUP) == nrow(Y) && !anyNA(GROUP)),
    is.null(X) || (is.matrix(X) && is.numeric(X) && nrow(X) == nrow(Y) && !any(is.infinite(X)))
  )
  group <- if (is.null(GROUP)) rep(0, nrow(Y)) else match(GROUP, unique(GROUP)) - 1
  if (is.null(X)) X <- matrix(0, nrow(Y), 0)
  storage.mode(X) <- "double"
  G <- max(group) + 1
  keep_x <- stats::complete.cases(X)
  keep_y <- MISSING != "listwise" | stats::complete.cases(Y)
  if (!all(keep_x)) message(sum(!keep_x), " observations with missing covariates removed.")
  if (!all(keep_y)) message(sum(!keep_y), " incomplete observations removed (MISSING = \"listwise\").")
  keep <- keep_x & keep_y
  if (!all(keep)) {
    Y <- Y[keep, , drop = FALSE]
    X <- X[keep, , drop = FALSE]
    if (!is.null(WEIGHTS)) WEIGHTS <- WEIGHTS[keep]
    if (!is.null(CLUSTER)) CLUSTER <- CLUSTER[keep]
    group <- group[keep]
  }
  stopifnot(nrow(Y) > 0, colSums(!is.na(Y)) > 0)
  cat <- apply(Y, 2, max, na.rm = TRUE) + 1
  cat[numeric] <- 1
  stopifnot(cat[!numeric] > 1)
  weights <- jweights <- cluster <- numeric(0)
  if (!is.null(WEIGHTS)) {
    weights <- if (WEIGHTS_NORMALIZATION == "group") {
      WEIGHTS / stats::ave(WEIGHTS, group, FUN = sum) * stats::ave(WEIGHTS, group, FUN = length)
    } else {
      WEIGHTS / sum(WEIGHTS) * nrow(Y)
    }
    jweights <- if (WEIGHTS_TYPE == "design") weights else sqrt(weights)
  }
  if (!is.null(CLUSTER)) {
    cluster <- as.numeric(match(paste(group, CLUSTER), unique(paste(group, CLUSTER))) - 1)
    if (any(tapply(cluster, group, function(x) length(unique(x))) < 2)) stop("Every group needs at least two clusters.")
  }
  freq_start_time <- Sys.time()
  uweights <- if (MISSING == "available.cases") rowSums(is.na(Y)) else rep(0, nrow(Y))
  uw <- if (length(weights)) uweights * weights else uweights
  st <- strata_index(X, group, Y, cat, if (length(weights)) weights else numeric(0))
  ufreq <- numeric(0)
  if (any(uw > 0)) {
    uoff <- cumsum(c(0, cat[-length(cat)]))
    sizes <- vapply(st$strata, nrow, numeric(1))
    ufreq <- numeric(sum(sizes) * sum(cat))
    for (j in seq_len(ncol(Y))) {
      seen <- which(!is.na(Y[, j]) & uw > 0)
      cells <- (cumsum(c(0, sizes[-G]))[group[seen] + 1] + st$stratum[seen]) * sum(cat) + uoff[j] + Y[seen, j] + 1
      added <- rowsum(uw[seen], cells)
      at <- as.numeric(rownames(added))
      ufreq[at] <- ufreq[at] + added[, 1]
    }
  }
  out <- list(
    Y = Y,
    n = nrow(Y),
    p = ncol(Y),
    cat = cat,
    numeric = numeric,
    freq = unlist(st$freq),
    missing = MISSING,
    uweights = uweights,
    ufreq = ufreq,
    weights = weights,
    jweights = jweights,
    cluster = cluster,
    group = group,
    ngroups = G,
    X = X,
    stratum = st$stratum,
    strata = st$strata,
    tables = st$tables,
    freqTime = as.numeric(difftime(Sys.time(), freq_start_time, units = "secs"))
  )
  class(out) <- "latpl_data"
  out
}

#' Distinct covariate rows and their frequency tables
#'
#' @param X Matrix of covariates.
#' @param GROUP Group of each row, from zero.
#' @param Y Data matrix.
#' @param CAT Number of categories of each item.
#' @param WEIGHTS Unit weights, or empty for unit weights.
#' @return A list with, for each group, the distinct rows of `X` in order of appearance `strata`, whether frequency
#'   tables are used `tables` (see [use_tables()]) and the tables `freq` (`NULL` when not used), and the index `stratum`
#'   (from zero) of each row among the strata of its group.
#' @seealso [set_latpl_data()]
#' @keywords internal
strata_index <- function(X, GROUP, Y, CAT, WEIGHTS) {
  G <- max(GROUP) + 1
  out <- list(strata = vector("list", G), tables = logical(G), freq = vector("list", G), stratum = numeric(nrow(X)))
  for (g in seq_len(G)) {
    rows <- GROUP == g - 1
    key <- do.call(paste, as.data.frame(X[rows, , drop = FALSE]))
    first <- if (ncol(X)) !duplicated(key) else which.max(rows) == which(rows)
    out$strata[[g]] <- X[rows, , drop = FALSE][first, , drop = FALSE]
    out$stratum[rows] <- if (ncol(X)) match(key, key[first]) - 1 else 0
    out$tables[g] <- use_tables(nrow(out$strata[[g]]), Y[rows, , drop = FALSE], CAT, out$stratum[rows])
    if (out$tables[g]) {
      out$freq[[g]] <- cpp_pairs_freq(Y[rows, , drop = FALSE], CAT, if (length(WEIGHTS)) WEIGHTS[rows] else numeric(0), out$stratum[rows])
    }
  }
  out
}

#' Whether frequency tables are cheaper than individual observations
#'
#' Compares the bivariate normal probabilities evaluated by each route over all pairs of ordinal items. A table
#' evaluates \eqn{(c_k - 1)(c_l - 1)} of them for each pair and stratum. An observation evaluates only the corners of its
#' cell: \eqn{a_k a_l}, where \eqn{a} is 1 for an extreme category, 2 for a middle one and 0 for a missing response,
#' and identical observations count once.
#'
#' @param STRATA Number of distinct covariate rows.
#' @param Y Data matrix.
#' @param CAT Number of categories of each item.
#' @param STRATUM Stratum of each row of `Y`.
#' @return `TRUE` if the tables should be used. They must also fit in memory.
#' @seealso [set_latpl_data()], [engine()]
#' @keywords internal
use_tables <- function(STRATA, Y, CAT, STRATUM) {
  ordinal <- CAT > 1
  if (STRATA == 1 || sum(ordinal) < 2) return(TRUE)
  pairs_of <- function(a, b) (a^2 - b) / 2
  cells <- pairs_of(sum(CAT[ordinal]), sum(CAT[ordinal]^2))
  if (STRATA * cells > 5e7) return(FALSE)
  rows <- if (any(!ordinal)) Y else unique(cbind(STRATUM, Y))[, -1, drop = FALSE]
  categories <- rows[, ordinal, drop = FALSE]
  a <- 2 - (categories == 0) - sweep(categories, 2, CAT[ordinal] - 1, "==")
  a[is.na(a)] <- 0
  STRATA * pairs_of(sum(CAT[ordinal] - 1), sum((CAT[ordinal] - 1)^2)) < sum(pairs_of(rowSums(a), rowSums(a^2)))
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
#' @param GAMMA \eqn{q*x}-dimensional matrix of regressions of the latent variables on the `x` covariates of
#'  [set_latpl_data()]. Elements set to `NA` are free, numerical values are fixed. If `NULL`, all zero.
#' @param KAPPA \eqn{p*x}-dimensional matrix of direct effects of the covariates on the latent responses.
#'  Elements set to `NA` are free, numerical values are fixed. If `NULL`, all zero.
#' @param NUMERIC Logical vector marking the numeric items, as in [set_latpl_data()]. Unless `INTERCEPTS`
#'  and `SCALES` say otherwise, their intercepts and residual standard deviations are free. For them,
#'  `SCALES` holds the residual standard deviations under both parameterizations.
#'
#' @return An object of class `latpl_model` with the identified loading constraints `CONSTRMAT`,
#'   the log standard deviations of the latent variables `CONSTRLOGSD` (`NA` when free),
#'   `CORRFLAG`, `STDLV`, `CONSTRAINTS`, `LATCOV`, `REGMAT`, `RESCOV`, `PARAMETERIZATION`, the
#'   intercepts `CONSTRINT` and latent means `CONSTRMEAN`, the log scale factors or log residual
#'   standard deviations `CONSTRLOGSCALE` (`NA` when free), the covariate effects `GAMMA` and `KAPPA`,
#'   the indicator of numeric items `NUMERIC`, and the numbers of items `p` and latent variables `q`.
#'
#' @keywords internal
set_latpl_model <- function(CONSTRMAT, CONSTRVAR = NULL, CORRFLAG = TRUE, STDLV = TRUE, CONSTRAINTS = NULL,
                           LATCOV = c("lkj", "direct"), REGMAT = NULL, RESCOV = NULL,
                           PARAMETERIZATION = c("delta", "theta"), INTERCEPTS = NULL, MEANS = NULL, SCALES = NULL,
                           GAMMA = NULL, KAPPA = NULL, NUMERIC = NULL) {
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
  if (is.null(NUMERIC)) NUMERIC <- rep(FALSE, p)
  stopifnot(is.logical(NUMERIC), length(NUMERIC) == p, !anyNA(NUMERIC))
  if (is.null(INTERCEPTS)) INTERCEPTS <- ifelse(NUMERIC, NA, 0)
  if (is.null(MEANS)) MEANS <- rep(0, q)
  stopifnot(length(INTERCEPTS) == p, !is.infinite(INTERCEPTS), length(MEANS) == q, !is.infinite(MEANS))
  if (PARAMETERIZATION == "theta") {
    if (!is.null(SCALES)) stop("SCALES requires PARAMETERIZATION = \"delta\"; set the residual variances in the diagonal of RESCOV.")
    SCALES <- if (length(RESCOV)) sqrt(diag(RESCOV)) else ifelse(NUMERIC, NA, 1)
  }
  if (is.null(SCALES)) SCALES <- ifelse(NUMERIC, NA, 1)
  stopifnot(length(SCALES) == p, !is.infinite(SCALES))
  if (any(SCALES <= 0, na.rm = TRUE)) stop(if (PARAMETERIZATION == "theta") "Residual variances must be positive." else "SCALES must be positive.")
  if (length(RESCOV)) {
    stopifnot(is.matrix(RESCOV), dim(RESCOV) == p, !is.infinite(RESCOV))
    diag(RESCOV) <- 0
    if (!identical(is.na(RESCOV), t(is.na(RESCOV))) || !isSymmetric(unname(RESCOV))) stop("RESCOV must be symmetric.")
  }
  nx <- max(NCOL(GAMMA) * !is.null(GAMMA), NCOL(KAPPA) * !is.null(KAPPA))
  if (is.null(GAMMA)) GAMMA <- matrix(0, q, nx)
  if (is.null(KAPPA)) KAPPA <- matrix(0, p, nx)
  stopifnot(is.matrix(GAMMA), is.matrix(KAPPA), dim(GAMMA) == c(q, nx), dim(KAPPA) == c(p, nx), !is.infinite(GAMMA), !is.infinite(KAPPA))
  storage.mode(GAMMA) <- storage.mode(KAPPA) <- "double"
  out <- list(CONSTRMAT = CONSTRMAT, CONSTRLOGSD = log(sqrt(CONSTRVAR)), CORRFLAG = as.logical(CORRFLAG), STDLV = STDLV,
              CONSTRAINTS = CONSTRAINTS, LATCOV = LATCOV, p = p, q = q, REGMAT = REGMAT, RESCOV = RESCOV,
              PARAMETERIZATION = PARAMETERIZATION, CONSTRINT = as.numeric(INTERCEPTS), CONSTRMEAN = as.numeric(MEANS),
              CONSTRLOGSCALE = log(as.numeric(SCALES)), GAMMA = GAMMA, KAPPA = KAPPA, NUMERIC = NUMERIC)
  class(out) <- "latpl_model"
  out
}

#' Multigroup factor model for pairwise likelihood estimation
#'
#' @param MODELS List of model objects from [set_latpl_model()], one per group in the order of the
#'   groups of [set_latpl_data()], with the same number of items.
#' @param CONSTRAINTS Linear equality constraints `A %*% theta = b` across groups, as
#'  `list(A = , b = )`, with one column of `A` per parameter in the order of [theta_names()]: the
#'  parameters of the first group, then those of the second, and so on. Each model's own
#'  `CONSTRAINTS` apply within its group.
#'
#' @return An object of class `latpl_groups` with `MODELS`, `CONSTRAINTS`, the number of items `p`
#'   and the number of groups `ngroups`, to be passed to [fit_latpl()] via `MODEL`.
#'
#' @keywords internal
set_latpl_groups <- function(MODELS, CONSTRAINTS = NULL) {
  stopifnot(
    is.list(MODELS), length(MODELS) > 0, all(vapply(MODELS, inherits, logical(1), "latpl_model")),
    length(unique(vapply(MODELS, function(x) x$p, numeric(1)))) == 1
  )
  if (!is.null(CONSTRAINTS)) {
    stopifnot(
      is.list(CONSTRAINTS), is.matrix(CONSTRAINTS$A), is.numeric(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$A)),
      is.numeric(CONSTRAINTS$b), length(CONSTRAINTS$b) == nrow(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$b))
    )
  }
  out <- list(MODELS = MODELS, CONSTRAINTS = CONSTRAINTS, p = MODELS[[1]]$p, ngroups = length(MODELS))
  class(out) <- "latpl_groups"
  out
}

#' Control options of the quasi-Newton method
#'
#' The `QN` method: the exact pairwise likelihood estimate, by the BFGS optimiser of the \code{ucminf} package.
#'
#' @param GRTOL Tolerance on the largest absolute gradient element.
#' @param XTOL Tolerance on the relative step length.
#' @param STEPMAX Initial maximal step length of the trust region.
#' @param MAXEVAL Maximum number of function evaluations.
#'
#' @return An object of class `latpl_qn_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_qn_control <- function(
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
  class(out) <- "latpl_qn_control"
  out
}

#' Control options of the pairwise stochastic approximation
#'
#' The `SA` method. It draws item pairs at each iteration, uniformly and independently across iterations, and runs in
#' three phases: Adam, preconditioning, and stochastic gradient updates with a frozen gain and Polyak-Ruppert averaging. An
#' epoch is the number of iterations needed to draw as many pairs as there are in the model, one pass over the pairs.
#'
#' @param PAIRS_PER_ITERATION Number of item pairs drawn at each iteration.
#' @param ALPHA Step length of the Adam phase, divided by the square root of the number of
#'   iterations per epoch.
#' @param ADAM_EPOCHS Maximum number of epochs of the Adam phase. If `NULL`, 5, or as many as make 150 iterations on small
#'   models, whose epochs are short.
#' @param PRECONDITION_EPOCHS Number of epochs of the preconditioning phase, which draws each pair once per epoch,
#'   recomputes the observed information as the gain after each epoch, and takes stochastic Newton steps with it. The
#'   gain is then frozen, unless the observed information is not positive definite. With `0`, the gain is computed at
#'   the end of the Adam phase.
#' @param MAX_PRECONDITION_EPOCHS Maximum number of epochs of the preconditioning phase while the observed information
#'   is not positive definite.
#' @param AVERAGE_EPOCHS Number of epochs of the averaged phase, at least 1.
#' @param TOL Absolute change of the pairwise log-likelihood at the epoch average, per free parameter, below which, at
#'   two consecutive epochs, the Adam or the preconditioning phase ends before its maximum.
#' @param DECAY Exponent of the step length decay in the averaged phase, in \eqn{(1/2, 1)}.
#' @param AVG_FROM Fraction of the averaged phase discarded before averaging.
#' @param GAIN Gain of the preconditioning and averaged phases: `"hessian"` (inverse observed information),
#'   `"diagonal"` (Adam's diagonal preconditioner, no preconditioning phase), or `"auto"` (`"hessian"` up to 2000 free
#'   parameters).
#' @param SEED Seed of the pair sampling.
#'
#' @return An object of class `latpl_sa_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_sa_control <- function(
  PAIRS_PER_ITERATION = 16,
  ALPHA = 0.2,
  ADAM_EPOCHS = NULL,
  PRECONDITION_EPOCHS = 3,
  MAX_PRECONDITION_EPOCHS = 20,
  AVERAGE_EPOCHS = 50,
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
    is.null(ADAM_EPOCHS) || (num(ADAM_EPOCHS) && ADAM_EPOCHS >= 1),
    num(PRECONDITION_EPOCHS), PRECONDITION_EPOCHS >= 0,
    num(MAX_PRECONDITION_EPOCHS), MAX_PRECONDITION_EPOCHS >= 0,
    num(AVERAGE_EPOCHS), AVERAGE_EPOCHS >= 0,
    num(TOL), TOL > 0,
    num(DECAY), DECAY > 0.5, DECAY < 1,
    num(AVG_FROM), AVG_FROM >= 0, AVG_FROM < 1,
    num(SEED), SEED > 0
  )
  out <- list(PAIRS_PER_ITERATION = as.integer(PAIRS_PER_ITERATION), ALPHA = ALPHA,
              ADAM_EPOCHS = if (!is.null(ADAM_EPOCHS)) as.integer(ADAM_EPOCHS), PRECONDITION_EPOCHS = as.integer(PRECONDITION_EPOCHS),
              MAX_PRECONDITION_EPOCHS = as.integer(MAX_PRECONDITION_EPOCHS), AVERAGE_EPOCHS = as.integer(AVERAGE_EPOCHS),
              TOL = TOL, DECAY = DECAY, AVG_FROM = AVG_FROM, GAIN = GAIN,
              SEED = as.integer(SEED))
  class(out) <- "latpl_sa_control"
  out
}

#' Control options of quasi-Newton from the burn-in of stochastic approximation
#'
#' The `SA-QN` method: the Adam and preconditioning phases of [set_sa_control()], without the averaged phase, followed by
#' [set_qn_control()] from their end point, with the inverse of the observed information as the initial inverse Hessian.
#' Above 2000 free parameters, it starts from the end of the Adam phase, without curvature.
#'
#' @inheritParams set_sa_control
#' @inheritParams set_qn_control
#' @param NSTARTS Number of full runs from the same starting values, with seeds `SEED`, `SEED + 1` and so on. The highest
#'   pairwise log-likelihood is kept, with a warning if the runs disagree. Each run adds its full cost.
#'
#' @return An object of class `latpl_saqn_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_saqn_control <- function(
  PAIRS_PER_ITERATION = 16,
  ALPHA = 0.2,
  ADAM_EPOCHS = NULL,
  PRECONDITION_EPOCHS = 3,
  MAX_PRECONDITION_EPOCHS = 20,
  TOL = 0.05,
  SEED = 123,
  NSTARTS = 1,
  GRTOL = 1e-6,
  XTOL = 1e-12,
  STEPMAX = 1,
  MAXEVAL = 500
) {
  stopifnot(is.numeric(NSTARTS), length(NSTARTS) == 1, is.finite(NSTARTS), NSTARTS >= 1)
  out <- c(
    set_sa_control(PAIRS_PER_ITERATION = PAIRS_PER_ITERATION, ALPHA = ALPHA, ADAM_EPOCHS = ADAM_EPOCHS,
                   PRECONDITION_EPOCHS = PRECONDITION_EPOCHS, MAX_PRECONDITION_EPOCHS = MAX_PRECONDITION_EPOCHS,
                   AVERAGE_EPOCHS = 0, TOL = TOL, SEED = SEED),
    set_qn_control(GRTOL = GRTOL, XTOL = XTOL, STEPMAX = STEPMAX, MAXEVAL = MAXEVAL),
    list(NSTARTS = as.integer(NSTARTS))
  )
  class(out) <- c("latpl_saqn_control", "latpl_qn_control", "latpl_sa_control")
  out
}

#' Method of a control object
#'
#' @param CONTROL Control object from [set_saqn_control()], [set_sa_control()] or [set_qn_control()].
#' @return `"SA-QN"`, `"SA"` or `"QN"`.
#' @seealso [fit_latpl()], [pmlcfa()]
#' @keywords internal
control_method <- function(CONTROL) {
  if (inherits(CONTROL, "latpl_saqn_control")) return("SA-QN")
  if (inherits(CONTROL, "latpl_sa_control")) "SA" else "QN"
}

