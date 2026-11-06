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
#'  [lavaanpl_groups()]). Groups are numbered in order of appearance. If `NULL`, a single group.
#'@param WEIGHTS_NORMALIZATION Normalisation of the sampling weights, as in lavaan: to the total
#'  number of observations (`"total"`) or to the number of observations of each group (`"group"`).
#'@param X Numeric matrix of covariates, one row per observation, on which the responses are
#'  conditioned. Observations with missing covariates are removed. If `NULL`, no covariates.
#'@param NUMERIC Logical vector marking the numeric items. If `NULL`, all items are ordinal.
#'
#' @return An object of class `lavaanpl_data` with the data matrix `Y`, its dimensions `n` and `p`, the number of item pairs
#'   over all groups `pairs`, the number of categories of each item `cat` (one for numeric items), the indicator of numeric
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
lavaanpl_data <- function(DATA, MISSING = c("listwise", "pairwise", "available.cases"),
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
    pairs = G * ncol(Y) * (ncol(Y) - 1) / 2,
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
  class(out) <- "lavaanpl_data"
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
#' @seealso [lavaanpl_data()]
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
#' @seealso [lavaanpl_data()], [engine()]
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
#' Each matrix argument is the model matrix of lavaan with the same name: elements set to `NA` are free, numerical
#' values are fixed.
#'
#' @param LAMBDA \eqn{p*q}-dimensional matrix of loadings. Latent variables with a free variance need a fixed nonzero
#'  loading.
#' @param PSI \eqn{q*q}-dimensional symmetric matrix of latent (disturbance) covariances. If `NULL`, unit variances and
#'  free covariances. Covariances all fixed to zero make the latent variables independent.
#' @param BETA \eqn{q*q}-dimensional matrix of structural regressions among the latent variables, with `BETA[i, j]` the
#'  effect of latent variable `j` on latent variable `i`. If `NULL`, none.
#' @param THETA \eqn{p*p}-dimensional symmetric matrix of residual covariances of the items. The diagonal holds the
#'  residual variances of the numeric items and, with `PARAMETERIZATION = "theta"`, of all items; it is ignored for the
#'  other items. If `NULL`, no residual covariances, free residual variances of the numeric items and unit residual
#'  variances of the others.
#' @param NU \eqn{p}-dimensional vector of item intercepts. If `NULL`, free for the numeric items and zero otherwise.
#' @param ALPHA \eqn{q}-dimensional vector of latent means (intercepts, with `BETA`). If `NULL`, all zero.
#' @param DELTA \eqn{p}-dimensional vector of positive scale factors of the latent responses, with
#'  `PARAMETERIZATION = "delta"`; ignored for the numeric items. If `NULL`, all one.
#' @param TAU Vector of thresholds, item after item, with the length set by the data. If `NULL`, all free.
#' @param GAMMA \eqn{q*x}-dimensional matrix of regressions of the latent variables on the `x` covariates of
#'  [lavaanpl_data()]. If `NULL`, none.
#' @param KAPPA \eqn{p*x}-dimensional matrix of direct effects of the covariates on the latent responses. If `NULL`,
#'  none.
#' @param CONSTRAINTS Linear equality constraints `A %*% parameter_vector = b` among free parameters, as
#'  `list(A = , b = )`, with one column of `A` per parameter in the order of [parameter_names()]. Each constraint
#'  involves a single block of [parameter_blocks()].
#' @param LATCOV Parametrisation of the latent correlations: `"lkj"` (Cholesky factor with Fisher-z angles, always
#'  positive definite, all covariances free) or `"direct"` (the latent covariances themselves, which can be fixed or
#'  constrained).
#' @param PARAMETERIZATION `"delta"` (latent responses with unit variances in the reference metric, scaled by `DELTA`)
#'  or `"theta"` (residual variances in the diagonal of `THETA`), as in lavaan.
#' @param NUMERIC Logical vector marking the numeric items, as in [lavaanpl_data()].
#'
#' @return An object of class `lavaanpl_model` with the matrices `LAMBDA`, `PSI`, `BETA`, `THETA` (with zero diagonal
#'   entries where they are ignored), `NU`, `ALPHA`, `DELTA` (one where it is ignored), `TAU`, `GAMMA` and `KAPPA`,
#'   `CONSTRAINTS`, `LATCOV` (`"none"` for independent latent variables), `PARAMETERIZATION`, `NUMERIC`, and the
#'   numbers of items `p` and latent variables `q`.
#'
#' @keywords internal
lavaanpl_model <- function(LAMBDA, PSI = NULL, BETA = NULL, THETA = NULL, NU = NULL, ALPHA = NULL, DELTA = NULL,
                            TAU = NULL, GAMMA = NULL, KAPPA = NULL, CONSTRAINTS = NULL, LATCOV = c("lkj", "direct"),
                            PARAMETERIZATION = c("delta", "theta"), NUMERIC = NULL) {
  LATCOV <- match.arg(LATCOV)
  PARAMETERIZATION <- match.arg(PARAMETERIZATION)
  stopifnot(is.matrix(LAMBDA))
  p <- nrow(LAMBDA)
  q <- ncol(LAMBDA)
  if (is.null(NUMERIC)) NUMERIC <- rep(FALSE, p)
  stopifnot(is.logical(NUMERIC), length(NUMERIC) == p, !anyNA(NUMERIC))
  rsd <- PARAMETERIZATION == "theta" | NUMERIC
  if (!is.null(DELTA) && PARAMETERIZATION == "theta") {
    stop("DELTA requires PARAMETERIZATION = \"delta\"; set the residual variances in the diagonal of THETA.")
  }
  if (is.null(PSI)) {
    PSI <- matrix(NA, q, q)
    diag(PSI) <- 1
  }
  if (is.null(BETA)) BETA <- matrix(0, q, q)
  if (is.null(THETA)) THETA <- diag(ifelse(NUMERIC, NA, 1), p)
  if (is.null(NU)) NU <- ifelse(NUMERIC, NA, 0)
  if (is.null(ALPHA)) ALPHA <- rep(0, q)
  if (is.null(DELTA)) DELTA <- rep(1, p)
  nx <- max(NCOL(GAMMA) * !is.null(GAMMA), NCOL(KAPPA) * !is.null(KAPPA))
  if (is.null(GAMMA)) GAMMA <- matrix(0, q, nx)
  if (is.null(KAPPA)) KAPPA <- matrix(0, p, nx)
  stopifnot(
    is.matrix(PSI), dim(PSI) == q, is.matrix(BETA), dim(BETA) == q, is.matrix(THETA), dim(THETA) == p,
    is.matrix(GAMMA), dim(GAMMA) == c(q, nx), is.matrix(KAPPA), dim(KAPPA) == c(p, nx),
    length(NU) == p, length(ALPHA) == q, length(DELTA) == p,
    !is.infinite(c(LAMBDA, PSI, BETA, THETA, NU, ALPHA, DELTA, TAU, GAMMA, KAPPA))
  )
  for (M in list(PSI = PSI, THETA = THETA)) {
    if (!identical(is.na(M), t(is.na(M))) || !isSymmetric(unname(M))) stop("PSI and THETA must be symmetric.")
  }
  if (any(diag(PSI) <= 0, na.rm = TRUE)) stop("Fixed latent variances must be positive.")
  if (all(PSI[upper.tri(PSI)] %in% 0)) {
    LATCOV <- "none"
  } else if (LATCOV == "lkj" && !all(is.na(PSI[upper.tri(PSI)]))) {
    stop("Fixed latent covariances require LATCOV = \"direct\".")
  }
  if (any(is.na(diag(BETA)) | diag(BETA) != 0)) stop("The diagonal of BETA must be zero.")
  if (!anyNA(BETA) && abs(det(diag(q) - BETA)) < 1e-8) stop("diag(q) - BETA must be invertible.")
  if (anyNA(diag(THETA)[!rsd])) stop("Free residual variances require PARAMETERIZATION = \"theta\".")
  if (any(diag(THETA)[rsd] <= 0, na.rm = TRUE)) stop("Residual variances must be positive.")
  if (any(DELTA[!rsd] <= 0, na.rm = TRUE)) stop("DELTA must be positive.")
  diag(THETA)[!rsd] <- 0
  DELTA[rsd] <- 1
  if (!is.null(CONSTRAINTS)) {
    stopifnot(
      is.list(CONSTRAINTS), is.matrix(CONSTRAINTS$A), is.numeric(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$A)),
      is.numeric(CONSTRAINTS$b), length(CONSTRAINTS$b) == nrow(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$b))
    )
  }
  storage.mode(LAMBDA) <- storage.mode(PSI) <- storage.mode(BETA) <- storage.mode(THETA) <- "double"
  storage.mode(GAMMA) <- storage.mode(KAPPA) <- "double"
  out <- list(LAMBDA = LAMBDA, PSI = PSI, BETA = BETA, THETA = THETA, NU = as.numeric(NU), ALPHA = as.numeric(ALPHA),
              DELTA = as.numeric(DELTA), TAU = if (!is.null(TAU)) as.numeric(TAU), GAMMA = GAMMA, KAPPA = KAPPA,
              CONSTRAINTS = CONSTRAINTS, LATCOV = LATCOV, PARAMETERIZATION = PARAMETERIZATION, NUMERIC = NUMERIC,
              p = p, q = q)
  class(out) <- "lavaanpl_model"
  out
}

#' Multigroup factor model for pairwise likelihood estimation
#'
#' @param MODELS List of model objects from [lavaanpl_model()], one per group in the order of the
#'   groups of [lavaanpl_data()], with the same number of items.
#' @param CONSTRAINTS Linear equality constraints `A %*% parameter_vector = b` across groups, as
#'  `list(A = , b = )`, with one column of `A` per parameter in the order of [parameter_names()]: the
#'  parameters of the first group, then those of the second, and so on. Each model's own
#'  `CONSTRAINTS` apply within its group.
#'
#' @return An object of class `lavaanpl_groups` with `MODELS`, `CONSTRAINTS`, the number of items `p`
#'   and the number of groups `ngroups`, to be passed to [fit_lavaanpl()] via `MODEL`.
#'
#' @keywords internal
lavaanpl_groups <- function(MODELS, CONSTRAINTS = NULL) {
  stopifnot(
    is.list(MODELS), length(MODELS) > 0, all(vapply(MODELS, inherits, logical(1), "lavaanpl_model")),
    length(unique(vapply(MODELS, function(x) x$p, numeric(1)))) == 1
  )
  if (!is.null(CONSTRAINTS)) {
    stopifnot(
      is.list(CONSTRAINTS), is.matrix(CONSTRAINTS$A), is.numeric(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$A)),
      is.numeric(CONSTRAINTS$b), length(CONSTRAINTS$b) == nrow(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$b))
    )
  }
  out <- list(MODELS = MODELS, CONSTRAINTS = CONSTRAINTS, p = MODELS[[1]]$p, ngroups = length(MODELS))
  class(out) <- "lavaanpl_groups"
  out
}

#' Control options of the quasi-Newton method
#'
#' The `QN` method: the exact pairwise likelihood estimate, by the BFGS optimiser of the \code{ucminf} package.
#'
#' @param grtol Tolerance on the largest absolute gradient element.
#' @param xtol Tolerance on the relative step length.
#' @param stepmax Initial maximal step length of the trust region.
#' @param maxeval Maximum number of function evaluations.
#'
#' @return An object of class `lavaanpl_qn_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_qn_control <- function(
  grtol = 1e-6,
  xtol = 1e-12,
  stepmax = 1,
  maxeval = 500
) {
  stopifnot(
    is.numeric(grtol), length(grtol) == 1, grtol > 0,
    is.numeric(xtol), length(xtol) == 1, xtol > 0,
    is.numeric(stepmax), length(stepmax) == 1, stepmax > 0,
    is.numeric(maxeval), length(maxeval) == 1, maxeval >= 1
  )
  out <- list(grtol = grtol, xtol = xtol, stepmax = stepmax, maxeval = as.integer(maxeval))
  class(out) <- "lavaanpl_qn_control"
  out
}

#' Control options of the pairwise stochastic approximation
#'
#' The `SA` method. It draws item pairs at each iteration, uniformly and independently across iterations, and runs in
#' three phases: Adam, preconditioning, and stochastic gradient updates with a frozen gain and Polyak-Ruppert averaging. An
#' epoch is the number of iterations needed to draw as many pairs as there are in the model, one pass over the pairs.
#'
#' @param pairs_per_iteration Number of item pairs drawn at each iteration.
#' @param alpha Step length of the Adam phase, divided by the square root of the number of
#'   iterations per epoch.
#' @param adam_epochs Maximum number of epochs of the Adam phase. If `NULL`, 5, or as many as make 150 iterations on small
#'   models, whose epochs are short.
#' @param precondition_epochs Number of epochs of the preconditioning phase, which draws each pair once per epoch,
#'   recomputes the observed information as the gain after each epoch, and takes stochastic Newton steps with it. The
#'   gain is then frozen, unless the observed information is not positive definite. With `0`, the gain is computed at
#'   the end of the Adam phase.
#' @param max_precondition_epochs Maximum number of epochs of the preconditioning phase while the observed information
#'   is not positive definite.
#' @param average_epochs Number of epochs of the averaged phase, at least 1.
#' @param tol Absolute change of the pairwise log-likelihood at the epoch average, per free parameter, below which, at
#'   two consecutive epochs, the Adam or the preconditioning phase ends before its maximum.
#' @param decay Exponent of the step length decay in the averaged phase, in \eqn{(1/2, 1)}.
#' @param avg_from Fraction of the averaged phase discarded before averaging.
#' @param gain Gain of the preconditioning and averaged phases: `"hessian"` (inverse observed information),
#'   `"diagonal"` (Adam's diagonal preconditioner, no preconditioning phase), or `"auto"` (`"hessian"` up to 2000 free
#'   parameters).
#' @param seed Seed of the pair sampling.
#'
#' @return An object of class `lavaanpl_sa_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_sa_control <- function(
  pairs_per_iteration = 16,
  alpha = 0.2,
  adam_epochs = NULL,
  precondition_epochs = 3,
  max_precondition_epochs = 20,
  average_epochs = 50,
  tol = 0.05,
  decay = 0.75,
  avg_from = 0.25,
  gain = c("auto", "hessian", "diagonal"),
  seed = 123
) {
  gain <- match.arg(gain)
  num <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  stopifnot(
    num(pairs_per_iteration), pairs_per_iteration >= 1,
    num(alpha), alpha > 0,
    is.null(adam_epochs) || (num(adam_epochs) && adam_epochs >= 1),
    num(precondition_epochs), precondition_epochs >= 0,
    num(max_precondition_epochs), max_precondition_epochs >= 0,
    num(average_epochs), average_epochs >= 0,
    num(tol), tol > 0,
    num(decay), decay > 0.5, decay < 1,
    num(avg_from), avg_from >= 0, avg_from < 1,
    num(seed), seed > 0
  )
  out <- list(pairs_per_iteration = as.integer(pairs_per_iteration), alpha = alpha,
              adam_epochs = if (!is.null(adam_epochs)) as.integer(adam_epochs), precondition_epochs = as.integer(precondition_epochs),
              max_precondition_epochs = as.integer(max_precondition_epochs), average_epochs = as.integer(average_epochs),
              tol = tol, decay = decay, avg_from = avg_from, gain = gain,
              seed = as.integer(seed))
  class(out) <- "lavaanpl_sa_control"
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
#' @param nstarts Number of full runs from the same starting values, with seeds `seed`, `seed + 1` and so on. The highest
#'   pairwise log-likelihood is kept, with a warning if the runs disagree. Each run adds its full cost.
#'
#' @return An object of class `lavaanpl_saqn_control` to be passed to [lavaan.pl()] via `control`.
#'
#' @export
set_saqn_control <- function(
  pairs_per_iteration = 16,
  alpha = 0.2,
  adam_epochs = NULL,
  precondition_epochs = 3,
  max_precondition_epochs = 20,
  tol = 0.05,
  seed = 123,
  nstarts = 1,
  grtol = 1e-6,
  xtol = 1e-12,
  stepmax = 1,
  maxeval = 500
) {
  stopifnot(is.numeric(nstarts), length(nstarts) == 1, is.finite(nstarts), nstarts >= 1)
  out <- c(
    set_sa_control(pairs_per_iteration = pairs_per_iteration, alpha = alpha, adam_epochs = adam_epochs,
                   precondition_epochs = precondition_epochs, max_precondition_epochs = max_precondition_epochs,
                   average_epochs = 0, tol = tol, seed = seed),
    set_qn_control(grtol = grtol, xtol = xtol, stepmax = stepmax, maxeval = maxeval),
    list(nstarts = as.integer(nstarts))
  )
  class(out) <- c("lavaanpl_saqn_control", "lavaanpl_qn_control", "lavaanpl_sa_control")
  out
}

#' Method of a control object
#'
#' @param CONTROL Control object from [set_saqn_control()], [set_sa_control()] or [set_qn_control()].
#' @return `"SA-QN"`, `"SA"` or `"QN"`.
#' @seealso [fit_lavaanpl()], [pmlcfa()]
#' @keywords internal
control_method <- function(CONTROL) {
  if (inherits(CONTROL, "lavaanpl_saqn_control")) return("SA-QN")
  if (inherits(CONTROL, "lavaanpl_sa_control")) "SA" else "QN"
}

