#' Ordinal and numeric data for pairwise likelihood estimation
#'
#' @param DATA \eqn{n*p} data matrix. Ordinal items are coded as integers from zero.
#' @param MISSING Treatment of missing values: `"listwise"` drops incomplete observations, `"pairwise"` uses the
#'   complete pairs and `"available.cases"` adds the univariate terms of the available cases.
#' @param WEIGHTS Sampling weights, one per observation, normalised to sum to the number of observations. If `NULL`,
#'   unit weights.
#' @param WEIGHTS_TYPE Type of sampling weights, as in lavaan: `"design"` or `"frequency"`.
#' @param CLUSTER Cluster of each observation, for cluster-robust standard errors. If `NULL`, independent observations.
#' @param GROUP Group of each observation, numbered in order of appearance. If `NULL`, a single group.
#' @param WEIGHTS_NORMALIZATION Normalisation of the sampling weights, as in lavaan: `"total"` or `"group"`.
#' @param X Matrix of covariates, one row per observation. Observations with missing covariates are removed. If
#'   `NULL`, no covariates.
#' @param NUMERIC Logical vector marking the numeric items. If `NULL`, all items are ordinal.
#'
#' @return An object of class `lavaanpl_data` with the data matrix `Y`, its dimensions `n` and `p`, the number of item
#'   pairs `pairs`, the categories of each item `cat` (one for numeric items), the indicator `numeric`, the pairwise
#'   frequency tables `freq`, the treatment of missing values `missing`, the available-cases weights `uweights` and
#'   univariate frequencies `ufreq`, the normalised sampling weights `weights` and their score factors `jweights`, the
#'   clusters `cluster`, the groups `group` (from zero) and their number `ngroups`, the covariates `X`, the stratum of
#'   each observation `stratum` (from zero), the distinct covariate rows of each group `strata`, whether each group uses
#'   frequency tables `tables`, and the time spent on the frequencies `freqTime`. The tables are stacked by group and,
#'   within each group, by stratum.
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
#' @param STRATA Number of distinct covariate rows.
#' @param Y Data matrix.
#' @param CAT Number of categories of each item.
#' @param STRATUM Stratum of each row of `Y`.
#' @return `TRUE` if the frequency tables fit in memory and need fewer bivariate normal evaluations than the individual
#'   observations.
#' @seealso [lavaanpl_data()], [cpp_inputs()]
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
#' @param LAMBDA \eqn{p*q} matrix of loadings, `NA` where free. A latent variable with a free variance needs a fixed
#'   nonzero loading.
#' @param PSI \eqn{q*q} symmetric matrix of latent (disturbance) covariances, `NA` where free. If `NULL`, unit
#'   variances and free covariances.
#' @param BETA \eqn{q*q} matrix of regressions among the latent variables, `NA` where free, with `BETA[i, j]` the
#'   effect of `j` on `i`. If `NULL`, none.
#' @param THETA \eqn{p*p} symmetric matrix of residual covariances, `NA` where free. The diagonal is used for the
#'   numeric items and, with `PARAMETERIZATION = "theta"`, for all items. If `NULL`, no residual covariances, free
#'   residual variances of the numeric items and unit ones otherwise.
#' @param NU \eqn{p} vector of intercepts, `NA` where free. If `NULL`, free for the numeric items and zero otherwise.
#' @param ALPHA \eqn{q} vector of latent means, `NA` where free. If `NULL`, all zero.
#' @param DELTA \eqn{p} vector of scale factors, `NA` where free, used with `PARAMETERIZATION = "delta"` for the
#'   ordinal items. If `NULL`, all one.
#' @param TAU Vector of thresholds, item after item, `NA` where free. If `NULL`, all free.
#' @param GAMMA \eqn{q*x} matrix of regressions of the latent variables on the covariates, `NA` where free. If `NULL`,
#'   none.
#' @param KAPPA \eqn{p*x} matrix of direct effects of the covariates on the items, `NA` where free. If `NULL`, none.
#' @param CONSTRAINTS Linear equality constraints `A %*% parameter_vector = b`, as `list(A = , b = )`, with one column
#'   of `A` per parameter of [parameter_names()]. Each constraint involves a single block of [parameter_blocks()].
#' @param LATCOV Parametrisation of the latent correlations: `"lkj"` (always positive definite, all covariances free)
#'   or `"direct"` (the covariances themselves, which can be fixed or constrained).
#' @param PARAMETERIZATION `"delta"` or `"theta"`, as in lavaan.
#' @param NUMERIC Logical vector marking the numeric items, as in [lavaanpl_data()].
#'
#' @return An object of class `lavaanpl_model` with the matrices `LAMBDA`, `PSI`, `BETA`, `THETA`, `NU`, `ALPHA`,
#'   `DELTA`, `TAU`, `GAMMA` and `KAPPA` (unused entries set to zero or one), `CONSTRAINTS`, `LATCOV` (`"none"` for
#'   independent latent variables), `PARAMETERIZATION`, `NUMERIC`, and the numbers of items `p` and latent variables
#'   `q`.
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
#' @param CONSTRAINTS Linear equality constraints `A %*% parameter_vector = b` across groups, as `list(A = , b = )`,
#'   with one column of `A` per parameter of [parameter_names()], group after group. Each model's own `CONSTRAINTS`
#'   apply within its group.
#'
#' @return An object of class `lavaanpl_groups` with `MODELS`, `CONSTRAINTS`, the number of items `p` and the number
#'   of groups `ngroups`, to be passed to [fit_lavaanpl()] via `MODEL`.
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

#' Control options of the estimation
#'
#' @param method Estimation method: `"SA-QN"` (quasi-Newton from the burn-in of stochastic approximation), `"SA"`
#'   (pairwise stochastic approximation, which needs latent variables scaled by `std.lv = TRUE`) or `"QN"`
#'   (quasi-Newton).
#' @param ... Options of `method`, each with a default. `"QN"` takes the quasi-Newton options, `"SA"` the stochastic
#'   approximation options, and `"SA-QN"` both, except `halvings`, `average_epochs`, `decay` and `avg_from`, plus its own options.
#'
#'   Quasi-Newton options:
#'   \describe{
#'     \item{`grtol`}{Tolerance on the largest absolute gradient element. Default `1e-6`.}
#'     \item{`xtol`}{Tolerance on the relative step length. Default `1e-12`.}
#'     \item{`stepmax`}{Initial maximal step length of the trust region. Default `1`.}
#'     \item{`maxeval`}{Maximum number of function evaluations. Default `500`.}
#'   }
#'
#'   Stochastic approximation options:
#'   \describe{
#'     \item{`pairs_per_iteration`}{Number of item pairs drawn at each iteration. Default `16`.}
#'     \item{`alpha`}{Step length of the Adam phase. Default `0.2`.}
#'     \item{`adam_epochs`}{Number of epochs (passes over the item pairs) of the Adam phase. With `"SA"`, the phase is
#'       split into `halvings + 1` levels of equal length, and the step is halved after each level but the last, from
#'       the best point so far. If `NULL`, the default, 5 epochs and at least 50 iterations per level.}
#'     \item{`halvings`}{Number of halvings of the Adam step. Default `2`.}
#'     \item{`average_epochs`}{Number of epochs of the averaged phase, at least 1. If `NULL`, the default, one epoch and
#'       at least 150 iterations. The averaged phase uses Adam's diagonal preconditioner, frozen.}
#'     \item{`decay`}{Exponent of the step length decay in the averaged phase, between 0.5 and 1. Default `0.75`.}
#'     \item{`avg_from`}{Fraction of the averaged phase discarded before averaging. Default `0.25`.}
#'     \item{`seed`}{Seed of the pair sampling. Default `123`.}
#'   }
#'
#'   `"SA-QN"` also takes:
#'   \describe{
#'     \item{`adam_epochs`}{With `"SA-QN"`, the maximum number of epochs of the Adam phase, which stops early when the
#'       pairwise log-likelihood changes by less than `tol` per free parameter in two consecutive epochs. If `NULL`, the
#'       default, enough for 5 epochs and 150 iterations.}
#'     \item{`tol`}{Change of the pairwise log-likelihood, per free parameter, below which the burn-in phases stop
#'       early. Default `0.05`.}
#'     \item{`precondition_epochs`}{Number of epochs of the Newton preconditioning phase after Adam, up to 2000 free
#'       parameters. With `0`, the inverse observed information is computed at the end of the Adam phase and starts
#'       quasi-Newton. Default `3`.}
#'     \item{`max_precondition_epochs`}{Maximum number of preconditioning epochs while the observed information is not
#'       positive definite. Default `20`.}
#'     \item{`nstarts`}{Number of runs, with seeds `seed`, `seed + 1` and so on. The highest pairwise log-likelihood is
#'       kept, with a warning if the runs disagree. Default `1`.}
#'   }
#'
#' @return A control object to be passed to [lavaan.pl()] via `control`.
#' @seealso [lavaan.pl()], [diagnostics()], [ucminf::ucminf()]
#' @examples
#' fit <- pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE, control = set_control("QN", maxeval = 200))
#' @export
set_control <- function(method = c("SA-QN", "SA", "QN"), ...) {
  method <- match.arg(method)
  switch(method, `SA-QN` = set_saqn_control(...), SA = set_sa_control(...), QN = set_qn_control(...))
}

#' @noRd
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

#' @noRd
set_sa_control <- function(
  pairs_per_iteration = 16,
  alpha = 0.2,
  adam_epochs = NULL,
  halvings = 2,
  average_epochs = NULL,
  decay = 0.75,
  avg_from = 0.25,
  seed = 123
) {
  num <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  stopifnot(
    num(pairs_per_iteration), pairs_per_iteration >= 1,
    num(alpha), alpha > 0,
    is.null(adam_epochs) || (num(adam_epochs) && adam_epochs >= 1),
    num(halvings), halvings >= 0, halvings == round(halvings),
    is.null(average_epochs) || (num(average_epochs) && average_epochs >= 1),
    num(decay), decay > 0.5, decay < 1,
    num(avg_from), avg_from >= 0, avg_from < 1,
    num(seed), seed > 0
  )
  out <- list(pairs_per_iteration = as.integer(pairs_per_iteration), alpha = alpha,
              adam_epochs = if (!is.null(adam_epochs)) as.integer(adam_epochs), halvings = as.integer(halvings),
              average_epochs = if (!is.null(average_epochs)) as.integer(average_epochs),
              decay = decay, avg_from = avg_from, seed = as.integer(seed))
  class(out) <- "lavaanpl_sa_control"
  out
}

#' @noRd
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
  stopifnot(
    is.numeric(precondition_epochs), length(precondition_epochs) == 1, is.finite(precondition_epochs), precondition_epochs >= 0,
    is.numeric(max_precondition_epochs), length(max_precondition_epochs) == 1, is.finite(max_precondition_epochs),
    max_precondition_epochs >= 0,
    is.numeric(tol), length(tol) == 1, is.finite(tol), tol > 0,
    is.numeric(nstarts), length(nstarts) == 1, is.finite(nstarts), nstarts >= 1
  )
  burn <- set_sa_control(pairs_per_iteration = pairs_per_iteration, alpha = alpha, adam_epochs = adam_epochs, seed = seed)
  burn$halvings <- NULL
  burn$average_epochs <- 0L
  burn$tol <- tol
  out <- c(
    burn,
    list(precondition_epochs = as.integer(precondition_epochs), max_precondition_epochs = as.integer(max_precondition_epochs)),
    set_qn_control(grtol = grtol, xtol = xtol, stepmax = stepmax, maxeval = maxeval),
    list(nstarts = as.integer(nstarts))
  )
  class(out) <- c("lavaanpl_saqn_control", "lavaanpl_qn_control", "lavaanpl_sa_control")
  out
}

#' Method of a control object
#'
#' @param CONTROL Control object from [set_control()].
#' @return `"SA-QN"`, `"SA"` or `"QN"`.
#' @seealso [fit_lavaanpl()], [pmlcfa()]
#' @keywords internal
control_method <- function(CONTROL) {
  if (inherits(CONTROL, "lavaanpl_saqn_control")) return("SA-QN")
  if (inherits(CONTROL, "lavaanpl_sa_control")) "SA" else "QN"
}

