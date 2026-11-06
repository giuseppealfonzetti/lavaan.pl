#' lavaan.pl inputs from a lavaan model
#'
#' @param template Unfitted [lavaan::lavaan] object, created with `do.fit = FALSE`.
#' @param start Starting values in lavaan's free-parameter order, or `NULL`.
#' @return A list with the data object `DATA`, the model object `MODEL` and the starting values `START` for
#'   [fit_lavaanpl()].
#' @seealso [lavaan.pl()], [create_lav_from_lavaanpl()]
#' @keywords internal
spec_from_lavaan <- function(template, start = NULL) {
  pt <- lavaan::parTable(template)
  lv <- lavaan::lavNames(template, "lv")
  ovx <- lavaan::lavNames(template, "ov.x")

  if (template@Data@nlevels > 1L) {
    stop("Multilevel models are not currently supported.", call. = FALSE)
  }
  # Check any ordinal data or not?
  if (length(ovx) && !isTRUE(template@Options$conditional.x)) {
    stop("conditional.x = FALSE is not currently supported.", call. = FALSE)
  }
  if (length(ovx) && !isTRUE(template@Options$fixed.x)) {
    stop("fixed.x = FALSE is not currently supported.", call. = FALSE)
  }
  numeric_items <- template@Data@ov$type[match(template@Data@ov.names[[1]], template@Data@ov$name)] != "ordered"
  if (any(pt$op == "<~")) {
    stop("Composites (`<~`) are not currently supported.", call. = FALSE)
  }
  if (any(pt$op == "~" & !(pt$rhs %in% ovx | (pt$lhs %in% lv & pt$rhs %in% lv)))) {
    stop("Regressions involving observed variables other than exogenous covariates are not currently supported.", call. = FALSE)
  }
  if (nrow(template@Model@cin.JAC) > 0L || any(is.finite(c(pt$lower[pt$free > 0], pt$upper[pt$free > 0])))) {
    stop("Inequality constraints and parameter bounds are not currently supported.", call. = FALSE)
  }
  thetapar <- template@Options$parameterization == "theta"
  if (!thetapar && any(pt$free > 0 & pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% template@Data@ov.names[[1]][!numeric_items])) {
    stop("Free residual variances require parameterization = \"theta\".", call. = FALSE)
  }
  other <- pt$free > 0 & !(pt$op %in% c("=~", "~~", "|", "~", "~1", "~*~"))
  if (any(other)) {
    stop(paste0("Free `", paste(unique(pt$op[other]), collapse = "`, `"), "` parameters are not currently supported."), call. = FALSE)
  }

  X <- lapply(template@Data@X, function(x) {
    x[, !numeric_items] <- x[, !numeric_items] - 1
    x
  })
  G <- length(X)
  missing <- template@Options$missing
  if (!missing %in% c("listwise", "pairwise", "available.cases")) {
    stop("missing = \"", missing, "\" is not currently supported; use \"listwise\", \"pairwise\" or \"available.cases\".", call. = FALSE)
  }
  nlev <- template@Data@ov$nlev[match(template@Data@ov.names[[1]], template@Data@ov$name)]
  if (any(vapply(X, function(x) any((apply(x, 2, max, na.rm = TRUE) + 1 != nlev)[!numeric_items]), logical(1)))) {
    stop("Every category of every item must be observed at least once in every group.", call. = FALSE)
  }
  if (isTRUE(template@Model@ceq.simple.only)) {
    stop("ceq.simple = TRUE is not currently supported.", call. = FALSE)
  }
  if (length(template@Model@ceq.nonlinear.idx) > 0L) {
    stop("Nonlinear equality constraints are not currently supported.", call. = FALSE)
  }

  weights <- NULL
  normalization <- template@Options$sampling.weights.normalization
  if (length(template@Data@sampling.weights) > 0L) {
    if (normalization == "none") {
      stop("sampling.weights.normalization = \"none\" is not currently supported.", call. = FALSE)
    }
    weights <- unlist(template@Data@weights)
  }
  cluster <- if (isTRUE(template@Options$.clustered)) unlist(lapply(template@Data@Lp, function(x) x$cluster.idx[[2]]))
  data_spec <- lavaanpl_data(Reduce(rbind, X), MISSING = missing, WEIGHTS = weights,
                              WEIGHTS_TYPE = template@Options$sampling.weights.type, CLUSTER = cluster,
                              GROUP = if (G > 1L) rep(template@Data@group.label, vapply(X, nrow, numeric(1))),
                              WEIGHTS_NORMALIZATION = if (normalization == "group") "group" else "total",
                              X = if (length(ovx)) unname(Reduce(rbind, template@Data@eXo)),
                              NUMERIC = numeric_items)

  frees <- lavaan_matrices(template, "free")
  ests <- lavaan_matrices(template, "est")
  jac <- template@Model@ceq.JAC
  rhs <- template@Model@ceq.rhs
  models <- Map(function(free, est) {
    mask <- Map(function(f, e) unname(replace(e, f > 0, NA)), free, est)

    # Latent correlations: LKJ when all free and unconstrained, orthogonal when all
    # fixed to zero, direct covariances otherwise
    off <- lower.tri(free$psi)
    cov_con <- nrow(jac) > 0L && any(jac[, free$psi[off][free$psi[off] > 0], drop = FALSE] != 0)
    lavaanpl_model(
      LAMBDA = mask$lambda,
      PSI = mask$psi,
      BETA = mask$beta,
      THETA = mask$theta,
      NU = as.vector(mask$nu),
      ALPHA = as.vector(mask$alpha),
      DELTA = if (!thetapar) as.vector(mask$delta),
      TAU = as.vector(mask$tau),
      GAMMA = if (ncol(mask$gamma)) mask$gamma,
      KAPPA = if (ncol(mask$kappa)) mask$kappa,
      LATCOV = if (!all(free$psi[off] > 0) || cov_con) "direct" else "lkj",
      PARAMETERIZATION = template@Options$parameterization,
      NUMERIC = numeric_items
    )
  }, frees, ests)

  # Linear equality constraints A parameter_vector = b, in lavaan.pl's parameter order
  ncat <- ifelse(numeric_items, 0, nlev - 1)
  symmetric <- function(M, diagonal, other) {
    ifelse(row(M) == col(M), sprintf(diagonal, row(M)), sprintf(other, pmin(row(M), col(M)), pmax(row(M), col(M))))
  }
  label <- character(template@Model@nx.free)
  for (g in seq_len(G)) {
    free <- frees[[g]]
    cell_names <- list(
      tau = sprintf("tau[%d,%d]", rep(seq_along(ncat), ncat), sequence(ncat)),
      lambda = sprintf("lambda[%d,%d]", row(free$lambda), col(free$lambda)),
      psi = symmetric(free$psi, "logsd[%d]", "psi[%d,%d]"),
      beta = if (!is.null(free$beta)) sprintf("beta[%d,%d]", row(free$beta), col(free$beta)),
      theta = symmetric(free$theta, "logrsd[%d]", "theta[%d,%d]"),
      nu = sprintf("nu[%d]", seq_along(free$nu)),
      alpha = sprintf("alpha[%d]", seq_along(free$alpha)),
      delta = if (!is.null(free$delta)) sprintf("logdelta[%d]", seq_along(free$delta)),
      gamma = sprintf("gamma[%d,%d]", row(free$gamma), col(free$gamma)),
      kappa = sprintf("kappa[%d,%d]", row(free$kappa), col(free$kappa))
    )
    for (m in names(cell_names)) {
      at <- free[[m]] > 0
      label[free[[m]][at]] <- paste0(cell_names[[m]][at], if (g > 1) paste0(".g", g))
    }
  }
  parameters <- parameter_names(data_spec, lavaanpl_groups(models))
  pos <- match(label, parameters)
  type <- sub("\\[.*", "", label)
  A <- matrix(0, 0, length(parameters))
  b <- numeric(0)
  for (r in seq_len(nrow(jac))) {
    nz <- which(jac[r, ] != 0)
    row <- numeric(ncol(A))
    logt <- type[nz] %in% c("logsd", "logrsd", "logdelta")
    if (any(type[nz] == "psi") && any(type[nz] == "logsd")) {
      stop("Constraints mixing latent covariances and variances are not currently supported.", call. = FALSE)
    } else if (all(logt) && length(unique(type[nz])) == 1L) {
      k <- if (length(nz) == 2L && rhs[r] == 0) -jac[r, nz[2]] / jac[r, nz[1]] else NA
      if (!isTRUE(k > 0)) {
        stop("Constraints on variances and scales must be equalities or ratios (v1 == k*v2).", call. = FALSE)
      }
      row[pos[nz]] <- c(1, -1)
      A <- rbind(A, row)
      b <- c(b, log(k) / if (type[nz[1]] == "logdelta") 1 else 2)
    } else if (any(logt)) {
      stop("Constraints mixing variances or scales with other parameters are not currently supported.", call. = FALSE)
    } else {
      row[pos[nz]] <- jac[r, nz]
      A <- rbind(A, row)
      b <- c(b, rhs[r])
    }
  }
  model <- lavaanpl_groups(models, if (nrow(A) > 0L) list(A = unname(A), b = b))

  # Note: the values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, latent variances, regressions, then residual covariances.
  init <- NULL
  if (!is.null(start)) {
    glist <- lavaan::lav_model_set_parameters(template@Model, as.numeric(start))@GLIST
    init <- matrices_to_vector(model, lavaan_matrices(template, unname(split(glist, rep(seq_len(G), template@Model@nmat)))))
  }

  list(DATA = data_spec, MODEL = model, START = init)
}

#' Model matrices of a lavaan model, for lavaan.pl
#'
#' @param template Unfitted [lavaan::lavaan] object.
#' @param WHAT `"free"`, `"est"` or `"start"` for [lavaan::lavInspect()], or a list of the model matrices of each group.
#' @return A list with the matrices of each group, without lavaan's phantom latent variables for items regressed on
#'   covariates, and with `gamma` (latent variables on covariates) and `kappa` (items on covariates).
#' @seealso [spec_from_lavaan()], [create_lav_from_lavaanpl()]
#' @keywords internal
lavaan_matrices <- function(template, WHAT) {
  groups <- if (is.character(WHAT)) lavaan::lavInspect(template, WHAT, drop.list.single.group = FALSE) else WHAT
  items <- template@Data@ov.names[[1]]
  lv <- colnames(lavaan::lavInspect(template, "free", drop.list.single.group = FALSE)[[1]]$lambda)
  ph <- which(lv %in% items)
  keep <- setdiff(seq_along(lv), ph)
  nx <- NCOL(template@Data@eXo[[1]])
  if (is.null(template@Data@eXo[[1]])) nx <- 0
  lapply(groups, function(x) {
    m <- lapply(x, unclass)
    if (is.null(m$gamma)) m$gamma <- matrix(0, length(lv), nx)
    m$kappa <- matrix(0, length(items), nx)
    for (j in ph) {
      i <- match(lv[j], items)
      if (
        any(m$beta[keep, j] != 0) ||
          any(m$beta[ph, j] != 0) ||
          any(m$beta[j, ph] != 0) ||
          any(m$psi[j, -j] != 0)
      ) {
        stop("Regressions involving observed variables other than exogenous covariates are not currently supported.", call. = FALSE)
      }
      m$lambda[i, keep] <- m$beta[j, keep]
      m$kappa[i, ] <- m$gamma[j, ]
      m$theta[i, i] <- m$psi[j, j]
      m$nu[i] <- m$nu[i] + m$alpha[j]
    }
    if (length(ph)) {
      m$lambda <- m$lambda[, keep, drop = FALSE]
      m$psi <- m$psi[keep, keep, drop = FALSE]
      m$beta <- m$beta[keep, keep, drop = FALSE]
      m$alpha <- m$alpha[keep, , drop = FALSE]
      m$gamma <- m$gamma[keep, , drop = FALSE]
    }
    m
  })
}
