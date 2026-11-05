#' lavaan.pl inputs from a lavaan model
#'
#' @param template Unfitted [lavaan::lavaan] object, created with `do.fit = FALSE`.
#' @param start Starting values in lavaan's free-parameter order, or `NULL`.
#' @param valdata Validation data passed to [set_latpl_data()] as `VALDATA`, or `NULL`.
#' @param valgroup Group of each row of `valdata`, passed to [set_latpl_data()] as `VALGROUP`, or `NULL`.
#' @param valx Covariates of `valdata`, passed to [set_latpl_data()] as `VALX`, or `NULL`.
#' @return A list with the data object `DATA`, the model object `MODEL` and the starting values `THETA` for
#'   [fit_latpl()].
#' @seealso [lavaan.pl()], [create_lav_from_latpl()]
#' @keywords internal
spec_from_lavaan <- function(template, start = NULL, valdata = NULL, valgroup = NULL, valx = NULL) {
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

  frees <- lavaan_matrices(template, "free")
  ests <- lavaan_matrices(template, "est")
  pos <- type <- rep(NA, template@Model@nx.free)
  fix <- fixval <- numeric(0)
  specs <- vector("list", G)
  shift <- 0
  for (g in seq_len(G)) {
    free <- frees[[g]]
    est <- ests[[g]]
    q <- ncol(free$lambda)

    # Latent correlations: LKJ when all free and unconstrained, orthogonal when all
    # fixed to zero, direct covariances otherwise
    off <- lower.tri(free$psi)
    cov_x <- free$psi[off][free$psi[off] > 0]
    cov_con <- nrow(template@Model@ceq.JAC) > 0L && any(template@Model@ceq.JAC[, cov_x, drop = FALSE] != 0)
    corrflag <- !all(free$psi[off] == 0 & est$psi[off] == 0)
    latcov <- if (corrflag && (!all(free$psi[off] > 0) || cov_con)) "direct" else "lkj"
    constrvar <- diag(est$psi)
    constrvar[diag(free$psi) > 0] <- NA
    if (any(constrvar <= 0, na.rm = TRUE)) {
      stop("Fixed latent variances must be positive.", call. = FALSE)
    }

    # Build the p x q loading constraints matrix
    constrmat <- est$lambda
    constrmat[free$lambda > 0] <- NA

    # Structural regressions and residual covariances: NA when free
    regmat <- rescov <- NULL
    if (!is.null(free$beta) && (any(free$beta > 0) || any(est$beta != 0))) {
      regmat <- est$beta
      regmat[free$beta > 0] <- NA
    }
    up <- upper.tri(free$theta)
    if (thetapar || any(free$theta[up] > 0) || any(est$theta[up] != 0)) {
      rescov <- est$theta
      rescov[free$theta > 0] <- NA
    }
    intercepts <- as.vector(est$nu)
    intercepts[free$nu > 0] <- NA
    means <- as.vector(est$alpha)
    means[free$alpha > 0] <- NA
    scales <- NULL
    sfree <- if (thetapar) diag(free$theta) else ifelse(numeric_items, diag(free$theta), as.vector(free$delta))
    if (!thetapar) {
      scales <- ifelse(numeric_items, sqrt(pmax(diag(est$theta), 0)), as.vector(est$delta))
      scales[sfree > 0] <- NA
    }

    nthr <- length(est$tau)
    nload <- sum(is.na(constrmat))
    ncorr <- if (corrflag) q * (q - 1) / 2 else 0
    vfree <- diag(free$psi) > 0
    bfree <- if (is.null(regmat)) numeric(0) else free$beta[free$beta > 0]
    tfree <- free$theta[up][free$theta[up] > 0]
    d0 <- nthr + nload + ncorr + sum(vfree)
    nufree <- free$nu[free$nu > 0]
    afree <- free$alpha[free$alpha > 0]
    scfree <- sfree[sfree > 0]
    d1 <- d0 + length(bfree) + length(tfree)
    d2 <- d1 + length(nufree) + length(afree) + length(scfree)
    gamma <- est$gamma
    kappa <- est$kappa
    gamma[free$gamma > 0] <- NA
    kappa[free$kappa > 0] <- NA
    gfree <- free$gamma[free$gamma > 0]
    kfree <- free$kappa[free$kappa > 0]
    pos[gfree] <- shift + d2 + seq_along(gfree)
    type[gfree] <- "gamma"
    pos[kfree] <- shift + d2 + length(gfree) + seq_along(kfree)
    type[kfree] <- "kappa"
    pos[free$tau[free$tau > 0]] <- shift + which(free$tau > 0)
    type[free$tau[free$tau > 0]] <- "tau"
    pos[free$lambda[free$lambda > 0]] <- shift + nthr + seq_len(nload)
    type[free$lambda[free$lambda > 0]] <- "lambda"
    pos[diag(free$psi)[vfree]] <- shift + nthr + nload + ncorr + seq_len(sum(vfree))
    type[diag(free$psi)[vfree]] <- "var"
    type[cov_x] <- "cov"
    pos[bfree] <- shift + d0 + seq_along(bfree)
    type[bfree] <- "beta"
    pos[tfree] <- shift + d0 + length(bfree) + seq_along(tfree)
    type[tfree] <- "rescov"
    pos[nufree] <- shift + d1 + seq_along(nufree)
    type[nufree] <- "nu"
    pos[afree] <- shift + d1 + length(nufree) + seq_along(afree)
    type[afree] <- "alpha"
    pos[scfree] <- shift + d1 + length(nufree) + length(afree) + seq_along(scfree)
    type[scfree] <- if (thetapar) "resvar" else ifelse(numeric_items[sfree > 0], "resvar", "delta")
    cpos <- matrix(NA, q, q)
    if (corrflag) {
      cpos[upper.tri(cpos)] <- shift + nthr + nload + seq_len(ncorr)
      cpos[lower.tri(cpos)] <- t(cpos)[lower.tri(cpos)]
      pos[cov_x] <- cpos[off][free$psi[off] > 0]
    }
    fix <- c(fix, shift + which(free$tau == 0))
    fixval <- c(fixval, est$tau[free$tau == 0])
    if (latcov == "direct") {
      fix <- c(fix, cpos[off & free$psi == 0])
      fixval <- c(fixval, est$psi[off & free$psi == 0])
    }
    specs[[g]] <- list(constrmat = constrmat, constrvar = constrvar, corrflag = corrflag, latcov = latcov,
                       regmat = regmat, rescov = rescov, intercepts = intercepts, means = means, scales = scales,
                       gamma = gamma, kappa = kappa)
    shift <- shift + d2 + length(gfree) + length(kfree)
  }

  # Linear equality constraints A theta = b, in lavaan.pl's parameter order
  d <- shift
  jac <- template@Model@ceq.JAC
  rhs <- template@Model@ceq.rhs
  A <- matrix(0, 0, d)
  b <- numeric(0)
  for (r in seq_len(nrow(jac))) {
    nz <- which(jac[r, ] != 0)
    row <- numeric(d)
    logt <- type[nz] %in% c("var", "resvar", "delta")
    if (any(type[nz] == "cov") && any(type[nz] == "var")) {
      stop("Constraints mixing latent covariances and variances are not currently supported.", call. = FALSE)
    } else if (all(logt) && length(unique(type[nz])) == 1L) {
      k <- if (length(nz) == 2L && rhs[r] == 0) -jac[r, nz[2]] / jac[r, nz[1]] else NA
      if (!isTRUE(k > 0)) {
        stop("Constraints on variances and scales must be equalities or ratios (v1 == k*v2).", call. = FALSE)
      }
      row[pos[nz]] <- c(1, -1)
      A <- rbind(A, row)
      b <- c(b, log(k) / if (type[nz[1]] == "delta") 1 else 2)
    } else if (any(logt)) {
      stop("Constraints mixing variances or scales with other parameters are not currently supported.", call. = FALSE)
    } else {
      row[pos[nz]] <- jac[r, nz]
      A <- rbind(A, row)
      b <- c(b, rhs[r])
    }
  }
  for (k in seq_along(fix)) {
    row <- numeric(d)
    row[fix[k]] <- 1
    A <- rbind(A, row)
    b <- c(b, fixval[k])
  }

  # Build constraint list
  con <- if (nrow(A) > 0L) list(A = unname(A), b = b)
  models <- lapply(specs, function(s) {
    m <- set_latpl_model(
      CONSTRMAT = s$constrmat,
      CONSTRVAR = s$constrvar,
      CORRFLAG = s$corrflag,
      STDLV = !anyNA(s$constrvar),
      CONSTRAINTS = if (G == 1L) con,
      LATCOV = s$latcov,
      REGMAT = unname(s$regmat),
      RESCOV = unname(s$rescov),
      PARAMETERIZATION = template@Options$parameterization,
      INTERCEPTS = s$intercepts,
      MEANS = s$means,
      SCALES = s$scales,
      GAMMA = if (ncol(s$gamma)) unname(s$gamma),
      KAPPA = if (ncol(s$kappa)) unname(s$kappa),
      NUMERIC = numeric_items
    )
    if (!identical(is.na(m$CONSTRMAT), is.na(s$constrmat))) {
      stop("This identification of the latent variables is not currently supported.", call. = FALSE)
    }
    m
  })
  model <- if (G == 1L) models[[1]] else set_latpl_groups(models, con)

  weights <- NULL
  normalization <- template@Options$sampling.weights.normalization
  if (length(template@Data@sampling.weights) > 0L) {
    if (normalization == "none") {
      stop("sampling.weights.normalization = \"none\" is not currently supported.", call. = FALSE)
    }
    weights <- unlist(template@Data@weights)
  }
  cluster <- if (isTRUE(template@Options$.clustered)) unlist(lapply(template@Data@Lp, function(x) x$cluster.idx[[2]]))
  data_spec <- set_latpl_data(Reduce(rbind, X), valdata, MISSING = missing, WEIGHTS = weights,
                              WEIGHTS_TYPE = template@Options$sampling.weights.type, CLUSTER = cluster,
                              GROUP = if (G > 1L) rep(template@Data@group.label, vapply(X, nrow, numeric(1))),
                              VALGROUP = valgroup,
                              WEIGHTS_NORMALIZATION = if (normalization == "group") "group" else "total",
                              X = if (length(ovx)) unname(Reduce(rbind, template@Data@eXo)), VALX = valx,
                              NUMERIC = numeric_items)

  # Note: the values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, latent variances, regressions, then residual covariances.
  init <- NULL
  if (!is.null(start)) {
    glist <- lavaan::lav_model_set_parameters(template@Model, as.numeric(start))@GLIST
    starts <- lavaan_matrices(template, unname(split(glist, rep(seq_len(G), template@Model@nmat))))
    pars <- lapply(starts, function(s) {
      list(thresholds = as.numeric(s$tau), loadings = s$lambda, latent_correlations = s$psi, regressions = s$beta,
           residual_covariances = s$theta, intercepts = as.vector(s$nu), means = as.vector(s$alpha),
           scales = as.vector(s$delta), gamma = s$gamma, kappa = s$kappa)
    })
    init <- par2theta(model, if (G == 1L) pars[[1]] else pars)
  }

  list(DATA = data_spec, MODEL = model, THETA = init)
}

#' Model matrices of a lavaan model, for lavaan.pl
#'
#' lavaan represents the regression of an item on a covariate with a phantom latent variable named after the item.
#' Here the phantom is removed: its row of `gamma` is the row of the item in `kappa`, its row of `beta` that of the
#' item in `lambda`, its variance the residual variance of the item and its mean the intercept of the item.
#'
#' @param template Unfitted [lavaan::lavaan] object.
#' @param WHAT `"free"`, `"est"` or `"start"` for [lavaan::lavInspect()], or a list of the model matrices of each group.
#' @return A list with the matrices of each group, including `gamma` (latent variables on covariates) and `kappa`
#'   (items on covariates).
#' @seealso [spec_from_lavaan()], [create_lav_from_latpl()]
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
