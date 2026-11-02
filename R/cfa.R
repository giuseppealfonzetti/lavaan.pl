#' Fit Confirmatory Factor Analysis and Structural Equation Models
#'
#' Fit a Confirmatory Factor Analysis (CFA) model, or a Structural Equation Model (SEM) with
#' regressions among the latent variables and residual covariances of the items, using the
#' `lavaan.pl()` engine within the [`lavaan`] framework. By default, a quasi-Newton BFGS optimiser
#' from the [`ucminf`] package is used. For large data sets, the stochastic
#' approximation algorithm can be used.
#'
#' @details Not all [`lavaan`] options can be used at present. Some options of
#'   interest are:
#' \describe{
#'   \item{\code{information}}{The information matrix to use. Only `"observed"` is currently supported.}
#'   \item{\code{se}}{The standard error method to use. Only `"robust.huber.white"` is currently supported.}
#'   \item{\code{test}}{Default `"none"`. `"standard"` and `"mean.var.adjusted"` are computed by [compute_plrt()].}
#'   \item{\code{missing}}{Default `"listwise"`. `"pairwise"` and `"available.cases"` use the complete pairs and the available cases of Katsikatsou et al. (2022).}
#'   \item{\code{sampling.weights}, \code{cluster}}{Sampling weights, of type `"design"` or `"frequency"`, and cluster-robust standard errors, as in lavaan.}
#' }
#'
#' @param model A description of the user-specified model. Typically, the model
#'   is described using the lavaan model syntax. See [`lavaan::model.syntax`]
#'   for more information. Alternatively, a parameter table (eg. the output of
#'   the [`lavaan::lavaanify()`] function) is also accepted.
#' @param data A data frame containing the observed variables used in the model.
#'   Variables must be declared as ordered factors.
#' @param std.lv If `TRUE`, the metric of each latent variable is determined by
#'   fixing their (residual) variances to 1.0. If `FALSE`, the metric of each
#'   latent variable is determined by fixing the factor loading of the first
#'   indicator to 1.0.
#' @param start A vector of starting values in lavaan's free-parameter order
#'   (as returned by `coef()`). If not provided, lavaan's starting values are
#'   used and refined according to `estimator.args$init`.
#' @param estimator The estimator is `PML`. If any other estimator is provided,
#'   then [`lavaan::cfa()`] (or [`lavaan::sem()`]) is used.
#' @param estimator.args A list of arguments for [fit_latpl()]--see the help file
#'   for more details. Possible options are:
#'   \describe{
#'     \item{\code{init}}{Starting values object from [set_init()], by default `set_init(METHOD = "SA")`. Its `THETA` is taken from `start` or from the lavaan starting values.}
#'     \item{\code{ncores}}{The number of cores to use for parallel computation.}
#'     \item{\code{valdata}}{Validation data.}
#'     \item{\code{valgroup}}{Group of each row of `valdata`, with the values of the `group` variable. Required with `valdata` and several groups.}
#'     \item{\code{computevar_numderiv}}{If `TRUE`, the asymptotic variance-covariance matrix is computed using numerical derivatives.}
#'   }
#' @param control Control object of the estimation algorithm: [set_ucminf_control()]
#'   for the quasi-Newton BFGS optimiser (default), [set_sa_control()] for stochastic
#'   approximation.
#' @param verbose If `TRUE`, print additional information during the estimation
#'   process.
#' @param ... Additional arguments to be passed to [lavaan()].
#'
#' @returns A [lavaan::lavaan] object. The lavaan.pl fit, from [fit_latpl()], its
#'   variance components, from [compute_var()], and its test, from
#'   [compute_plrt()], are stored in `@external$lavaan.pl`, `@external$vars` and
#'   `@external$plrt`. Methods for `lavaan` objects work as usual.
#'
#' @references
#'
#' Katsikatsou, M., Moustaki, I., Yang-Wallentin, F., & Jöreskog, K. G. (2012).
#' Pairwise likelihood estimation for factor analysis models with ordinal data.
#' *Computational Statistics & Data Analysis*, *56*(12), 4243--4258.
#' <https://doi.org/10.1016/j.csda.2012.04.010>
#'
#' Katsikatsou, M., Moustaki, I., & Jamil, H. (2022). Pairwise likelihood
#' estimation for confirmatory factor analysis models with categorical variables
#' and data that are missing at random. *British Journal of Mathematical and
#' Statistical Psychology*, *75*(1), 23--45. <https://doi.org/10.1111/bmsp.12243>
#'
#' Alfonzetti, G., Bellio, R., Chen, Y., & Moustaki, I. (2025). Pairwise
#' stochastic approximation for confirmatory factor analysis of categorical
#' data. *British Journal of Mathematical and Statistical Psychology*, *78*(1),
#' 22--43. <https://doi.org/10.1111/bmsp.12347>
#'
#' @examples
#'
#' # A simple binary factor model using the LSAT data
#' fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
#' summary(fit)
#'
#' @export
cfa <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(
      init = set_init(METHOD = "SA"),
      ncores = 1,
      valdata = NULL,
      computevar_numderiv = FALSE
    ),
    start = NULL,
    control = NULL,
    verbose = FALSE,
    ...
) {
  fit_lavaan_pml("cfa", match.call(), model, data, std.lv, estimator, estimator.args, start, control, verbose, ...)
}

#' @rdname cfa
#' @export
sem <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(
      init = set_init(METHOD = "SA"),
      ncores = 1,
      valdata = NULL,
      computevar_numderiv = FALSE
    ),
    start = NULL,
    control = NULL,
    verbose = FALSE,
    ...
) {
  fit_lavaan_pml("sem", match.call(), model, data, std.lv, estimator, estimator.args, start, control, verbose, ...)
}

#' Fit a lavaan model with lavaan.pl
#'
#' @param LAVFUN Name of the lavaan function, `"cfa"` or `"sem"`.
#' @param CALL Call to store in the output.
#' @inheritParams cfa
#' @return An object of class `lavaan`.
#' @seealso [cfa()], [sem()], [spec_from_lavaan()], [as_lavaan()]
#' @keywords internal
fit_lavaan_pml <- function(LAVFUN, CALL, model, data, std.lv, estimator, estimator.args, start, control, verbose, ...) {

  # Validate arguments ---------------------------------------------------------
  if (!is.null(estimator.args$method)) {
    stop("`estimator.args$method` is no longer supported; the method follows `control`, e.g. `control = set_sa_control()`.", call. = FALSE)
  }
  if (any(c("init_method", "init_control", "cpp_control_init") %in% names(estimator.args))) {
    stop("`estimator.args$init_method`, `init_control` and `cpp_control_init` are no longer supported; use `init = set_init(...)`.", call. = FALSE)
  }
  if (is.null(control)) control <- set_ucminf_control()
  if (!inherits(control, c("latpl_ucminf_control", "latpl_sa_control"))) {
    stop("`control` must be created with `set_ucminf_control()` or `set_sa_control()`.", call. = FALSE)
  }
  init <- estimator.args$init
  if (is.null(init)) init <- set_init(METHOD = "SA")
  if (!inherits(init, "latpl_init")) {
    stop("`estimator.args$init` must be created with `set_init()`.", call. = FALSE)
  }
  if (!is.null(init$THETA)) {
    stop("Supply starting values through `start`, not `estimator.args$init`.", call. = FALSE)
  }
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  valdata <- estimator.args$valdata
  valgroup <- estimator.args$valgroup
  computevar_numderiv <- estimator.args$computevar_numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE

  # Initialise {lavaan} model object -------------------------------------------
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$std.lv <- std.lv
  lavargs$estimator <- estimator
  lavargs$do.fit <- FALSE

  if (estimator == "PML") {
    if ("information" %in% names(lavargs)) {
      if (lavargs$information != "observed") {
        stop("Only 'observed' information is currently supported.", call. = FALSE)
      }
    }
    lavargs$information <- "observed"
    if ("se" %in% names(lavargs)) {
      if (!(lavargs$se %in% c("robust.huber.white", "none"))) {
        message("Only 'robust.huber.white' is currently supported.")
        lavargs$se <- "robust.huber.white"
      }
    } else {
      lavargs$se <- "robust.huber.white"
    }
    if (!("test" %in% names(lavargs))) {
      lavargs$test <- "none"
    }
    template <- do.call(get(LAVFUN, envir = asNamespace("lavaan")), lavargs)
  } else{
    # In case user selects "DWLS" or "WLSMV", revert to lavaan::cfa()
    lavargs$do.fit <- TRUE
    fit <- do.call(get(LAVFUN, envir = asNamespace("lavaan")), lavargs)
    return(fit)
  }

  # Fit lavaan.pl -------------------------------------------------------------------
  spec <- spec_from_lavaan(template, start, valdata, valgroup)
  init$THETA <- spec$THETA
  # Send to lavaan.pl
  fit <- fit_latpl(
    DATA = spec$DATA,
    MODEL = spec$MODEL,
    CONTROL = control,
    INIT = init,
    NCORES = ncores,
    VERBOSE = verbose
  )

  # Compute standard errors ----------------------------------------------------
  vars <- NULL
  if (lavargs$se != "none") {
    vars <- compute_var(
      OBJ = fit,
      NUMDERIV = computevar_numderiv,
      VERBOSE = verbose
    )
  }

  out <- as_lavaan(fit, vars, template)
  out@call <- CALL
  out@call[[1]] <- call("::", quote(lavaan.pl), as.name(LAVFUN))
  out
}

#' lavaan.pl inputs from a lavaan model
#'
#' @param template Unfitted [lavaan::lavaan] object, created with `do.fit = FALSE`.
#' @param start Starting values in lavaan's free-parameter order, or `NULL`.
#' @param valdata Validation data passed to [set_latpl_data()] as `VALDATA`, or `NULL`.
#' @param valgroup Group of each row of `valdata`, passed to [set_latpl_data()] as `VALGROUP`, or `NULL`.
#' @return A list with the data object `DATA`, the model object `MODEL` and the starting values `THETA` for
#'   [fit_latpl()].
#' @seealso [fit_lavaan_pml()], [as_lavaan()]
#' @keywords internal
spec_from_lavaan <- function(template, start = NULL, valdata = NULL, valgroup = NULL) {
  pt <- lavaan::parTable(template)
  lv <- lavaan::lavNames(template, "lv")

  if (template@Data@nlevels > 1L) {
    stop("Multilevel models are not currently supported.", call. = FALSE)
  }
  # Check any ordinal data or not?
  if (!all(template@Data@ov$type == "ordered")) {
    stop("All measurement items must be declared as ordered factors.", call. = FALSE)
  }
  if (any(pt$op == "<~")) {
    stop("Composites (`<~`) are not currently supported.", call. = FALSE)
  }
  if (any(pt$op == "~" & !(pt$lhs %in% lv & pt$rhs %in% lv))) {
    stop("Regressions involving observed variables are not currently supported.", call. = FALSE)
  }
  if (nrow(template@Model@cin.JAC) > 0L || any(is.finite(c(pt$lower[pt$free > 0], pt$upper[pt$free > 0])))) {
    stop("Inequality constraints and parameter bounds are not currently supported.", call. = FALSE)
  }
  thetapar <- template@Options$parameterization == "theta"
  if (!thetapar && any(pt$free > 0 & pt$op == "~~" & pt$lhs == pt$rhs & !(pt$lhs %in% lv))) {
    stop("Free residual variances require parameterization = \"theta\".", call. = FALSE)
  }
  other <- pt$free > 0 & !(pt$op %in% c("=~", "~~", "|", "~", "~1", "~*~"))
  if (any(other)) {
    stop(paste0("Free `", paste(unique(pt$op[other]), collapse = "`, `"), "` parameters are not currently supported."), call. = FALSE)
  }

  X <- lapply(template@Data@X, function(x) x - 1)
  G <- length(X)
  missing <- template@Options$missing
  if (!missing %in% c("listwise", "pairwise", "available.cases")) {
    stop("missing = \"", missing, "\" is not currently supported; use \"listwise\", \"pairwise\" or \"available.cases\".", call. = FALSE)
  }
  nlev <- template@Data@ov$nlev[match(template@Data@ov.names[[1]], template@Data@ov$name)]
  if (any(vapply(X, function(x) any(apply(x, 2, max, na.rm = TRUE) + 1 != nlev), logical(1)))) {
    stop("Every category of every item must be observed at least once in every group.", call. = FALSE)
  }
  if (isTRUE(template@Model@ceq.simple.only)) {
    stop("ceq.simple = TRUE is not currently supported.", call. = FALSE)
  }
  if (length(template@Model@ceq.nonlinear.idx) > 0L) {
    stop("Nonlinear equality constraints are not currently supported.", call. = FALSE)
  }

  frees <- lapply(lavaan::lavInspect(template, "free", drop.list.single.group = FALSE), function(x) lapply(x, unclass))
  ests <- lapply(lavaan::lavInspect(template, "est", drop.list.single.group = FALSE), function(x) lapply(x, unclass))
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
    sfree <- if (thetapar) diag(free$theta) else as.vector(free$delta)
    if (!thetapar) {
      scales <- as.vector(est$delta)
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
    type[scfree] <- if (thetapar) "resvar" else "delta"
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
                       regmat = regmat, rescov = rescov, intercepts = intercepts, means = means, scales = scales)
    shift <- shift + d1 + length(nufree) + length(afree) + length(scfree)
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
      SCALES = s$scales
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
                              WEIGHTS_NORMALIZATION = if (normalization == "group") "group" else "total")

  # Note: the values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, latent variances, regressions, then residual covariances.
  if (is.null(start) && G > 1L) {
    init <- init_theta(set_init(), data_spec, latpl_dims(data_spec, model), model, FALSE)
    sgn <- unlist(Map(function(s, f) s$lambda[f$lambda > 0], lavaan::lavInspect(template, "start"), frees))
    lam <- theta_block(model, sum(data_spec$cat - 1)) == "lambda"
    init[lam] <- ifelse(sgn < 0, -1, 1) * init[lam]
  } else {
    if (is.null(start)) {
      starts <- list(lapply(lavaan::lavInspect(template, "start"), unclass))
      starts[[1]]$tau[] <- unlist(lapply(seq_len(ncol(X[[1]])), function(j) {
        stats::qnorm(cumsum(tabulate(X[[1]][, j] + 1, nlev[j]))[-nlev[j]] / sum(!is.na(X[[1]][, j])))
      }))
    } else {
      glist <- lavaan::lav_model_set_parameters(template@Model, as.numeric(start))@GLIST
      starts <- unname(split(glist, rep(seq_len(G), template@Model@nmat)))
    }
    pars <- lapply(starts, function(s) {
      list(thresholds = as.numeric(s$tau), loadings = s$lambda, latent_correlations = s$psi, regressions = s$beta,
           residual_covariances = s$theta, intercepts = as.vector(s$nu), means = as.vector(s$alpha),
           scales = as.vector(s$delta))
    })
    init <- par2theta(model, if (G == 1L) pars[[1]] else pars)
  }

  list(DATA = data_spec, MODEL = model, THETA = init)
}

#' Convert a lavaan.pl fit into a lavaan object
#'
#' @param fit Object of class `latpl_fit` returned by [fit_latpl()].
#' @param vars Output of [compute_var()] for `fit`, or `NULL` to skip standard errors.
#' @param template Unfitted [lavaan::lavaan] object for the same model and data, created with `do.fit = FALSE`.
#' @return A [lavaan::lavaan] object.
#' @export
as_lavaan <- function(fit, vars = NULL, template) {
  frees <- lapply(lavaan::lavInspect(template, "free", drop.list.single.group = FALSE), function(x) lapply(x, unclass))
  par_lav <- function(theta) {
    fit$theta <- theta
    pars <- get_par(fit, OPTION = "list")
    if (!inherits(fit$model, "latpl_groups")) pars <- list(pars)
    x <- numeric(template@Model@nx.free)
    for (g in seq_along(frees)) {
      free <- frees[[g]]
      par <- pars[[g]]
      x[free$tau[free$tau > 0]] <- par$thresholds[free$tau > 0]
      x[free$lambda[free$lambda > 0]] <- par$loadings[free$lambda > 0]
      x[free$psi[free$psi > 0]] <- par$latent_correlations[free$psi > 0]
      if (!is.null(free$beta)) x[free$beta[free$beta > 0]] <- par$regressions[free$beta > 0]
      x[free$theta[free$theta > 0]] <- par$residual_covariances[free$theta > 0]
      x[free$nu[free$nu > 0]] <- par$intercepts[free$nu > 0]
      x[free$alpha[free$alpha > 0]] <- par$means[free$alpha > 0]
      if (!is.null(par$scales)) x[free$delta[free$delta > 0]] <- par$scales[free$delta > 0]
    }
    x
  }

  # Change Model and implied slots ---------------------------------------------
  pt <- lavaan::parTable(template)
  idx <- pt$free > 0
  pt$est[idx] <- pt$start[idx] <- par_lav(fit$theta)[pt$free[idx]]
  options <- template@Options
  options[c("do.fit", "optim.method", "se", "test", "baseline")] <- list(TRUE, "none", "none", "none", FALSE)
  out <- lavaan::lavaan(
    model = pt,
    slot_options = options,
    slot_sample_stats = template@SampleStats,
    slot_data = template@Data,
    slot_cache = template@Cache,
    sloth1 = template@h1
  )
  out@Options[c("optim.method", "se", "test", "baseline")] <-
    list(if (fit$method == "SA") "SA" else "ucminf", template@Options$se, template@Options$test, template@Options$baseline)

  # Change Fit slot (depends whether it is numFit or stoFit) -------------------
  if (fit$method == "ucminf") {
    out@optim$iterations <- out@Fit@iterations <- as.integer(fit$numFit$info[["neval"]])
    out@optim$converged <- out@Fit@converged <- fit$numFit$convergence %in% c(1L, 2L) ||
      (fit$numFit$convergence == 4L && fit$numFit$info[["maxgradient"]] < template@Options$optim.dx.tol)
  } else {
    out@optim$iterations <- out@Fit@iterations <- as.integer(fit$stoFit$iterations)
    out@optim$converged <- out@Fit@converged <- fit$stoFit$converged
  }

  # Change vcov slot -----------------------------------------------------------
  vcov <- NULL
  if (!is.null(vars)) {
    delta <- numDeriv::jacobian(function(x) par_lav(as.vector(fit$lincon$K %*% x + fit$lincon$c)),
                                fit$theta[fit$lincon$free_idx])
    hinv <- solve(vars$H)
    vcov <- delta %*% hinv %*% vars$J %*% hinv %*% t(delta) / fit$dims$n
    attr(vcov, "E.inv") <- delta %*% hinv %*% t(delta)
    if (!is.null(vars$noise)) vcov <- vcov + delta %*% vars$noise %*% t(delta)
    out@vcov <- list(se = out@Options$se, information = out@Options$information[1], vcov = vcov)
    out@ParTable$se <- out@Fit@se <- lavaan___lav_model_vcov_se(out@Model, out@ParTable, vcov = vcov)
  }

  # Change test slot -----------------------------------------------------------
  plrt <- NULL
  if (!identical(out@Options$test, "none")) {
    plrt <- compute_plrt(fit, if (is.null(vars)) compute_var(fit) else vars)
    tests <- setdiff(out@Options$test, "standard")
    for (t in setdiff(tests, "mean.var.adjusted")) {
      warning(sprintf("test option %s not available for estimator PML", t), call. = FALSE)
    }
    pml_test <- function(x) {
      sg <- pmax(x$stat.group, 0)
      std <- sum(sg)
      out_test <- list(standard = list(test = "standard", stat = std, stat.group = sg, df = x$df,
                                       refdistr = "unknown", pvalue = NA))
      if ("mean.var.adjusted" %in% tests) {
        out_test$mean.var.adjusted <- list(
          test = "mean.var.adjusted", stat = x$stat.scaled, stat.group = sg / x$scaling.factor,
          df = x$df.scaled, pvalue = x$pvalue.scaled, scaling.factor = x$scaling.factor,
          scaled.test.stat = std, scaled.test = "standard", label = "mean+var adjusted correction (PML)",
          shift.parameter = NA, trace.UGamma = NA, trace.UGamma4 = NA, trace.UGamma2 = NA,
          UGamma.eigenvalues = NA
        )
      }
      attr(out_test, "info") <- lavaan___lav_test_info_attr(lavdata = out@Data, lavoptions = out@Options)
      out_test
    }
    out@test <- pml_test(plrt$model)
    if (isTRUE(out@Options$baseline)) out@baseline <- list(test = pml_test(plrt$baseline))
  }

  # Change timing slot ---------------------------------------------------------
  out@timing$optim <- fit$RTime
  out@timing$vcov <- if (is.null(vars)) 0 else vars$RTime
  out@timing$total <- out@timing$total + fit$RTime + out@timing$vcov
  # Include the entire output of lavaan.pl fit --------------------------------------
  out@external <- list(lavaan.pl = fit, vars = vars, plrt = plrt)
  # Output ---------------------------------------------------------------------
  out@call <- template@call
  out
}

#' Likelihood ratio test for nested models
#'
#' @param object Object of class `lavaan`.
#' @param ... Further objects of class `lavaan`, and arguments for [lavaan::lavTestLRT()].
#' @param model_names Names of the models.
#' @return An object of class `anova`, as [lavaan::lavTestLRT()]. When all models are fitted by
#'   [cfa()], the difference tests come from [compute_plrt()].
#' @export
lavTestLRT <- function(object, ..., model_names = NULL) {
  dots <- list(...)
  modp <- vapply(dots, inherits, logical(1), "lavaan")
  mods <- c(list(object), dots[modp])
  if (is.null(model_names)) {
    mc <- as.list(match.call(expand.dots = TRUE))[-1]
    model_names <- vapply(mc[c(TRUE, modp)], function(x) deparse1(x), "")
  }
  ours <- length(mods) > 1 &&
    all(vapply(mods, function(x) !is.null(x@external$lavaan.pl), logical(1))) &&
    (is.null(dots$method) || dots$method %in% c("default", "mean.var.adjusted.PLRT")) &&
    (is.null(dots$type) || tolower(dots$type) == "chisq")
  if (!ours) {
    return(lavaan::lavTestLRT(object, ..., model_names = model_names))
  }

  names(mods) <- make.unique(model_names)
  mods <- mods[order(vapply(mods, function(x) x@test[[1]]$df, numeric(1)))]
  fits <- lapply(mods, function(x) x@external$lavaan.pl)
  dim_theta <- vapply(mods, function(x) {
    vars <- if (is.null(x@external$vars)) compute_var(x@external$lavaan.pl) else x@external$vars
    sum(diag(solve(vars$H, vars$J)))
  }, numeric(1))
  logl <- vapply(mods, function(x) x@optim$logl, numeric(1))
  stat_delta <- df_delta <- rep(NA, length(mods))
  for (m in seq_len(length(mods) - 1)) {
    out <- compute_plrt(fits[[m + 1]], H1 = fits[[m]])$model
    stat_delta[m + 1] <- out$stat.scaled
    df_delta[m + 1] <- out$df.scaled
  }

  val <- data.frame(
    Df = vapply(mods, function(x) x@test[[1]]$df, numeric(1)),
    PL_AIC = -2 * logl + 2 * dim_theta,
    PL_BIC = -2 * logl + dim_theta * log(fits[[1]]$dims$n),
    Chisq = vapply(mods, function(x) x@test[[1]]$stat, numeric(1)),
    `Chisq diff` = stat_delta,
    `Df diff` = df_delta,
    `Pr(>Chisq)` = stats::pchisq(stat_delta, df_delta, lower.tail = FALSE),
    row.names = names(mods),
    check.names = FALSE
  )
  attr(val, "heading") <- paste0("\nMean and Variance Adjusted Chi-Squared Difference Test (method = ",
                                 dQuote("mean.var.adjusted.PLRT"), ", lavaan.pl)\n")
  class(val) <- c("anova", class(val))
  val
}

#' Analysis of variance for nested models
#'
#' @param object A fitted model object.
#' @param ... Further fitted model objects, and arguments for [lavTestLRT()] or [stats::anova()].
#' @return An object of class `anova`. For objects of class `lavaan`, the output of [lavTestLRT()];
#'   otherwise, the output of `lavaan::anova()`.
#' @export
anova <- function(object, ...) {
  if (!inherits(object, "lavaan")) return(lavaan::anova(object, ...))
  mc <- as.list(match.call(expand.dots = TRUE))[-1]
  modp <- vapply(list(...), inherits, logical(1), "lavaan")
  lavTestLRT(object, ..., model_names = vapply(mc[c(TRUE, modp)], function(x) deparse1(x), ""))
}

# Export lavaan functions
lavaan___lav_test_info_attr <-
  utils::getFromNamespace("lav_test_info_attr", "lavaan")
lavaan___lav_model_vcov_se <-
  utils::getFromNamespace("lav_model_vcov_se", "lavaan")
