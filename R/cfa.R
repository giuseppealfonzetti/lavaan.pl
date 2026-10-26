#' Fit Confirmatory Factor Analysis Models
#'
#' Fit a Confirmatory Factor Analysis (CFA) model using the `lavaan.pl()` engine
#' using the [`lavaan`] framework. By default, a quasi-Newton BFGS optimiser
#' from the [`ucminf`] package is used. For large data sets, the stochastic
#' approximation algorithm can be used.
#'
#' @details Not all [`lavaan`] options can be used at present. Some options of
#'   interest are:
#' \describe{
#'   \item{\code{information}}{The information matrix to use. Only `"observed"` is currently supported.}
#'   \item{\code{se}}{The standard error method to use. Only `"robust.huber.white"` is currently supported.}
#'   \item{\code{test}}{Default `"none"`. `"standard"` and `"mean.var.adjusted"` are computed by [computePLRT()].}
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
#'   then [`lavaan::cfa()`] is used.
#' @param estimator.args A list of arguments for [fit_plFA()]--see the help file
#'   for more details. Possible options are:
#'   \describe{
#'     \item{\code{init}}{Starting values object from [set_init()], by default `set_init(METHOD = "SA")`. Its `THETA` is taken from `start` or from the lavaan starting values.}
#'     \item{\code{ncores}}{The number of cores to use for parallel computation.}
#'     \item{\code{valdata}}{Validation data.}
#'     \item{\code{computevar_numderiv}}{If `TRUE`, the asymptotic variance-covariance matrix is computed using numerical derivatives.}
#'   }
#' @param control Control object of the estimation algorithm: [set_ucminf_control()]
#'   for the quasi-Newton BFGS optimiser (default), [set_sa_control()] for stochastic
#'   approximation.
#' @param verbose If `TRUE`, print additional information during the estimation
#'   process.
#' @param ... Additional arguments to be passed to [lavaan()].
#'
#' @returns A [lavaan::lavaan] object. The lavaan.pl fit, from [fit_plFA()], its
#'   variance components, from [computeVar()], and its test, from
#'   [computePLRT()], are stored in `@external$lavaan.pl`, `@external$vars` and
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

  # Validate arguments ---------------------------------------------------------
  if (!is.null(estimator.args$method)) {
    stop("`estimator.args$method` is no longer supported; the method follows `control`, e.g. `control = set_sa_control()`.", call. = FALSE)
  }
  if (any(c("init_method", "init_control", "cpp_control_init") %in% names(estimator.args))) {
    stop("`estimator.args$init_method`, `init_control` and `cpp_control_init` are no longer supported; use `init = set_init(...)`.", call. = FALSE)
  }
  if (is.null(control)) control <- set_ucminf_control()
  if (!inherits(control, c("plfa_ucminf_control", "plfa_sa_control"))) {
    stop("`control` must be created with `set_ucminf_control()` or `set_sa_control()`.", call. = FALSE)
  }
  init <- estimator.args$init
  if (is.null(init)) init <- set_init(METHOD = "SA")
  if (!inherits(init, "plfa_init")) {
    stop("`estimator.args$init` must be created with `set_init()`.", call. = FALSE)
  }
  if (!is.null(init$THETA)) {
    stop("Supply starting values through `start`, not `estimator.args$init`.", call. = FALSE)
  }
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  valdata <- estimator.args$valdata
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
    template <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)
  } else{
    # In case user selects "DWLS" or "WLSMV", revert to lavaan::cfa()
    lavargs$do.fit <- TRUE
    fit <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)
    return(fit)
  }

  # Fit lavaan.pl -------------------------------------------------------------------
  spec <- spec_from_lavaan(template, start, valdata)
  init$THETA <- spec$THETA
  # Send to lavaan.pl
  fit <- fit_plFA(
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
    vars <- computeVar(
      OBJ = fit,
      NUMDERIV = computevar_numderiv,
      VERBOSE = verbose
    )
  }

  out <- as_lavaan(fit, vars, template)
  out@call <- match.call()
  out@call[[1]] <- quote(lavaan.pl::cfa)
  out
}

spec_from_lavaan <- function(template, start = NULL, valdata = NULL) {
  pt <- lavaan::parTable(template)
  lv <- lavaan::lavNames(template, "lv")

  if (template@Data@ngroups > 1L) {
    stop("Multigroup analysis is not currently supported.", call. = FALSE)
  }
  if (template@Data@nlevels > 1L) {
    stop("Multilevel models are not currently supported.", call. = FALSE)
  }
  # Check any ordinal data or not?
  if (!all(template@Data@ov$type == "ordered")) {
    stop("All measurement items must be declared as ordered factors.", call. = FALSE)
  }
  if (template@Options$parameterization != "delta") {
    stop("Only the delta parameterization is currently supported.", call. = FALSE)
  }
  if (any(pt$op %in% c("~", "<~")) || any(pt$op == "=~" & pt$rhs %in% lv)) {
    stop("Structural regressions and higher-order factors are not currently supported.", call. = FALSE)
  }
  if (nrow(template@Model@cin.JAC) > 0L || any(is.finite(c(pt$lower[pt$free > 0], pt$upper[pt$free > 0])))) {
    stop("Inequality constraints and parameter bounds are not currently supported.", call. = FALSE)
  }
  if (any(pt$free > 0 & pt$op == "~~" & !(pt$lhs %in% lv & pt$rhs %in% lv))) {
    stop("Free residual variances or covariances are not currently supported.", call. = FALSE)
  }
  other <- pt$free > 0 & !(pt$op %in% c("=~", "~~", "|"))
  if (any(other)) {
    stop(paste0("Free `", paste(unique(pt$op[other]), collapse = "`, `"), "` parameters are not currently supported."), call. = FALSE)
  }

  data <- template@Data@X[[1]] - 1
  missing <- template@Options$missing
  if (!missing %in% c("listwise", "pairwise", "available.cases")) {
    stop("missing = \"", missing, "\" is not currently supported; use \"listwise\", \"pairwise\" or \"available.cases\".", call. = FALSE)
  }
  nlev <- template@Data@ov$nlev[match(template@Data@ov.names[[1]], template@Data@ov$name)]
  if (any(apply(data, 2, max, na.rm = TRUE) + 1 != nlev)) {
    stop("Every category of every item must be observed at least once.", call. = FALSE)
  }

  free <- lapply(lavaan::lavInspect(template, "free"), unclass)
  est <- lapply(lavaan::lavInspect(template, "est"), unclass)
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

  # Linear equality constraints A theta = b, in lavaan.pl's parameter order
  if (isTRUE(template@Model@ceq.simple.only)) {
    stop("ceq.simple = TRUE is not currently supported.", call. = FALSE)
  }
  if (length(template@Model@ceq.nonlinear.idx) > 0L) {
    stop("Nonlinear equality constraints are not currently supported.", call. = FALSE)
  }
  nthr <- length(est$tau)
  nload <- sum(is.na(constrmat))
  ncorr <- if (corrflag) q * (q - 1) / 2 else 0
  vfree <- diag(free$psi) > 0
  d <- nthr + nload + ncorr + sum(vfree)
  pos <- type <- rep(NA, template@Model@nx.free)
  pos[free$tau[free$tau > 0]] <- which(free$tau > 0)
  type[free$tau[free$tau > 0]] <- "tau"
  pos[free$lambda[free$lambda > 0]] <- nthr + seq_len(nload)
  type[free$lambda[free$lambda > 0]] <- "lambda"
  pos[diag(free$psi)[vfree]] <- nthr + nload + ncorr + seq_len(sum(vfree))
  type[diag(free$psi)[vfree]] <- "var"
  type[cov_x] <- "cov"
  cpos <- matrix(NA, q, q)
  if (corrflag) {
    cpos[upper.tri(cpos)] <- nthr + nload + seq_len(ncorr)
    cpos[lower.tri(cpos)] <- t(cpos)[lower.tri(cpos)]
    pos[cov_x] <- cpos[off][free$psi[off] > 0]
  }
  jac <- template@Model@ceq.JAC
  rhs <- template@Model@ceq.rhs
  A <- matrix(0, 0, d)
  b <- numeric(0)
  for (r in seq_len(nrow(jac))) {
    nz <- which(jac[r, ] != 0)
    row <- numeric(d)
    if (any(type[nz] == "cov") && any(type[nz] == "var")) {
      stop("Constraints mixing latent covariances and variances are not currently supported.", call. = FALSE)
    } else if (all(type[nz] == "var")) {
      k <- if (length(nz) == 2L && rhs[r] == 0) -jac[r, nz[2]] / jac[r, nz[1]] else NA
      if (!isTRUE(k > 0)) {
        stop("Constraints on latent variances must be equalities or ratios (v1 == k*v2).", call. = FALSE)
      }
      row[pos[nz]] <- c(1, -1)
      A <- rbind(A, row)
      b <- c(b, log(k) / 2)
    } else if (any(type[nz] == "var")) {
      stop("Constraints mixing latent variances with other parameters are not currently supported.", call. = FALSE)
    } else {
      row[pos[nz]] <- jac[r, nz]
      A <- rbind(A, row)
      b <- c(b, rhs[r])
    }
  }
  for (k in which(free$tau == 0)) {
    row <- numeric(d)
    row[k] <- 1
    A <- rbind(A, row)
    b <- c(b, est$tau[k])
  }
  if (latcov == "direct") {
    for (k in which(off & free$psi == 0)) {
      row <- numeric(d)
      row[cpos[k]] <- 1
      A <- rbind(A, row)
      b <- c(b, est$psi[k])
    }
  }

  # Build constraint list
  model <- set_plfa_model(
    CONSTRMAT = constrmat,
    CONSTRVAR = constrvar,
    CORRFLAG = corrflag,
    STDLV = !anyNA(constrvar),
    CONSTRAINTS = if (nrow(A) > 0L) list(A = unname(A), b = b),
    LATCOV = latcov
  )
  if (!identical(is.na(model$CONSTRMAT), is.na(constrmat))) {
    stop("This identification of the latent variables is not currently supported.", call. = FALSE)
  }

  # Note: the values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, then latent variances.
  if (is.null(start)) {
    start_glist <- lapply(lavaan::lavInspect(template, "start"), unclass)
    start_glist$tau[] <- unlist(lapply(seq_len(ncol(data)), function(j) {
      stats::qnorm(cumsum(tabulate(data[, j] + 1, nlev[j]))[-nlev[j]] / sum(!is.na(data[, j])))
    }))
  } else {
    start_glist <- lavaan::lav_model_set_parameters(template@Model, as.numeric(start))@GLIST
  }
  init <- c(
    as.numeric(start_glist$tau),
    cpp_loadings_mat2vec(start_glist$lambda, constrmat, sum(is.na(constrmat))),
    cpp_latvar_mat2vec(start_glist$psi, log(sqrt(constrvar)), if (corrflag) q * (q - 1) / 2 else 0, sum(is.na(constrvar)),
                       corrflag_code(model))
  )

  weights <- NULL
  if (length(template@Data@sampling.weights) > 0L) {
    if (template@Options$sampling.weights.normalization == "none") {
      stop("sampling.weights.normalization = \"none\" is not currently supported.", call. = FALSE)
    }
    weights <- template@Data@weights[[1]]
  }
  cluster <- if (isTRUE(template@Options$.clustered)) template@Data@Lp[[1]]$cluster.idx[[2]] else NULL
  data_spec <- set_plfa_data(data, valdata, MISSING = missing, WEIGHTS = weights,
                             WEIGHTS_TYPE = template@Options$sampling.weights.type, CLUSTER = cluster)

  list(DATA = data_spec, MODEL = model, THETA = init)
}

#' Convert a lavaan.pl fit into a lavaan object
#'
#' @param fit Object of class `PlFaFit` returned by [fit_plFA()].
#' @param vars Output of [computeVar()] for `fit`, or `NULL` to skip standard errors.
#' @param template Unfitted [lavaan::lavaan] object for the same model and data, created with `do.fit = FALSE`.
#' @return A [lavaan::lavaan] object.
#' @export
as_lavaan <- function(fit, vars = NULL, template) {
  free <- lapply(lavaan::lavInspect(template, "free"), unclass)
  par_lav <- function(theta) {
    fit$theta <- theta
    par <- getPar(fit, OPTION = "list")
    x <- numeric(template@Model@nx.free)
    x[free$tau[free$tau > 0]] <- par$thresholds[free$tau > 0]
    x[free$lambda[free$lambda > 0]] <- par$loadings[free$lambda > 0]
    x[free$psi[free$psi > 0]] <- par$latent_correlations[free$psi > 0]
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
    out@optim$iterations <- out@Fit@iterations <- as.integer(fit$stoFit$last_iter)
    out@optim$converged <- out@Fit@converged <- fit$stoFit$convergence$convergence_full == 1 ||
      (!fit$control$CHECK_TOL && fit$stoFit$convergence$convergence_full == 0)
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
    plrt <- computePLRT(fit, if (is.null(vars)) computeVar(fit) else vars)
    tests <- setdiff(out@Options$test, "standard")
    for (t in setdiff(tests, "mean.var.adjusted")) {
      warning(sprintf("test option %s not available for estimator PML", t), call. = FALSE)
    }
    pml_test <- function(x) {
      std <- max(x$stat, 0)
      out_test <- list(standard = list(test = "standard", stat = std, stat.group = std, df = x$df,
                                       refdistr = "unknown", pvalue = NA))
      if ("mean.var.adjusted" %in% tests) {
        out_test$mean.var.adjusted <- list(
          test = "mean.var.adjusted", stat = x$stat.scaled, stat.group = std / x$scaling.factor,
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
#'   [cfa()], the difference tests come from [computePLRT()].
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
    vars <- if (is.null(x@external$vars)) computeVar(x@external$lavaan.pl) else x@external$vars
    sum(diag(solve(vars$H, vars$J)))
  }, numeric(1))
  logl <- vapply(mods, function(x) x@optim$logl, numeric(1))
  stat_delta <- df_delta <- rep(NA, length(mods))
  for (m in seq_len(length(mods) - 1)) {
    out <- computePLRT(fits[[m + 1]], H1 = fits[[m]])$model
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
