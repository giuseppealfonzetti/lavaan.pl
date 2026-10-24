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
#'   \item{\code{test}}{Default `"none"`. `"standard"` and `"mean.var.adjusted"` are computed by lavaan from the lavaan.pl estimates.}
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
#' @returns A [lavaan::lavaan] object. The lavaan.pl fit, from [fit_plFA()], and
#'   its variance components, from [computeVar()], are stored in
#'   `@external$lavaan.pl` and `@external$vars`. Methods for `lavaan` objects
#'   work as usual.
#'
#' @references
#'
#' Katsikatsou, M., Moustaki, I., Yang-Wallentin, F., & Jöreskog, K. G. (2012).
#' Pairwise likelihood estimation for factor analysis models with ordinal data.
#' *Computational Statistics & Data Analysis*, *56*(12), 4243--4258.
#' <https://doi.org/10.1016/j.csda.2012.04.010>
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
  if (any(pt$op %in% c("<", ">"))) {
    stop("Inequality constraints are not currently supported.", call. = FALSE)
  }
  if (any(pt$free > 0 & pt$op == "~~" & !(pt$lhs %in% lv & pt$rhs %in% lv))) {
    stop("Free residual variances or covariances are not currently supported.", call. = FALSE)
  }
  other <- pt$free > 0 & !(pt$op %in% c("=~", "~~", "|"))
  if (any(other)) {
    stop(paste0("Free `", paste(unique(pt$op[other]), collapse = "`, `"), "` parameters are not currently supported."), call. = FALSE)
  }
  if (any(pt$free[pt$op == "|"] == 0)) {
    stop("Fixed thresholds are not currently supported.", call. = FALSE)
  }

  data <- template@Data@X[[1]] - 1
  if (anyNA(data)) {
    stop("Missing data are not currently supported.", call. = FALSE)
  }
  nlev <- template@Data@ov$nlev[match(template@Data@ov.names[[1]], template@Data@ov$name)]
  if (any(apply(data, 2, max) + 1 != nlev)) {
    stop("Every category of every item must be observed at least once.", call. = FALSE)
  }

  free <- lapply(lavaan::lavInspect(template, "free"), unclass)
  est <- lapply(lavaan::lavInspect(template, "est"), unclass)
  q <- ncol(free$lambda)

  # CORRFLAG
  off <- lower.tri(free$psi)
  if (all(free$psi[off] > 0)) {
    corrflag <- TRUE
  } else if (all(free$psi[off] == 0 & est$psi[off] == 0)) {
    corrflag <- FALSE
  } else {
    stop("Latent covariances must be either all free or all fixed to zero.", call. = FALSE)
  }
  constrvar <- diag(est$psi)
  constrvar[diag(free$psi) > 0] <- NA
  if (any(constrvar <= 0, na.rm = TRUE)) {
    stop("Fixed latent variances must be positive.", call. = FALSE)
  }

  # Build the p x q loading constraints matrix
  constrmat <- est$lambda
  constrmat[free$lambda > 0] <- NA

  # Linear loadings constraints
  llc <- NULL
  jac <- template@Model@ceq.JAC
  if (nrow(jac) > 0L) {
    if (length(template@Model@ceq.nonlinear.idx) > 0L || any(template@Model@ceq.rhs != 0)) {
      stop("Only homogeneous linear equality constraints are currently supported.", call. = FALSE)
    }
    load_idx <- free$lambda[free$lambda > 0]
    if (any(jac[, -load_idx] != 0)) {
      stop("Equality constraints are only supported on factor loadings.", call. = FALSE)
    }
    jac <- jac[, load_idx, drop = FALSE]
    piv <- integer(0)
    for (j in seq_len(ncol(jac))) {
      r <- length(piv) + 1L
      if (r > nrow(jac)) break
      i <- r - 1L + which.max(abs(jac[r:nrow(jac), j]))
      if (abs(jac[i, j]) < 1e-10) next
      jac[c(r, i), ] <- jac[c(i, r), ]
      jac[r, ] <- jac[r, ] / jac[r, j]
      jac[-r, ] <- jac[-r, , drop = FALSE] - outer(jac[-r, j], jac[r, ])
      piv <- c(piv, j)
    }
    coords <- t(vapply(load_idx, function(i) as.numeric(which(free$lambda == i, arr.ind = TRUE)), numeric(2)))
    llc <- lapply(seq_along(piv), function(r) {
      refs <- setdiff(which(abs(jac[r, ]) > 1e-10), piv)
      c(list(coords[piv[r], ]), lapply(refs, function(j) c(-jac[r, j], coords[j, ])))
    })
    constrmat[coords[piv, , drop = FALSE]] <- 1
  }

  # Build constraint list
  model <- set_plfa_model(
    CONSTRMAT = constrmat,
    CONSTRVAR = constrvar,
    CORRFLAG = corrflag,
    STDLV = !anyNA(constrvar),
    LLC = llc
  )
  if (!identical(is.na(model$CONSTRMAT), is.na(constrmat))) {
    stop("This identification of the latent variables is not currently supported.", call. = FALSE)
  }

  # Note: the values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, then latent variances.
  if (is.null(start)) {
    start_glist <- lapply(lavaan::lavInspect(template, "start"), unclass)
  } else {
    start_glist <- lavaan::lav_model_set_parameters(template@Model, as.numeric(start))@GLIST
  }
  init <- c(
    as.numeric(start_glist$tau),
    cpp_loadings_mat2vec(start_glist$lambda, constrmat, sum(is.na(constrmat))),
    cpp_latvar_mat2vec(start_glist$psi, log(sqrt(constrvar)), if (corrflag) q * (q - 1) / 2 else 0, sum(is.na(constrvar)))
  )

  list(DATA = set_plfa_data(data, valdata), MODEL = model, THETA = init)
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
    list("nlminb", template@Options$se, template@Options$test, template@Options$baseline)

  # Change Fit slot (depends whether it is numFit or stoFit) -------------------
  if (fit$method == "ucminf") {
    out@optim$iterations <- out@Fit@iterations <- as.integer(fit$numFit$info[["neval"]])
    out@optim$converged <- out@Fit@converged <- fit$numFit$convergence %in% c(1L, 2L) ||
      (fit$numFit$convergence == 4L && fit$numFit$info[["maxgradient"]] < template@Options$optim.dx.tol)
  } else {
    out@optim$iterations <- out@Fit@iterations <- as.integer(fit$stoFit$last_iter)
    out@optim$converged <- out@Fit@converged <- fit$stoFit$convergence$convergence_full == 1
  }

  # Change vcov slot -----------------------------------------------------------
  vcov <- NULL
  if (!is.null(vars)) {
    delta <- numDeriv::jacobian(par_lav, fit$theta)
    hinv <- solve(vars$H)
    vcov <- delta %*% hinv %*% vars$J %*% hinv %*% t(delta) / fit$dims$n
    attr(vcov, "E.inv") <- delta %*% hinv %*% t(delta)
    out@vcov <- list(se = out@Options$se, information = out@Options$information[1], vcov = vcov)
    out@ParTable$se <- out@Fit@se <- lavaan___lav_model_vcov_se(out@Model, out@ParTable, vcov = vcov)
  }

  # Change test slot -----------------------------------------------------------
  if (!identical(out@Options$test, "none")) {
    x <- out@optim$x
    attr(x, "fx") <- lavaan___lav_model_objective(
      lavmodel = out@Model,
      lavsamplestats = out@SampleStats,
      lavdata = out@Data,
      lavcache = out@Cache
    )
    out@test <- lavaan___lav_model_test(
      lavmodel = out@Model,
      lavpartable = out@ParTable,
      lavsamplestats = out@SampleStats,
      lavimplied = out@implied,
      lavh1 = out@h1,
      lavoptions = out@Options,
      x = x,
      vcov_1 = vcov,
      lavcache = out@Cache,
      lavdata = out@Data,
      lavloglik = out@loglik
    )
    out@baseline <- lavaan___lav_step15_baseline(
      lavoptions = out@Options,
      lavsamplestats = out@SampleStats,
      lavdata = out@Data,
      lavcache = out@Cache,
      lavh1 = out@h1,
      lavpartable = out@ParTable
    )
  }

  # Change timing slot ---------------------------------------------------------
  out@timing$optim <- fit$RTime
  out@timing$vcov <- if (is.null(vars)) 0 else vars$RTime
  out@timing$total <- out@timing$total + fit$RTime + out@timing$vcov
  # Include the entire output of lavaan.pl fit --------------------------------------
  out@external <- list(lavaan.pl = fit, vars = vars)
  # Output ---------------------------------------------------------------------
  out@call <- template@call
  out
}

# Export lavaan functions
lavaan___lav_step15_baseline <-
  utils::getFromNamespace("lav_step15_baseline", "lavaan")
lavaan___lav_model_test <-
  utils::getFromNamespace("lav_model_test", "lavaan")
lavaan___lav_model_vcov_se <-
  utils::getFromNamespace("lav_model_vcov_se", "lavaan")
lavaan___lav_model_objective <-
  utils::getFromNamespace("lav_model_objective", "lavaan")
