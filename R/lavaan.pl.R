#' Fit CFA and SEM Models by Pairwise Likelihood
#'
#' Fit a Confirmatory Factor Analysis (CFA) model using the `lavaan.pl()` engine
#' using the [`lavaan`] framework.
#'
#' @details Not all [`lavaan`] options can be used at present. Some options of
#'   interest are:
#' \describe{
#'   \item{\code{information}}{The information matrix to use. Only `"observed"` is currently supported.}
#'   \item{\code{se}}{The standard error method to use. Only `"robust.huber.white"` is currently supported.}
#' }
#'
#' @param model A description of the user-specified model. Typically, the model
#'   is described using the lavaan model syntax. See [`lavaan::model.syntax`]
#'   for more information. Alternatively, a parameter table (eg. the output of
#'   the [`lavaan::lavaanify()`] function) is also accepted.
#' @param data A data frame containing the observed variables used in the model.
#'   Ordinal items must be declared as ordered factors.
#' @param model.type The lavaan function whose defaults are used, `"sem"` or `"cfa"`.
#' @param std.lv If `TRUE`, the metric of each latent variable is determined by
#'   fixing their (residual) variances to 1.0. If `FALSE`, the metric of each
#'   latent variable is determined by fixing the factor loading of the first
#'   indicator to 1.0.
#' @param estimator The estimator, which must be `"PML"`.
#' @param estimator.args A list with `ncores` and `computevar_numderiv`.
#' @param start A vector of starting values in lavaan's free-parameter order.
#' @param control Control object from [set_saqn_control()], [set_sa_control()] or [set_qn_control()], or `NULL`.
#' @param verbose If `TRUE`, print additional information during the estimation
#'   process.
#' @param ... Additional arguments for [lavaan::lavaan()], including `optim.method` (`"SA-QN"`, `"SA"` or `"QN"`)
#'   and `samplestats = FALSE`, which skips lavaan's sample statistics (faster, for simulations).
#' @return An object of class [`lavaan.pl`][lavaan.pl-class].
#' @seealso [lavaan.pl-class]
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
#' fit <- pmlcfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
#' summary(fit)
#'
#' @export
lavaan.pl <- function(
    model,
    data,
    model.type = c("sem", "cfa"),
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = NULL,
    verbose = FALSE,
    ...
) {
  mc <- match.call()
  mc[[1L]] <- quote(lavaan.pl::lavaan.pl)
  model.type <- match.arg(model.type)

  # Validate arguments ---------------------------------------------------------
  lavargs <- list(...)
  names(lavargs) <- sub("^(std|optim|ceq|do)_", "\\1.", names(lavargs))
  if (!identical(estimator, "PML")) {
    stop("`estimator` must be \"PML\"; use lavaan::", model.type, "() for other estimators.", call. = FALSE)
  }
  if (!is.null(lavargs$std.lv)) std.lv <- lavargs$std.lv
  unknown <- setdiff(names(estimator.args), c("ncores", "computevar_numderiv"))
  if (length(unknown)) {
    stop("Unsupported `estimator.args`: ", paste(unknown, collapse = ", "), "; use `optim.method` for the method and `start` for starting values.", call. = FALSE)
  }
  method <- lavargs$optim.method
  if (!is.null(method) && !(length(method) == 1 && method %in% c("SA-QN", "SA", "QN"))) {
    stop("`optim.method` must be \"SA-QN\", \"SA\" or \"QN\".", call. = FALSE)
  }
  if (is.null(control)) {
    control <- switch(if (is.null(method)) "SA-QN" else method, `SA-QN` = set_saqn_control(), SA = set_sa_control(), QN = set_qn_control())
  }
  if (!inherits(control, c("latpl_qn_control", "latpl_sa_control"))) {
    stop("`control` must be created with `set_saqn_control()`, `set_sa_control()` or `set_qn_control()`.", call. = FALSE)
  }
  if (!is.null(method) && !identical(method, control_method(control))) {
    stop("`control` is for ", control_method(control), ", but `optim.method` is ", method, ".", call. = FALSE)
  }
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  computevar_numderiv <- estimator.args$computevar_numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$std.lv <- std.lv
  lavargs$estimator <- "PML"
  lavargs$do.fit <- FALSE
  lavargs$optim.method <- NULL
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
  if (is.null(lavargs$ceq.simple)) lavargs$ceq.simple <- FALSE
  if (isFALSE(lavargs$samplestats)) lavargs$slot_cache <- list()
  template <- do.call(get(model.type, envir = asNamespace("lavaan")), lavargs)

  # Fit lavaan.pl -------------------------------------------------------------------
  spec <- spec_from_lavaan(template, start)
  # Send to lavaan.pl
  fit <- fit_latpl(
    DATA = spec$DATA,
    MODEL = spec$MODEL,
    CONTROL = control,
    START = spec$START,
    NCORES = ncores,
    VERBOSE = verbose
  )

  # Compute standard errors ----------------------------------------------------
  vars <- NULL
  if (lavargs$se != "none") {
    vars <- tryCatch(
      compute_var(fit, NUMDERIV = computevar_numderiv),
      latpl_not_pd = function(e) if (fit$method == "SA") NULL else stop(e)
    )
  }

  out <- create_lav_from_latpl(fit, vars, template)
  out@call <- mc
  out
}

#' @rdname lavaan.pl
#' @export
pmlcfa <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = NULL,
    verbose = FALSE,
    ...
) {
  sc <- sys.call()
  sc[["model.type"]] <- "cfa"
  sc[[1L]] <- quote(lavaan.pl::lavaan.pl)
  eval(sc, parent.frame())
}

#' @rdname lavaan.pl
#' @export
pmlsem <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = NULL,
    verbose = FALSE,
    ...
) {
  sc <- sys.call()
  sc[["model.type"]] <- "sem"
  sc[[1L]] <- quote(lavaan.pl::lavaan.pl)
  eval(sc, parent.frame())
}
