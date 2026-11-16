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
#'   Items are ordinal unless `ordered` says otherwise.
#' @param model.type The lavaan function whose defaults are used, `"sem"` or `"cfa"`.
#' @param std.lv If `TRUE`, the metric of each latent variable is determined by
#'   fixing their (residual) variances to 1.0. If `FALSE`, the metric of each
#'   latent variable is determined by fixing the factor loading of the first
#'   indicator to 1.0.
#' @param ordered `TRUE` to treat all items as ordinal, `FALSE` for the lavaan default (only the ordered factors of
#'   `data`), or the names of the ordinal items. With `TRUE`, numeric columns with non-integer values or more than 10
#'   categories are an error.
#' @param estimator The estimator, which must be `"PML"`.
#' @param estimator.args A list with `ncores` and `computevar_numderiv`.
#' @param start A vector of starting values in lavaan's free-parameter order.
#' @param control Control object from [set_control()], which also selects the method.
#' @param verbose If `TRUE`, print additional information during the estimation
#'   process.
#' @param ... Additional arguments for [lavaan::lavaan()], including
#'   for faster fits in simulations, `test = "none"`, which skips the mean and variance adjusted test, and
#'   `samplestats = FALSE`, which skips lavaan's sample statistics.
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
    ordered = TRUE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = set_control(),
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
    stop("Unsupported `estimator.args`: ", paste(unknown, collapse = ", "), "; use `control = set_control()` for the method and `start` for starting values.", call. = FALSE)
  }
  if (!is.null(lavargs$optim.method)) {
    stop("`optim.method` is not supported: choose the method with `control = set_control()`.", call. = FALSE)
  }
  if (!inherits(control, c("lavaanpl_qn_control", "lavaanpl_sa_control"))) {
    stop("`control` must be created with `set_control()`.", call. = FALSE)
  }
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  computevar_numderiv <- estimator.args$computevar_numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE

  if (isTRUE(ordered)) {
    items <- intersect(lavaan::lavNames(if (is.character(model)) lavaan::lavaanify(model) else model, "ov.nox"), names(data))
    values <- Filter(is.numeric, as.data.frame(data)[items])
    continuous <- names(values)[vapply(values, function(x) any(x != round(x), na.rm = TRUE) || length(unique(x[!is.na(x)])) > 10, logical(1))]
    if (length(continuous)) {
      stop("With `ordered = TRUE`, ", paste(continuous, collapse = ", "), " would be ordinal, but ", if (length(continuous) > 1) "they have" else "it has",
           " non-integer values or more than 10 categories: use `ordered = FALSE` for numeric items, or list the ordinal items in `ordered`.", call. = FALSE)
    }
  }

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$std.lv <- std.lv
  lavargs$ordered <- if (isFALSE(ordered)) NULL else ordered
  lavargs$estimator <- "PML"
  lavargs$do.fit <- FALSE
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
    lavargs$test <- "mean.var.adjusted"
  }
  if (is.null(lavargs$ceq.simple)) lavargs$ceq.simple <- FALSE
  if (isFALSE(lavargs$samplestats)) lavargs$slot_cache <- list()
  template <- do.call(get(model.type, envir = asNamespace("lavaan")), lavargs)

  # Fit lavaan.pl -------------------------------------------------------------------
  spec <- spec_from_lavaan(template, start)
  # Send to lavaan.pl
  fit <- fit_lavaanpl(
    DATA = spec$DATA,
    MODEL = spec$MODEL,
    CONTROL = control,
    START = spec$START,
    NCORES = ncores,
    VERBOSE = verbose,
    CHECK = lavargs$se != "none" || !identical(lavargs$test, "none")
  )

  # Compute standard errors ----------------------------------------------------
  vars <- NULL
  if (lavargs$se != "none") {
    vars <- tryCatch(
      compute_var(fit, NUMDERIV = computevar_numderiv),
      lavaanpl_not_pd = function(e) if (fit$method == "SA") NULL else stop(e)
    )
  }

  out <- create_lav_from_lavaanpl(fit, vars, template)
  out@call <- mc
  out
}

#' @rdname lavaan.pl
#' @export
pmlcfa <- function(
    model,
    data,
    std.lv = FALSE,
    ordered = TRUE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = set_control(),
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
    ordered = TRUE,
    estimator = "PML",
    estimator.args = list(ncores = 1, computevar_numderiv = FALSE),
    start = NULL,
    control = set_control(),
    verbose = FALSE,
    ...
) {
  sc <- sys.call()
  sc[["model.type"]] <- "sem"
  sc[[1L]] <- quote(lavaan.pl::lavaan.pl)
  eval(sc, parent.frame())
}
