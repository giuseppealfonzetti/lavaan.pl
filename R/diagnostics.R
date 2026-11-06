#' Fit behind a stochastic approximation object
#'
#' @param object Object of class lavaanpl_internal, or a `lavaan.pl` object from [lavaan.pl()].
#' @return The object of class lavaanpl_internal.
#' @seealso [sa_diagnostics()], [plot_sa()], [sa_stability()]
#' @keywords internal
sa_fit <- function(object) {
  fit <- if (inherits(object, "lavaan.pl")) get_lavaanpl_internal(object) else object
  if (!inherits(fit, "lavaanpl_internal")) stop("`object` must come from pmlcfa() or pmlsem().", call. = FALSE)
  if (fit$method != "SA") stop("`object` was not fitted by stochastic approximation.", call. = FALSE)
  fit
}

#' Parameters of a stochastic approximation fit
#'
#' @param FIT Object of class lavaanpl_internal, from stochastic approximation.
#' @return A data frame with one row per free parameter: `name`, `block`, `estimate`, standard error `se`, and, when the
#'   observed information is positive definite, standard deviation of the optimisation noise `noise`, its ratio to the
#'   standard error `ratio` and the distance from the optimum in standard errors `distance`.
#' @seealso [sa_diagnostics()]
#' @keywords internal
sa_parameters <- function(FIT) {
  free <- FIT$lincon$free_idx
  out <- data.frame(name = parameter_names(FIT$data, FIT$model), block = FIT$blocks, estimate = transform_parameters(FIT$model, FIT$parameter_vector),
                    se = NA, noise = NA, ratio = NA, distance = NA)
  vars <- tryCatch(compute_var(FIT), lavaanpl_not_pd = function(e) NULL)
  if (!is.null(vars)) {
    out$se <- sqrt(vars$asymptotic_variance / FIT$data$n)
    out$noise <- sqrt(pmax(vars$optimisation_noise, 0))
    out$ratio <- ifelse(out$se > 0, out$noise / out$se, 0)
    out$distance <- vars$distance
  }
  out[free, ]
}

#' Diagnostics of the stochastic approximation
#'
#' Reports whether the run behaved, and how precise its estimate is. The optimisation noise is the variability due to the
#' pairs sampled, on top of the statistical variability; the distance is how far the estimate is from the exact pairwise
#' likelihood estimate, predicted from the gradient without fitting it.
#'
#' @param object A `lavaan.pl` object from [lavaan.pl()], fitted by stochastic approximation.
#' @param target Largest ratio of the optimisation noise to the standard error used to suggest a number of epochs.
#' @return An object of class `lavaanpl_sa_diagnostics`: a list with `converged`, the checks `checks`, the summary of the
#'   phases `phases`, the settings `settings`, the table of free parameters `parameters` from `sa_parameters()`, its
#'   summary by block `blocks`, and the value of `average_epochs` predicted to reach `target` in `epochs_needed`.
#' @export
sa_diagnostics <- function(object, target = 0.1) {
  fit <- sa_fit(object)
  stopifnot(is.numeric(target), length(target) == 1, target > 0)
  st <- fit$stoFit
  pars <- sa_parameters(fit)
  phases <- data.frame(
    phase = c("Adam", "Preconditioning", "Averaged"),
    epochs = c(st$adam_epochs, st$precondition_epochs, fit$control$average_epochs),
    iterations = c(st$adam_epochs, st$precondition_epochs, fit$control$average_epochs) * st$epoch_length,
    seconds = unname(c(st$time[["adam"]], st$time[["precondition"]], st$time[["averaging"]]))
  )
  summarise <- function(x) {
    data.frame(block = x$block[1], parameters = nrow(x), median_ratio = stats::median(x$ratio), max_ratio = max(x$ratio),
               max_distance = max(x$distance))
  }
  worst <- suppressWarnings(max(pars$ratio, na.rm = TRUE))
  out <- list(
    converged = st$converged,
    checks = c(plateau = st$plateau, skipped = st$skipped_averaging, rejected = st$rejections, positive_definite = st$check$pd,
               distance = st$check$distance),
    phases = phases,
    settings = list(method = fit$method, pairs = st$pairs_per_iteration, total_pairs = fit$data$pairs,
                    epoch_length = st$epoch_length, gain = st$gain, averaged = st$averaged, seconds = fit$RTime),
    parameters = pars,
    blocks = do.call(rbind, lapply(split(pars, factor(pars$block, unique(pars$block))), summarise)),
    worst_ratio = worst,
    epochs_needed = if (is.finite(worst) && worst > 0) ceiling(fit$control$average_epochs * (worst / target)^2) else NA,
    target = target
  )
  rownames(out$blocks) <- NULL
  class(out) <- "lavaanpl_sa_diagnostics"
  out
}

#' @export
print.lavaanpl_sa_diagnostics <- function(x, ...) {
  s <- x$settings
  cat(s$method, "fit:", if (isTRUE(x$converged)) "converged" else "NOT converged", "\n\n")
  cat(" Pairs per iteration:", s$pairs, "of", s$total_pairs, "| epoch:", s$epoch_length, "iterations | gain:", s$gain, "\n")
  print(x$phases, row.names = FALSE, digits = 3)
  cat("\n Last burn-in phase ended at a plateau:", as.logical(x$checks[["plateau"]]), "\n",
      "Skipped steps in the averaged phase:", x$checks[["skipped"]], "\n",
      "Rejected steps:", x$checks[["rejected"]], "\n",
      "Observed information positive definite:", as.logical(x$checks[["positive_definite"]]), "\n",
      "Largest predicted distance from exact PML:", signif(x$checks[["distance"]], 2), "standard errors\n")
  if (all(is.na(x$parameters$ratio))) {
    cat("\n No standard errors, noise or distances: use set_qn_control() or other starting values.\n")
    return(invisible(x))
  }
  cat("\n Optimisation noise relative to the standard error:\n")
  print(x$blocks, row.names = FALSE, digits = 2)
  if (!is.na(x$epochs_needed) && x$worst_ratio > x$target) {
    cat("\n average_epochs =", x$epochs_needed, "would bring the largest ratio to about", x$target, "\n")
  }
  invisible(x)
}

#' Plots of the stochastic approximation
#'
#' @param object A `lavaan.pl` object from [lavaan.pl()], fitted by stochastic approximation.
#' @param type `"trajectory"` for the estimates across iterations by parameter block, `"objective"` for the pairwise
#'   likelihood along the burn-in (it is not monitored in the averaged phase), `"noise"` for the optimisation noise
#'   relative to the standard error of each parameter.
#' @param blocks Parameter blocks to show, such as `"lambda"` or `"tau"`. If `NULL`, all of them.
#' @return The plotted data, invisibly. Without the tinyplot package, a warning and the data.
#' @export
plot_sa <- function(object, type = c("trajectory", "objective", "noise"), blocks = NULL) {
  fit <- sa_fit(object)
  type <- match.arg(type)
  path <- fit$stoFit$path
  block <- fit$blocks
  ends <- tapply(path$iter, path$phase, max)[c("1", "2")]
  ends <- ends[!is.na(ends)]
  long <- switch(type,
    trajectory = do.call(rbind, Map(function(iter, vector) {
      data.frame(iter = iter, estimate = transform_parameters(fit$model, vector),
                 id = seq_along(vector), block = block)
    }, path$iter, path$parameter_vector)),
    objective = data.frame(iter = path$iter, objective = path$nll)[path$phase < 3, ],
    noise = sa_parameters(fit)
  )
  if (!is.null(blocks) && type != "objective") long <- long[long$block %in% blocks, ]
  if (!requireNamespace("tinyplot", quietly = TRUE)) {
    warning("{tinyplot} must be installed to plot; the data are returned.", call. = FALSE)
    return(invisible(long))
  }
  phases <- bquote(abline(v = .(ends), lty = 2, col = "lightgrey"))
  switch(type,
    trajectory = do.call(tinyplot::tinyplot, list(
      estimate ~ iter | factor(id), data = long, facet = ~block, type = "l", legend = FALSE, col = "black",
      xlab = "Iterations", ylab = "Estimates", draw = phases
    )),
    objective = do.call(tinyplot::tinyplot, list(
      objective ~ iter, data = long, type = "b", pch = 19, xlab = "Iterations", ylab = "Negative pairwise log-likelihood per unit",
      draw = phases
    )),
    noise = do.call(tinyplot::tinyplot, list(
      ratio ~ seq_along(ratio), data = long, facet = ~block, type = "h", xlab = "Parameter",
      ylab = "Optimisation noise / standard error"
    ))
  )
  invisible(long)
}

#' Stability of the stochastic approximation across seeds
#'
#' Fits the same model and data with other seeds of the pair sampling, and compares the spread of the estimates with the
#' optimisation noise predicted from a single fit.
#'
#' @param object A `lavaan.pl` object from [lavaan.pl()], fitted by stochastic approximation.
#' @param seeds Seeds of the pair sampling.
#' @return An object of class `lavaanpl_sa_stability`: a list with the table `parameters` of each free parameter (`name`,
#'   `block`, standard error `se`, standard deviation across seeds `observed`, predicted noise `predicted`), the estimates
#'   across seeds `estimates` and the number of converged fits `converged`.
#' @export
sa_stability <- function(object, seeds = 1:10) {
  fit <- sa_fit(object)
  stopifnot(is.numeric(seeds), length(seeds) > 1, !anyDuplicated(seeds))
  fits <- lapply(seeds, function(seed) {
    control <- fit$control
    control$seed <- as.integer(seed)
    suppressWarnings(fit_lavaanpl(fit$data, fit$model, CONTROL = control, START = fit$init, NCORES = fit$cores))
  })
  pars <- sa_parameters(fit)
  estimates <- vapply(fits, function(x) transform_parameters(x$model, x$parameter_vector), numeric(length(fit$parameter_vector)))[fit$lincon$free_idx, , drop = FALSE]
  out <- list(
    parameters = data.frame(name = pars$name, block = pars$block, se = pars$se, observed = apply(estimates, 1, stats::sd),
                            predicted = pars$noise),
    estimates = estimates,
    converged = sum(vapply(fits, function(x) x$stoFit$converged, logical(1)))
  )
  class(out) <- "lavaanpl_sa_stability"
  out
}

#' @export
print.lavaanpl_sa_stability <- function(x, ...) {
  p <- x$parameters
  keep <- p$observed > 0 & p$predicted > 0
  cat("Stochastic approximation across", ncol(x$estimates), "seeds:", x$converged, "converged\n\n",
      "Standard deviation across seeds / standard error: median", signif(stats::median(p$observed / p$se), 2),
      ", largest", signif(max(p$observed / p$se), 2), "\n",
      "Predicted noise / observed: median", if (any(keep)) signif(stats::median((p$predicted / p$observed)[keep]), 2) else NA, "\n")
  invisible(x)
}

#' @export
plot.lavaanpl_sa_stability <- function(x, ...) {
  if (!requireNamespace("tinyplot", quietly = TRUE)) {
    warning("{tinyplot} must be installed to plot; the data are returned.", call. = FALSE)
    return(invisible(x$parameters))
  }
  do.call(tinyplot::tinyplot, list(
    observed ~ predicted | block, data = x$parameters, pch = 19, xlab = "Predicted optimisation noise",
    ylab = "Standard deviation across seeds", draw = quote(abline(0, 1, lty = 2, col = "lightgrey")), ...
  ))
  invisible(x$parameters)
}
