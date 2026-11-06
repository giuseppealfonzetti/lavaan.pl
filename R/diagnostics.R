#' Free parameters of a fit, with their precision
#'
#' @param object A `lavaan.pl` object from [lavaan.pl()].
#' @return A data frame with one row per free parameter: `name`, `block`, `estimate`, standard error `se` and the
#'   distance from the optimum in standard errors `distance`, and with stochastic approximation the standard deviation of
#'   the optimisation noise `noise` and its ratio to the standard error `ratio`. Without a positive definite observed
#'   information, all but the estimates are `NA`.
#' @seealso [diagnostics()]
#' @keywords internal
fit_parameters <- function(object) {
  fit <- get_lavaanpl_internal(object)
  vars <- object@external$vars
  if (is.null(vars)) {
    vars <- tryCatch(compute_var(fit), lavaanpl_not_pd = function(e) NULL)
  }
  out <- data.frame(
    name = parameter_names(fit$data, fit$model),
    block = fit$blocks,
    estimate = transform_parameters(fit$model, fit$parameter_vector),
    se = NA,
    distance = NA
  )
  if (fit$method == "SA") {
    out[c("noise", "ratio")] <- NA
  }
  if (!is.null(vars)) {
    out$se <- sqrt(vars$asymptotic_variance / fit$data$n)
    out$distance <- vars$distance
    if (fit$method == "SA") {
      out$noise <- sqrt(pmax(vars$optimisation_noise, 0))
      out$ratio <- ifelse(out$se > 0, out$noise / out$se, 0)
    }
  }
  out[fit$lincon$free_idx, ]
}

#' Diagnostics of the optimisation
#'
#' Reports how the optimiser behaved and whether its estimate is at the optimum of the pairwise likelihood, for every
#' method. The distance is how far the estimate is from that optimum, in standard errors, predicted from the gradient. With
#' stochastic approximation, the optimisation noise is the variability due to the pairs sampled, on top of the
#' statistical variability.
#'
#' @param object A `lavaan.pl` object from [lavaan.pl()].
#' @return An object of class `lavaanpl_diagnostics`: a list with the `method`, `converged` as in
#'   `lavInspect(object, "converged")`, the run time `seconds`, the table of the phases `phases` (epochs and seconds; for
#'   quasi-Newton, each evaluation of the pairwise likelihood counts as one epoch), the table of the checks `checks`
#'   (`check`, `value`, and `ok`, `NA` when there is no rule), and the table of the free parameters `parameters`.
#' @seealso [plot.lavaan.pl()]
#' @export
diagnostics <- function(object) {
  fit <- get_lavaanpl_internal(object)
  pars <- fit_parameters(object)
  st <- fit$stoFit
  qn <- fit$numFit
  sa <- fit$method != "QN"
  num <- fit$method != "SA"
  burn <- if (sa) unname(st$time[c("adam", "precondition")]) else NULL
  last <- if (num) {
    fit$RTime - sum(burn, if (sa) st$time[["monitor"]])
  } else {
    st$time[["averaging"]]
  }
  phases <- data.frame(
    phase = c(
      if (sa) c("Adam", "Preconditioning"),
      if (num) "QN" else "Averaged"
    ),
    epochs = c(
      if (sa) c(st$adam_epochs, st$precondition_epochs),
      if (num) qn$info[["neval"]] else fit$control$average_epochs
    ),
    seconds = c(burn, last)
  )
  pd <- !all(is.na(pars$se))
  worst <- function(x) {
    if (pd) {
      paste0(format(signif(max(x), 2)), " (", pars$name[which.max(x)], ")")
    } else {
      "not available"
    }
  }
  check <- function(name, value, ok = NA) {
    data.frame(check = name, value = as.character(value), ok = ok)
  }
  checks <- rbind(
    check("Observed information positive definite", pd, pd),
    check(
      "Largest coordinate one-Newton-step/standard error",
      worst(pars$distance),
      if (pd && !num) max(pars$distance) <= max(1, 4 * max(pars$ratio)) else NA
    )
  )
  if (!num) {
    checks <- rbind(
      checks,
      check("Largest optimisation noise / standard error", worst(pars$ratio)),
      check("Burn-in ended at a plateau", st$plateau),
      check(
        "Skipped steps in the averaged phase",
        st$skipped_averaging,
        st$skipped_averaging == 0
      ),
      check("Rejected steps", st$rejections)
    )
  } else {
    checks <- rbind(
      checks,
      check("QN", qn$message),
      check(
        "QN evaluations",
        paste(qn$info[["neval"]], "out of", fit$control$maxeval),
        qn$info[["neval"]] < fit$control$maxeval
      ),
      check(
        "Largest gradient per observation",
        format(signif(qn$info[["maxgradient"]], 2)),
        qn$info[["maxgradient"]] < 1e-3
      )
    )
    if (length(st$start_logl) > 1) {
      spread <- diff(range(st$start_logl))
      checks <- rbind(
        checks,
        check(
          "Pairwise log-likelihood spread across starts",
          format(signif(spread, 2)),
          spread <= 1e-3
        )
      )
    }
  }
  out <- list(
    method = fit$method,
    converged = object@optim$converged,
    seconds = fit$RTime,
    phases = phases,
    checks = checks,
    parameters = pars
  )
  class(out) <- "lavaanpl_diagnostics"
  out
}

#' @export
print.lavaanpl_diagnostics <- function(x, ...) {
  cat(
    x$method,
    "fit:",
    if (isTRUE(x$converged)) "converged" else "NOT converged",
    paste0("(", signif(x$seconds, 3), " s)"),
    "\n\n"
  )
  print(x$phases, row.names = FALSE, digits = 3)
  cat(
    "\n",
    paste0(
      ifelse(x$checks$ok %in% FALSE, "! ", "  "),
      x$checks$check,
      ": ",
      x$checks$value,
      "\n"
    ),
    sep = ""
  )
  invisible(x)
}
