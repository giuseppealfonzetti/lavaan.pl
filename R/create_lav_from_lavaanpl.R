#' Convert a lavaan.pl fit into a lavaan.pl object
#'
#' @param fit Object of class `lavaanpl_internal` returned by [fit_lavaanpl()].
#' @param vars Output of [compute_var()] for `fit`, or `NULL` to skip standard errors.
#' @param template Unfitted [lavaan::lavaan] object for the same model and data, created with `do.fit = FALSE`.
#' @return An object of class `lavaan.pl`.
#' @keywords internal
create_lav_from_lavaanpl <- function(fit, vars = NULL, template) {
  frees <- lavaan_matrices(template, "free")
  par_lav <- function(vector) {
    pars <- vector_to_matrices(fit$model, vector)
    x <- numeric(template@Model@nx.free)
    for (g in seq_along(frees)) {
      for (m in intersect(names(frees[[g]]), names(pars[[g]]))) {
        free <- frees[[g]][[m]]
        x[free[free > 0]] <- pars[[g]][[m]][free > 0]
      }
    }
    x
  }

  # Change Model and implied slots ---------------------------------------------
  out <- new("lavaan.pl", template)
  x <- par_lav(fit$parameter_vector)
  logl <- cpp_pl(fit$parameter_vector, fit$cpp_inputs$model, fit$cpp_inputs$data)$pl
  out@Model <- lavaan::lav_model_set_parameters(template@Model, x)
  out@implied <- lavaan::lav_model_implied(out@Model)
  out@ParTable$est <- lavaan::lav_model_get_parameters(out@Model, type = "user", extra = TRUE)
  out@Options$do.fit <- TRUE
  out@Options$optim.method <- fit$method
  out@optim[c("x", "npar", "logl", "fx")] <- list(x, length(x), logl, -logl / fit$data$n)
  out@optim$fx.group <- out@optim$logl.group <- NA

  # Change optim slot (depends whether it is numFit or stoFit) -----------------
  if (fit$method != "SA") {
    out@optim$iterations <- as.integer(fit$numFit$info[["neval"]])
    out@optim$converged <- fit$numFit$convergence %in% c(1L, 2L) ||
      (fit$numFit$convergence == 4L && fit$numFit$info[["maxgradient"]] < template@Options$optim.dx.tol)
  } else {
    out@optim$iterations <- as.integer(fit$stoFit$iterations)
    out@optim$converged <- fit$stoFit$converged
  }
  if (isTRUE(out@Options$check.post)) lavaan::lavInspect(out, "post.check")

  # Change vcov slot -----------------------------------------------------------
  vcov <- NULL
  if (!is.null(vars)) {
    lat <- fit$blocks %in% c("corr", "logsd", "scale")
    delta <- matrix(0, length(x), length(fit$parameter_vector))
    delta[, !lat] <- numDeriv::jacobian(function(t) par_lav(replace(fit$parameter_vector, !lat, t)), fit$parameter_vector[!lat], method = "simple")
    if (any(lat)) delta[, lat] <- numDeriv::jacobian(function(t) par_lav(replace(fit$parameter_vector, lat, t)), fit$parameter_vector[lat])
    delta <- delta %*% fit$lincon$K
    hinv <- solve(vars$H)
    vcov <- delta %*% hinv %*% vars$J %*% hinv %*% t(delta) / fit$data$n
    attr(vcov, "E.inv") <- delta %*% hinv %*% t(delta)
    if (!is.null(vars$noise)) vcov <- vcov + delta %*% vars$noise %*% t(delta)
    out@vcov <- list(se = out@Options$se, information = out@Options$information[1], vcov = vcov)
    out@ParTable$se <- lavaan::lav_model_vcov_se(out@Model, out@ParTable, vcov = vcov)
  }

  # Change test slot -----------------------------------------------------------
  plrt <- NULL
  if (!identical(out@Options$test, "none") && !identical(fit$stoFit$check$pd, FALSE)) {
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
      attr(out_test, "info") <- list(
        ngroups = out@Data@ngroups, group.label = out@Data@group.label, information = out@Options$information,
        h1.information = out@Options$h1.information, observed.information = out@Options$observed.information
      )
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
  out@external <- list(lavaanpl_internal = fit, vars = vars, plrt = plrt)
  # Output ---------------------------------------------------------------------
  out@call <- template@call
  out
}
