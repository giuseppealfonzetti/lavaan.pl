#' Print method for `lavaanpl_internal` objects
#'
#' @param x Object of class `lavaanpl_internal` from [fit_lavaanpl()].
#' @param ... Ignored.
#' @keywords internal
#' @export
print.lavaanpl_internal <- function(x, ...){

  cat( "- Dimensions:\n",
       "  - Sample size:", x$data$n, "\n",
       "  - Items:", x$data$p,  "(", x$data$pairs, " pairs)\n",
       "  - Groups:", x$data$ngroups, "\n",
       "  - Latent traits:", x$dims$q,"\n\n",
       "- Free parameters:\n",
       "  - Thresholds:", x$dims$nthr, "\n",
       "  - Loadings:", x$dims$nload, "\n",
       "  - Latent correlations:", x$dims$ncorr, "\n",
       "  - Latent variances:", x$dims$nvar, "\n",
       "  - Regressions:", x$dims$nreg, "\n",
       "  - Residual covariances:", x$dims$nres, "\n",
       "  - Intercepts:", x$dims$nint, "\n",
       "  - Latent means:", x$dims$nmean, "\n",
       "  - Scales:", x$dims$nscale, "\n",
       "  - Effects of covariates:", x$dims$ngam + x$dims$nkap, "\n",
       "  - Total:", x$dims$npar, "\n\n")

  if(x$method != 'SA'){

    if(x$method == 'SA-QN') cat(
      "- Burn-in epochs: Adam", x$stoFit$adam_epochs, "| preconditioning", x$stoFit$precondition_epochs, "\n"
    )
    cat(
      "- Estimate from the", x$method, "method (ucminf).\n\n",
      "   Estimation time:", round(x$RTime, 2), "s\n",
      "   Data reduction:", round(x$data$freqTime,2), "s (in lavaanpl_data())\n",
      "   Cores used:", x$cores, "\n",

      "\n- Use vector_to_matrices() to extract the model matrices."

    )
  }else{
    cat(
      "- Estimate from the ", x$method, " method\n\n",
      "    Pairs per iteration:", x$stoFit$pairs_per_iteration, " out of ", x$data$pairs, "\n",
      "    Epochs: Adam", x$stoFit$adam_epochs, "| preconditioning", x$stoFit$precondition_epochs, "| averaged", x$control$average_epochs,
      "(", x$stoFit$epoch_length, "iterations each, gain:", x$stoFit$gain, ")\n",
      "    Rejected steps:", x$stoFit$rejections, "(", x$stoFit$rejections_averaging, "in the averaged phase )\n",
      "    Converged:", x$stoFit$converged, "\n",
      "    Largest distance from the optimum:", if (isTRUE(x$stoFit$check$pd)) signif(x$stoFit$check$distance, 2) else "not available", "standard errors\n",
      "    Estimation time:", round(x$RTime, 2), "s\n",
      "    Data reduction:", round(x$data$freqTime,2), "s (in lavaanpl_data())\n",
      "    Cores used:", x$cores, "\n",


      "\n- Use vector_to_matrices() to extract the model matrices."
    )
  }
  invisible(x)
}

#' Print method for `lavaanpl_data` objects
#'
#' @param x A `lavaanpl_data` object.
#' @param ... Ignored.
#'
#' @rdname lavaanpl_data
#' @export
print.lavaanpl_data <- function(x, ...) {
  cat(
    "<lavaanpl_data>", x$n, "observations,", x$p, "items (", x$pairs, "pairs )\n",
    " Numeric items:", sum(x$numeric), "\n",
    " Categories per ordinal item:", if (all(x$numeric)) "none" else paste(range(x$cat[!x$numeric]), collapse = " to "), "\n",
    " Missing values:", sum(is.na(x$Y)), "(", x$missing, ")\n",
    " Sampling weights:", if (length(x$weights)) "yes" else "none", "\n",
    " Clusters:", if (length(x$cluster)) length(unique(x$cluster)) else "none", "\n",
    " Groups:", x$ngroups, "\n",
    " Covariates:", ncol(x$X), "(", sum(vapply(x$strata, nrow, numeric(1))), "strata )\n"
  )
  invisible(x)
}

#' Print method for `lavaanpl_model` objects
#'
#' @param x A `lavaanpl_model` object.
#' @param ... Ignored.
#'
#' @rdname lavaanpl_model
#' @export
print.lavaanpl_model <- function(x, ...) {
  cat(
    "<lavaanpl_model>", x$p, "items,", x$q, "latent variables\n",
    " Free loadings:", sum(is.na(x$LAMBDA)), "\n",
    " Latent correlations:", if (x$q == 1) "none" else if (x$LATCOV == "none") "fixed to zero" else if (x$LATCOV == "direct") "direct" else "LKJ", "\n",
    " Latent variances:", if (anyNA(diag(x$PSI))) "free" else "fixed", "\n",
    " Free regressions:", sum(is.na(x$BETA)), "\n",
    " Free residual covariances:", sum(is.na(x$THETA[upper.tri(x$THETA)])), "\n",
    " Parameterization:", x$PARAMETERIZATION, "\n",
    " Free intercepts, latent means, scales:", sum(is.na(x$NU)), sum(is.na(x$ALPHA)), sum(is.na(item_scales(x))), "\n",
    " Free effects of covariates:", sum(is.na(x$GAMMA)) + sum(is.na(x$KAPPA)), "\n",
    " Linear constraints:", NROW(x$CONSTRAINTS$A), "\n"
  )
  invisible(x)
}

#' Variance of the estimates
#'
#' @param OBJ Object of class `lavaanpl_internal` from [fit_lavaanpl()].
#' @param NUMDERIV If `TRUE`, the Hessian by [numDeriv::jacobian()] of the gradient.
#' @return A list with the Jacobian `trJacob` of [transform_parameters()], `H`, `J`, the diagonal of the sandwich
#'   variance `asymptotic_variance` on the transformed scale, the Newton step to the optimum `step` and its size in
#'   standard errors `distance`, the estimates `parameter_vector` and the timing `RTime`; for SA fits, also the
#'   optimisation noise `optimisation_noise` on the transformed scale and `noise` on the optimisation scale.
#' @seealso [check_sa()], [compute_plrt()]
#' @keywords internal
compute_var <- function(OBJ, NUMDERIV = FALSE) {
  if (!NUMDERIV && identical(OBJ$vars$parameter_vector, OBJ$parameter_vector) && !is.null(OBJ$parameter_vector)) return(OBJ$vars)
  start_time <- Sys.time()
  lc <- OBJ$lincon
  full <- function(x) as.vector(lc$K %*% x + lc$c)
  xf <- OBJ$parameter_vector[lc$free_idx]
  lat <- OBJ$blocks %in% c("corr", "logsd", "scale")
  trJacob <- diag(length(OBJ$parameter_vector))
  if (any(lat)) {
    trJacob[lat, lat] <- numDeriv::jacobian(function(x) {
      transform_parameters(OBJ$model, replace(OBJ$parameter_vector, lat, x))[lat]
    }, OBJ$parameter_vector[lat])
  }
  trJacob <- trJacob %*% lc$K
  RcppParallel::setThreadOptions(numThreads = OBJ$cores)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  HJ <- cpp_sandwich(OBJ$parameter_vector, OBJ$cpp_inputs$model, OBJ$cpp_inputs$data, OBJ$method == "SA")
  H <- crossprod(lc$K, HJ$H %*% lc$K)
  J <- crossprod(lc$K, HJ$J %*% lc$K)
  if (NUMDERIV) {
    H <- numDeriv::jacobian(function(x) -as.vector(crossprod(lc$K, cpp_pl(full(x), OBJ$cpp_inputs$model, OBJ$cpp_inputs$data, 1L)$gradient)) / OBJ$data$n, xf)
  }
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R)) {
    stop(structure(class = c("lavaanpl_not_pd", "error", "condition"), list(message = "H not positive definite.", call = NULL)))
  }
  invH <- chol2inv(R)
  vcov <- trJacob %*% invH %*% J %*% invH %*% t(trJacob)
  opt_noise <- NA
  noise <- NULL
  if (OBJ$method == "SA") {
    P <- OBJ$data$pairs
    ppi <- OBJ$stoFit$pairs_per_iteration
    n <- OBJ$data$n
    G <- crossprod(lc$K, HJ$gradient) / n
    V <- (P - ppi) / (ppi * (P - 1)) * (P * crossprod(lc$K, HJ$pairvar %*% lc$K) / n^2 - tcrossprod(G))
    noise <- invH %*% V %*% invH / OBJ$stoFit$averaged
    opt_noise <- diag(trJacob %*% noise %*% t(trJacob))
  }
  step <- as.vector(invH %*% crossprod(lc$K, HJ$gradient)) / OBJ$data$n
  se <- sqrt(diag(vcov) / OBJ$data$n)
  list(trJacob = trJacob, H = H, J = J, asymptotic_variance = diag(vcov),
       optimisation_noise = opt_noise, noise = noise, step = step,
       distance = ifelse(se > 0, abs(as.vector(trJacob %*% step)) / se, 0), parameter_vector = OBJ$parameter_vector,
       RTime = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
}
