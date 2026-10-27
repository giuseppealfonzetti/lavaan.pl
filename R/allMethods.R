#' Show method for the output of fit_latpl()
#' @param x Object of class latpl_fit
#' @param ... Additional arguments
#' @export
print.latpl_fit <- function(x, ...){

  cat( "- Dimensions:\n",
       "  - Sample size:", x$dims$n, "\n",
       "  - Items:", x$dims$p,  "(", x$dims$pairs, " pairs)\n",
       "  - Latent traits:", x$dims$q,"\n\n",
       "- Free parameters:\n",
       "  - Thresholds:", x$dims$nthr, "\n",
       "  - Loadings:", x$dims$nload, "\n",
       "  - Latent correlations:", x$dims$ncorr, "\n",
       "  - Latent variances:", x$dims$nvar, "\n",
       "  - Regressions:", x$dims$nreg, "\n",
       "  - Residual covariances:", x$dims$nres, "\n",
       "  - Total:", x$dims$npar, "\n\n")

  if(x$method == 'ucminf'){

    cat(
      "- Numerical estimate obtained via ucminf.\n\n",
      "   Estimation time:", round(x$RTime, 2), "s\n",
      "   Data reduction:", round(x$data$freqTime,2), "s (in set_latpl_data())\n",
      "   Cores used:", x$cores, "\n",

      "\n- Use get_par() to extract parameter estimates."

    )
  }else{
    cat(
      "- Stochastic estimate\n\n",
      "    Pairs per iteration:", x$control$PAIRS_PER_ITERATION, " out of ", x$dims$pairs, "\n",
      "    Iterations:", x$stoFit$last_iter, "\n",
      "    Estimation time:", round(x$RTime, 2), "s\n",
      "    Data reduction:", round(x$data$freqTime,2), "s (in set_latpl_data())\n",
      "    Cores used:", x$cores, "\n",


      "\n- Use get_par() to extract parameter estimates."
    )
  }
  invisible(x)
}

#' Print method for `latpl_data` objects
#'
#' @param x A `latpl_data` object.
#' @param ... Ignored.
#'
#' @rdname set_latpl_data
#' @export
print.latpl_data <- function(x, ...) {
  cat(
    "<latpl_data>", x$n, "observations,", x$p, "items (", x$p * (x$p - 1) / 2, "pairs )\n",
    " Categories per item:", paste(range(x$cat), collapse = " to "), "\n",
    " Missing values:", sum(is.na(x$Y)), "(", x$missing, ")\n",
    " Sampling weights:", if (length(x$weights)) "yes" else "none", "\n",
    " Clusters:", if (length(x$cluster)) length(unique(x$cluster)) else "none", "\n",
    " Validation data:", if (is.null(x$valfreq)) "none" else "yes", "\n"
  )
  invisible(x)
}

#' Print method for `latpl_model` objects
#'
#' @param x A `latpl_model` object.
#' @param ... Ignored.
#'
#' @rdname set_latpl_model
#' @export
print.latpl_model <- function(x, ...) {
  cat(
    "<latpl_model>", x$p, "items,", x$q, "latent variables\n",
    " Free loadings:", sum(is.na(x$CONSTRMAT)), "\n",
    " Latent correlations:", if (x$q == 1) "none" else if (!x$CORRFLAG) "fixed to zero" else if (identical(x$LATCOV, "direct")) "direct" else "LKJ", "\n",
    " Latent variances:", if (anyNA(x$CONSTRLOGSD)) "free" else "fixed", "\n",
    " Free regressions:", sum(is.na(x$REGMAT)), "\n",
    " Free residual covariances:", sum(is.na(x$RESCOV[upper.tri(x$RESCOV)])), "\n",
    " Linear constraints:", NROW(x$CONSTRAINTS$A), "\n",
    " Identification:", if (x$STDLV) "standardised latent variables" else "marker loadings", "\n"
  )
  invisible(x)
}

#' Extract parameters
#'
#' @param OBJ Object of class latpl_fit.
#' @param ... Additional arguments
#' @export
get_par <- function(OBJ, ...) UseMethod('get_par')

#' Extract parameters
#'
#' @param OBJ Object of class latpl_fit.
#' @param OPTION Can take values \tabular{ll}{
#'    \code{'raw'} \tab The parameter vector. \cr
#'    \tab \cr
#'    \code{'transformed'} \tab Transformed vector of parameters, with interpretable latent correlations and standard deviations. \cr
#'     \tab \cr
#'    \code{'list'} \tab The model matrices, as returned by [theta2par()]. \cr
#' }
#' @param ... Additional arguments
#' @export
get_par.latpl_fit <- function(OBJ, OPTION = c("list", "raw", "transformed"), ...) {
  OPTION <- match.arg(OPTION)
  if (OPTION == "raw") return(OBJ$theta)
  par <- theta2par(OBJ$model, OBJ$theta)
  if (OPTION == "list") return(par)
  block <- theta_block(OBJ$model, OBJ$dims$nthr)
  R <- stats::cov2cor(par$latent_correlations)
  c(OBJ$theta[block %in% c("tau", "lambda")], if (any(block == "corr")) R[upper.tri(R)],
    exp(OBJ$theta[block == "logsd"]), OBJ$theta[block %in% c("beta", "rescov")])
}

#' Compute variance of the estimates
#' @param OBJ Object of class latpl_fit
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param VERBOSE TRUE if the function must print messages
#' @return A list with the Jacobian `trJacob` of the transformation to
#'   interpretable parameters, `H`, `J`, the sandwich `vcov` and its diagonal
#'   `asymptotic_variance` on the transformed scale, `optimisation_noise` (SA
#'   fits only), the optimisation noise matrix `noise` on the optimisation scale
#'   (SA fits only) and the timing `RTime`.
#' @export
compute_var <- function(OBJ, NUMDERIV = FALSE, VERBOSE = FALSE) UseMethod('compute_var')

#' @rdname compute_var
#' @export
compute_var.latpl_fit <- function(OBJ, NUMDERIV = FALSE, VERBOSE = FALSE) {
  start_time <- Sys.time()
  lc <- OBJ$lincon
  full <- function(x) as.vector(lc$K %*% x + lc$c)
  xf <- OBJ$theta[lc$free_idx]
  if (VERBOSE) message("- Computing H and J...")
  lat <- theta_block(OBJ$model, OBJ$dims$nthr) %in% c("corr", "logsd")
  trJacob <- diag(length(OBJ$theta))
  if (any(lat)) {
    trJacob[lat, lat] <- numDeriv::jacobian(function(x) {
      get_par(replace(OBJ, "theta", list(replace(OBJ$theta, lat, x))), OPTION = "transformed")[lat]
    }, OBJ$theta[lat])
  }
  trJacob <- trJacob %*% lc$K
  RcppParallel::setThreadOptions(numThreads = OBJ$cores)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  HJ <- cpp_HJ(OBJ$theta, OBJ$engine$model, OBJ$engine$data)
  H <- crossprod(lc$K, HJ$H %*% lc$K)
  J <- crossprod(lc$K, HJ$J %*% lc$K)
  if (NUMDERIV) {
    H <- numDeriv::jacobian(function(x) -as.vector(crossprod(lc$K, cpp_pl(full(x), OBJ$engine$model, OBJ$engine$data, 1L)$gradient)) / OBJ$dims$n, xf)
  }
  if (!(det(H) > 0)) stop("H not invertible.")
  invH <- solve(H)
  vcov <- trJacob %*% invH %*% J %*% invH %*% t(trJacob)
  opt_noise <- NA
  noise <- NULL
  if (OBJ$method == "SA") {
    P <- OBJ$dims$pairs
    ppi <- OBJ$control$PAIRS_PER_ITERATION
    a1 <- P * (P - ppi) / (ppi * (P - 1))
    a2 <- (P - ppi) / (ppi * (P - 1))
    noise <- invH %*% (a1 * H - a2 * J) %*% invH / (OBJ$dims$n * (OBJ$stoFit$last_iter - OBJ$stoFit$burnt))
    eig <- eigen((noise + t(noise)) / 2, symmetric = TRUE)
    noise <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))
    opt_noise <- diag(trJacob %*% noise %*% t(trJacob))
  }
  list(trJacob = trJacob, H = H, invH = trJacob %*% invH, J = J, vcov = vcov, asymptotic_variance = diag(vcov),
       optimisation_noise = opt_noise, noise = noise, RTime = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
}

plot_traj <- function(OBJ) UseMethod('plot_traj')
#' @export
plot_traj.latpl_fit <- function(OBJ){

  if(OBJ$method=="ucminf"){
    warning('Trajectories not available for numerical fit.')
    return(NULL)
  }

  if (!requireNamespace("tinyplot", quietly = TRUE)) {
    warning('{tinyplot} must be installed to plot trajectories.\n
            The trajectories list will be returned.')
    return(list(iters=OBJ$stoFit$path$iters, theta=OBJ$stoFit$path$avtheta))
  }

  type <- rep(c("thresholds", "loadings", "latent correlations", "latent variances"),
              c(OBJ$dims$nthr, OBJ$dims$nload, OBJ$dims$ncorr, OBJ$dims$nvar))
  long <- do.call(rbind, Map(function(iter, theta) {
    data.frame(iter = iter,
               par = get_par(replace(OBJ, "theta", list(theta)), OPTION = "transformed"),
               id = seq_along(theta),
               type = type)
  }, OBJ$stoFit$path$iters, OBJ$stoFit$path$avtheta))

  do.call(tinyplot::tinyplot, list(
    par ~ iter | factor(id), data = long, facet = ~type, type = "l",
    legend = FALSE, col = "black", xlab = "Iterations", ylab = "Estimates",
    draw = bquote(abline(v = .(OBJ$stoFit$burnt), lty = 2, col = "lightgrey"))
  ))
  invisible(long)
}

#' @export
plot_traj.lavaan <- function(OBJ) {
  plot_traj(OBJ@external$lavaan.pl)
}
