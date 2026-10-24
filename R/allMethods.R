#' Show method for the output of fit_plFA()
#' @param x Object of class PlFaFit
#' @param ... Additional arguments
#' @export
print.PlFaFit <- function(x, ...){

  cat( "- Dimensions:\n",
       "  - Sample size:", x$dims$n, "\n",
       "  - Items:", x$dims$p,  "(", x$dims$pairs, " pairs)\n",
       "  - Latent traits:", x$dims$q,"\n\n",
       "- Free parameters:\n",
       "  - Thresholds:", x$dims$nthr, "\n",
       "  - Loadings:", x$dims$nload, "\n",
       "  - Latent correlations:", x$dims$ncorr, "\n",
       "  - Latent variances:", x$dims$nvar, "\n",
       "  - Total:", x$dims$npar, "\n\n")

  if(x$method == 'ucminf'){

    cat(
      "- Numerical estimate obtained via ucminf.\n\n",
      "   Estimation time:", round(x$RTime, 2), "s\n",
      "   Data reduction:", round(x$data$freqTime,2), "s (in set_plfa_data())\n",
      "   Cores used:", x$cores, "\n",

      "\n- Use getPar() to extract parameter estimates."

    )
  }else{
    cat(
      "- Stochastic estimate\n\n",
      "    Pairs per iteration:", x$control$PAIRS_PER_ITERATION, " out of ", x$dims$pairs, "\n",
      "    Iterations:", x$stoFit$last_iter, "\n",
      "    Estimation time:", round(x$RTime, 2), "s\n",
      "    Data reduction:", round(x$data$freqTime,2), "s (in set_plfa_data())\n",
      "    Cores used:", x$cores, "\n",


      "\n- Use getPar() to extract parameter estimates."
    )
  }
  invisible(x)
}

#' Print method for `plfa_data` objects
#'
#' @param x A `plfa_data` object.
#' @param ... Ignored.
#'
#' @rdname set_plfa_data
#' @export
print.plfa_data <- function(x, ...) {
  cat(
    "<plfa_data>", x$n, "observations,", x$p, "items (", x$p * (x$p - 1) / 2, "pairs )\n",
    " Categories per item:", paste(range(x$cat), collapse = " to "), "\n",
    " Validation data:", if (is.null(x$valfreq)) "none" else "yes", "\n"
  )
  invisible(x)
}

#' Print method for `plfa_model` objects
#'
#' @param x A `plfa_model` object.
#' @param ... Ignored.
#'
#' @rdname set_plfa_model
#' @export
print.plfa_model <- function(x, ...) {
  cat(
    "<plfa_model>", x$p, "items,", x$q, "latent variables\n",
    " Free loadings:", sum(is.na(x$CONSTRMAT)), "\n",
    " Latent correlations:", if (x$q == 1) "none" else if (x$CORRFLAG) "free" else "fixed to zero", "\n",
    " Latent variances:", if (anyNA(x$CONSTRLOGSD)) "free" else "fixed", "\n",
    " Linear loading constraints:", length(x$LLC), "\n",
    " Identification:", if (x$STDLV) "standardised latent variables" else "marker loadings", "\n"
  )
  invisible(x)
}

#' Extract parameters
#'
#' @param OBJ Object of class PlFaFit. Also usable on raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
#' @export
getPar <- function(OBJ, ...) UseMethod('getPar')
#' Extract parameters
#'
#' @param OBJ Object of class PlFaFit.
#' @param ... Additional arguments
#' @export
getPar.PlFaFit <- function(OBJ, ...) extract_par(THETA       = OBJ$theta,
                                                              CONSTRMAT   = OBJ$model$CONSTRMAT,
                                                              CONSTRLOGSD = OBJ$model$CONSTRLOGSD,
                                                              LLC         = OBJ$model$LLC,
                                                              NTHR        = OBJ$dims$nthr,
                                                              NLOAD       = OBJ$dims$nload,
                                                              NCORR       = OBJ$dims$ncorr,
                                                              NVAR        = OBJ$dims$nvar,
                                                              ...)
#' Extract parameters
#'
#' @param OBJ Raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
#' @export
getPar.default <- function(OBJ, ...) extract_par(THETA = OBJ, ...)

#' Extract parameters
#'
#' @param THETA Unconstrained parameter vector.
#' @param OPTION Can take values \tabular{ll}{
#'    \code{'raw'} \tab Returns \code{THETA}. \cr
#'    \tab \cr
#'    \code{'transformed'} \tab Transformed vector of parameters, with interpretable latent correlations. \cr
#'     \tab \cr
#'    \code{'list'} \tab List containing the ordered thresholds vector, the loading matrix and the latent correlation matrix. \cr
#' }
#' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA` refers to free loading parameters. Elements set to numerical values denote fixed values constraints.
#' @param CONSTRLOGSD \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
#' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
#' @param NTHR Number of thresholds parameters.
#' @param NLOAD Number of free loadings parameters
#' @param NCORR Number of free latent correlations parameters.
#' @param NVAR Number of free latent variance parameters.
extract_par <- function(THETA, OPTION = c('list', 'raw', 'transformed'),
                        CONSTRMAT, CONSTRLOGSD, LLC, NTHR, NLOAD, NCORR, NVAR){
  OPTION <- match.arg(OPTION)

  q <- ncol(CONSTRMAT)
  if(OPTION == 'raw'){
    return(THETA)
    }else if(OPTION == 'transformed'){

      vec <- THETA[1:(NTHR+NLOAD)]
      if(NCORR>0){
        vec <- c(vec, get_corr(THETA = THETA,
                               NTHR  = NTHR,
                               NLOAD = NLOAD,
                               NCORR = NCORR,
                               NVAR  = NVAR,
                               Q     = q))
      }
      if(NVAR>0){
        vec <- c(vec, exp(THETA[(NTHR+NLOAD+NCORR+1):(NTHR+NLOAD+NCORR+NVAR)]))
      }

      return(vec)


    }else if(OPTION == 'list') {
    out <- list()
    thr = THETA[1:NTHR]
    ld <- cpp_loadings_theta2mat(
      THETA     = THETA,
      CONSTRMAT = CONSTRMAT,
      LLC       = LLC,
      NTHR      = NTHR,
      NLOAD     = NLOAD)

    S <- cpp_latvar_theta2mat(
      THETA       = THETA,
      CONSTRLOGSD = CONSTRLOGSD,
      NTHR        = NTHR,
      NLOAD       = NLOAD,
      NCORR       = NCORR,
      NVAR        = NVAR,
      Q           = q
    )

    out <- list(thresholds = thr, loadings = ld, latent_correlations = S)
    return(out)
  }
}

#' Compute variance of the estimates
#' @param OBJ Object of class PlFaFit
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param VERBOSE TRUE if the function must print messages
#' @return A list with the Jacobian `trJacob` of the transformation to
#'   interpretable parameters, `H`, `J`, the sandwich `vcov` and its diagonal
#'   `asymptotic_variance` on the transformed scale, `optimisation_noise` (SA
#'   fits only) and the timing `RTime`.
#' @export
computeVar <- function(OBJ, NUMDERIV = F, VERBOSE = FALSE) UseMethod('computeVar')
#' Compute variance of the estimates
#' @param OBJ Object of class PlFaFit
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param VERBOSE TRUE if the function must print messages
#' @return A list with the Jacobian `trJacob` of the transformation to
#'   interpretable parameters, `H`, `J`, the sandwich `vcov` and its diagonal
#'   `asymptotic_variance` on the transformed scale, `optimisation_noise` (SA
#'   fits only) and the timing `RTime`.
#' @export
computeVar.PlFaFit <- function(OBJ, NUMDERIV = F, VERBOSE = FALSE){ if(OBJ$method == 'ucminf'){
            compute_var(THETA      = OBJ$theta,
                       C_VEC       = OBJ$dims$cat,
                       N           = OBJ$dims$n,
                       IT          = NULL,
                       PAIRS       = NULL,
                       PPI         = NULL,
                       CONSTRMAT   = OBJ$model$CONSTRMAT,
                       CONSTRLOGSD = OBJ$model$CONSTRLOGSD,
                       LLC         = OBJ$model$LLC,
                       NTHR        = OBJ$dims$nthr,
                       NLOAD       = OBJ$dims$nload,
                       NCORR       = OBJ$dims$ncorr,
                       NVAR        = OBJ$dims$nvar,
                       FREQ        = OBJ$data$freq,
                       DATA        = OBJ$data$Y,
                       METHOD      = OBJ$method,
                       NUMDERIV    = NUMDERIV,
                       NCORES      = OBJ$cores,
                       VERBOSE     = VERBOSE)
          }else{
            compute_var(THETA      = OBJ$theta,
                       C_VEC       = OBJ$dims$cat,
                       N           = OBJ$dims$n,
                       IT          = OBJ$stoFit$last_iter - OBJ$control$BURN,
                       PAIRS       = OBJ$dims$pairs,
                       PPI         = OBJ$control$PAIRS_PER_ITERATION,
                       CONSTRMAT   = OBJ$model$CONSTRMAT,
                       CONSTRLOGSD = OBJ$model$CONSTRLOGSD,
                       LLC         = OBJ$model$LLC,
                       NTHR        = OBJ$dims$nthr,
                       NLOAD       = OBJ$dims$nload,
                       NCORR       = OBJ$dims$ncorr,
                       NVAR        = OBJ$dims$nvar,
                       FREQ        = OBJ$data$freq,
                       DATA        = OBJ$data$Y,
                       METHOD      = OBJ$method,
                       NUMDERIV    = NUMDERIV,
                       NCORES      = OBJ$cores,
                       VERBOSE     = VERBOSE)
          }}

compute_var <- function(THETA, C_VEC, N, IT = NULL, PAIRS = NULL, PPI = NULL,
                        CONSTRMAT, CONSTRLOGSD, LLC, NTHR, NLOAD, NCORR, NVAR,
                        FREQ, DATA, METHOD, NUMDERIV = F, NCORES=1,
                        VERBOSE = FALSE){

  start_time <- Sys.time()



  # Get jacobian of reparameterisation
  if (isTRUE(VERBOSE)) message('- Computing Delta...')
  Rwr_getPar <- function(x){
    getPar(OBJ         = x,
           OPTION      = "transformed",
           CONSTRMAT   = CONSTRMAT,
           CONSTRLOGSD = CONSTRLOGSD,
           LLC         = LLC,
           NTHR        = NTHR,
           NLOAD       = NLOAD,
           NCORR       = NCORR,
           NVAR        = NVAR)
  }
  trJacob <- numDeriv::jacobian(Rwr_getPar, x=THETA)
  end_time <- Sys.time()
  totaltime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
  if (isTRUE(VERBOSE)) message('Done! (', round(totaltime, 2),' secs)')

  # Sample estimation of of J (and H)
  if (isTRUE(VERBOSE)) message('- Estimating H and J...')
  RcppParallel::setThreadOptions(numThreads = NCORES)
  HJ <- cpp_sample_estimators_HJ(
    THETA = THETA,
    FREQ = FREQ,
    DATA = DATA,
    C_VEC = C_VEC,
    CONSTRMAT = CONSTRMAT,
    CONSTRLOGSD = CONSTRLOGSD,
    LLC = LLC,
    N = nrow(DATA),
    CORRFLAG = NCORR>0,
    NTHR = NTHR,
    NLOAD = NLOAD,
    NCORR= NCORR,
    NVAR = NVAR
  )
  Jhat <- HJ$J

  end_time <- Sys.time()
  totaltime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
  if (isTRUE(VERBOSE)) message('Done! (', round(totaltime, 2),' secs)')

  # Compute H if needed
  opt_noise <- NA
  Hhat <- HJ$H
  if(NUMDERIV){
    if (isTRUE(VERBOSE)) message('- Computing H numerically...')
    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = N,
        C_VEC      = C_VEC,
        CONSTRMAT  = CONSTRMAT,
        CONSTRLOGSD= CONSTRLOGSD,
        LLC        = LLC,
        FREQ       = FREQ,
        THETA      = par_vec,
        CORRFLAG   = (NCORR>0),
        NTHR       = NTHR,
        NLOAD      = NLOAD,
        NCORR      = NCORR,
        NVAR       = NVAR,
        GRFLAG     = 1,
        SILENTFLAG = 1
      )
      mod$iter_ngradient/N
    }

    Hhat <- numDeriv::jacobian(Rwr_ngr, THETA)
    end_time <- Sys.time()
    totaltime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
    if (isTRUE(VERBOSE)) message('Done! (', round(totaltime, 2),' secs)')
  }

  if (isTRUE(VERBOSE)) message('- Inverting H...')
  if(!(det(Hhat)>0)) stop('H not invertible.')
  invH <- solve(Hhat)


  vcov <- trJacob %*% invH %*% Jhat %*% invH %*% t(trJacob)

  if(METHOD =='ucminf'){
    asy_var <- diag(vcov)
  } else {
    asy_var <- diag(vcov)
    a1 <- PAIRS*(PAIRS - PPI)/(PPI*(PAIRS-1))
    a2 <- (PAIRS-PPI)/(PPI*(PAIRS-1))
    opt_noise <- diag( trJacob %*% (invH%*%(a1*Hhat - a2*Jhat)%*%invH/(N*IT) ) %*% t(trJacob) )
  }

  end_time <- Sys.time()
  totaltime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
  if (isTRUE(VERBOSE)) message('Done! (', round(totaltime, 2),' secs)')

  return(
    list(
      trJacob = trJacob,
      H = Hhat,
      invH = trJacob%*%invH,
      J = Jhat,
      vcov = vcov,
      asymptotic_variance = asy_var,
      optimisation_noise = opt_noise,
      RTime = totaltime
    )
  )

}

plotTraj <- function(OBJ) UseMethod('plotTraj')
#' @export
plotTraj.PlFaFit <- function(OBJ){

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
               par = getPar(replace(OBJ, "theta", list(theta)), OPTION = "transformed"),
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
plotTraj.lavaan <- function(OBJ) {
  plotTraj(OBJ@external$lavaan.pl)
}
