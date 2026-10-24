#' Show method for the output of fit_plFA()
#' @param object Object of class StoFit
setMethod('show', 'PlFaFit', function(object){

  cat( "- Dimensions:\n",
       "  - Sample size:", object@dims@n, "\n",
       "  - Items:", object@dims@p,  "(", object@dims@pairs, " pairs)\n",
       "  - Latent traits:", object@dims@q,"\n\n",
       "- Free parameters:\n",
       "  - Thresholds:", object@dims@nthr, "\n",
       "  - Loadings:", object@dims@nload, "\n",
       "  - Latent correlations:", object@dims@ncorr, "\n",
       "  - Latent variances:", object@dims@nvar, "\n",
       "  - Total:", object@dims@npar, "\n\n")

  if(object@method == 'ucminf'){

    cat(
      "- Numerical estimate obtained via ucminf.\n\n",
      "   Total time:", round(object@RTime, 2), "s (Data reduction:", round(object@freqTime,2), "s)\n",
      "   Cores used:", object@cores, "\n",

      "\n- Use getPar() to extract parameter estimates."

    )
  }else{
    cat(
      "- Stochastic estimate\n\n",
      "    Pairs per iteration:", object@stoFit@control$PAIRS_PER_ITERATION, " out of ", object@dims@pairs, "\n",
      "    Iterations:", object@stoFit@last_iter, "\n",
      "    Total time:", round(object@RTime, 2), "s (Data reduction:", round(object@freqTime,2), "s)\n",
      "    Cores used:", object@cores, "\n",


      "\n- Use getPar() to extract parameter estimates."
    )
  }
})

#' Extract parameters
#'
#' @param OBJ Object of class StoFit, PlFaFit. Also usable on raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
#' @export
setGeneric('getPar', function(OBJ, ...) standardGeneric('getPar'))
#' Extract parameters
#'
#' @param OBJ Object of class PlFaFit.
#' @param ... Additional arguments
setMethod('getPar', 'PlFaFit', function(OBJ, ...) extract_par(THETA       = OBJ@theta,
                                                              CONSTRMAT   = OBJ@cnstr@loadings,
                                                              CONSTRLOGSD = OBJ@cnstr@loglatsd,
                                                              LLC         = OBJ@cnstr@llc,
                                                              NTHR        = OBJ@dims@nthr,
                                                              NLOAD       = OBJ@dims@nload,
                                                              NCORR       = OBJ@dims@ncorr,
                                                              NVAR        = OBJ@dims@nvar,
                                                              ...))
#' Extract parameters
#'
#' @param OBJ Raw theta vector to extract \code{'transformed'} and \code{'list'}.
#' @param ... Additional arguments
setMethod('getPar', 'vector', function(OBJ, ...) extract_par(THETA = OBJ, ...))

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
#' @param OBJ Object of class plFaFit
#' @param DATA Original data
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param OPTION `transformed` if correlations are of interest. `raw` for inference on the reparametrised level.
#' @param VERBOSE TRUE if the function must print messages
#' @export
setGeneric('computeVar', function(OBJ, DATA, NUMDERIV = F, OPTION = 'transformed', VERBOSE = FALSE) standardGeneric('computeVar'), signature = 'OBJ')
#' Compute variance of the estimates
#' @param OBJ Object of class plFaFit
#' @param DATA Original data
#' @param NUMDERIV TRUE if the hessian must be computed using \link[numDeriv]{jacobian}
#' @param OPTION `transformed` if correlations are of interest. `raw` for inference on the reparametrised level.
#' @param VERBOSE TRUE if the function must print messages
setMethod('computeVar', 'PlFaFit',
          function(OBJ, DATA, NUMDERIV, OPTION, VERBOSE){ if(OBJ@method == 'ucminf'){
            compute_var(THETA      = OBJ@theta,
                       C_VEC       = OBJ@dims@cat,
                       N           = OBJ@dims@n,
                       IT          = NULL,
                       PAIRS       = NULL,
                       PPI         = NULL,
                       CONSTRMAT   = OBJ@cnstr@loadings,
                       CONSTRLOGSD = OBJ@cnstr@loglatsd,
                       LLC         = OBJ@cnstr@llc,
                       NTHR        = OBJ@dims@nthr,
                       NLOAD       = OBJ@dims@nload,
                       NCORR       = OBJ@dims@ncorr,
                       NVAR        = OBJ@dims@nvar,
                       FREQ        = OBJ@freq,
                       DATA        = DATA,
                       METHOD      = OBJ@method,
                       NUMDERIV    = NUMDERIV,
                       NCORES      = OBJ@cores,
                       VERBOSE     = VERBOSE)
          }else{
            compute_var(THETA      = OBJ@theta,
                       C_VEC       = OBJ@dims@cat,
                       N           = OBJ@dims@n,
                       IT          = OBJ@stoFit@last_iter - OBJ@stoFit@control$BURN,
                       PAIRS       = OBJ@dims@pairs,
                       PPI         = OBJ@stoFit@control$PAIRS_PER_ITERATION,
                       CONSTRMAT   = OBJ@cnstr@loadings,
                       CONSTRLOGSD = OBJ@cnstr@loglatsd,
                       LLC         = OBJ@cnstr@llc,
                       NTHR        = OBJ@dims@nthr,
                       NLOAD       = OBJ@dims@nload,
                       NCORR       = OBJ@dims@ncorr,
                       NVAR        = OBJ@dims@nvar,
                       FREQ        = OBJ@freq,
                       DATA        = DATA,
                       METHOD      = OBJ@method,
                       NUMDERIV    = NUMDERIV,
                       NCORES      = OBJ@cores,
                       VERBOSE     = VERBOSE)
          }}
)

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

get_AIC <- function(NLL, INVH, J){
  return(NLL+sum(diag(J%*%INVH)))
}

get_BIC <- function(NLL, INVH, J, N){
  return(2*NLL+sum(diag(J%*%INVH))*log(N))
}


#' @importFrom rlang .data
setGeneric('plotTraj', function(OBJ) standardGeneric('plotTraj'))
setMethod('plotTraj', 'PlFaFit', function(OBJ){

  if(OBJ@method=="ucminf"){
    warning('Trajectories not available for numerical fit.')
    return(NULL)
  }

  if (!requireNamespace(c("dplyr", "tidyr", "purrr", "ggplot2"), quietly = TRUE)) {
    warning('{dplyr}, {tidyr}, {purrr}, {ggplot2} must be installed to plot trajectories.\n
            The trajectories list will be returned.')
    return(list(iters=OBJ@stoFit@path$iters, theta=OBJ@stoFit@path$avtheta))
  }

  dplyr::tibble(iter=OBJ@stoFit@path$iters, theta=OBJ@stoFit@path$avtheta) |>
    dplyr::mutate(theta=purrr::map(.data$theta, ~{dplyr::tibble( par=extract_par(
      THETA=.x, OPTION = c('transformed'),
      CONSTRMAT = OBJ@cnstr@loadings,
      CONSTRLOGSD=OBJ@cnstr@loglatsd,
      LLC=OBJ@cnstr@llc,
      NTHR=OBJ@dims@nthr,
      NLOAD=OBJ@dims@nload,
      NCORR=OBJ@dims@ncorr,
      NVAR=OBJ@dims@nvar),
      id=1:length(.x),
      type=c(rep("thresholds", OBJ@dims@nthr),
             rep("loadings", OBJ@dims@nload),
             rep("latent correlations", OBJ@dims@ncorr),
             rep("latent variances", OBJ@dims@nvar)))})) |>
    tidyr::unnest("theta") |>
    ggplot2::ggplot(ggplot2::aes(x=.data$iter, y = .data$par, group = .data$id))+
    ggplot2::facet_wrap(~.data$type) +
    ggplot2::geom_vline(xintercept = OBJ@stoFit@burnt, linetype="dashed", col="lightgrey")+
    ggplot2::geom_line() +
    ggplot2::labs(x="Iterations", y="Estimates")+
    ggplot2::theme_bw()

})

