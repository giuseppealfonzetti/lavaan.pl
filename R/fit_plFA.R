#'Fit factor models for ordinal data with pairwise likelihood methods
#'
#'@description
#'
#'Fit latent variable models for ordinal variable using pairwise likelihood
#'methods running via quasi-Newton BFGS or stochastic approximations.
#'
#'@param DATA Data object from [set_plfa_data()].
#'@param MODEL Model object from [set_plfa_model()].
#'@param CONTROL Control object of the main optimiser: [set_ucminf_control()]
#'  for the quasi-Newton BFGS optimiser from the \code{ucminf} package,
#'  [set_sa_control()] for Stochastic Approximations.
#'@param INIT Starting values object from [set_init()].
#'@param NCORES Integer value setting the number of threads to carry out the
#'  estimation.
#'@param VERBOSE `TRUE` for verbose output
#'
#'@return Object of class `PlFaFit`: a list with the inputs `data`, `model` and
#'  `control`, the starting values `init`, the estimates `theta` (use
#'  [getPar()]), the dimensions `dims`, the optimiser output `numFit` (ucminf)
#'  or `stoFit` (SA), and the timing `RTime`.
#'@export
fit_plFA <- function(
    DATA,
    MODEL,
    CONTROL = set_ucminf_control(),
    INIT = set_init(),
    NCORES = 1,
    VERBOSE = FALSE
){

  start_time <- Sys.time()

  stopifnot(
    inherits(DATA, "plfa_data"),
    inherits(MODEL, "plfa_model"),
    inherits(CONTROL, c("plfa_ucminf_control", "plfa_sa_control")),
    inherits(INIT, "plfa_init"),
    DATA$p == MODEL$p
  )
  METHOD <- if (inherits(CONTROL, "plfa_sa_control")) "SA" else "ucminf"

  # Identify model dimensions
  dims <- check_dims(DATA$Y, MODEL)
  if(METHOD=="SA" && dims$p<40) warning("SA defaults are tuned for models with at least 40 items: check the estimates against set_ucminf_control().", call. = FALSE)
  INIT_METHOD <- check_init_method(dims, INIT$METHOD, INIT$THETA)

  tmp <- structure(
    list(data = DATA,
         model = MODEL,
         control = CONTROL,
         dims = list(n     = dims$n,
                     p     = dims$p,
                     q     = dims$q,
                     cat   = dims$cat,
                     pairs = dims$pairs,
                     nthr  = dims$nthr,
                     nload = dims$nload,
                     ncorr = dims$ncorr,
                     nvar  = dims$nvar,
                     npar  = dims$d),
         method = METHOD),
    class = "PlFaFit")

  # Version
  tmp$version <- as.character(utils::packageVersion("lavaan.pl"))

  # Set up multi-threads computations
  if (is.null(NCORES)) NCORES <- 1L
  tmp$cores <- NCORES

  RcppParallel::setThreadOptions(numThreads = NCORES)

  # Check Initialisation
  tmp$init <- check_init(INIT$THETA, DATA$freq, dims, MODEL, INIT_METHOD, VERBOSE, INIT_CONTROL = INIT$CONTROL, UFREQ = DATA$ufreq)

  # Numerical optimisation
  if(METHOD == 'ucminf'){

    if(VERBOSE) message('- Optimising with ucminf...')
    tmp$method <- "ucminf"

    # Compute frequency table bivariate patterns

    Rwr_ncl <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = MODEL$CONSTRMAT,
        CONSTRLOGSD= MODEL$CONSTRLOGSD,
        LLC        = MODEL$LLC,
        FREQ       = DATA$freq,
        THETA      = par_vec,
        CORRFLAG   = MODEL$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        SILENTFLAG = 1,
        UFREQ      = DATA$ufreq
      )
      out <- mod$iter_nll/dims$n
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = MODEL$CONSTRMAT,
        CONSTRLOGSD= MODEL$CONSTRLOGSD,
        LLC        = MODEL$LLC,
        FREQ       = DATA$freq,
        THETA      = par_vec,
        CORRFLAG   = MODEL$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        GRFLAG     = 1,
        SILENTFLAG = 1,
        UFREQ      = DATA$ufreq
      )

      out <- mod$iter_ngradient/dims$n
      return(out)
    }

    # list of ucminf args
    args <- list(
      'par' = tmp$init,
      'fn' = Rwr_ncl,
      'gr' = Rwr_ngr,
      'control' = list(grtol = CONTROL$GRTOL, xtol = CONTROL$XTOL, stepmax = CONTROL$STEPMAX, maxeval = CONTROL$MAXEVAL),
      'hessian' = 0)

    # optimisation
    opt <- do.call(ucminf::ucminf, args)


    end_time <- Sys.time()
    tmp$numFit <- opt
    tmp$theta <- opt$par
    tmp$RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])


    if(VERBOSE) message('Done! (', round(tmp$RTime,2),' secs)')
    RcppParallel::setThreadOptions(numThreads = 1)
    return(tmp)
  }

  if(METHOD=="SA"){

    if(VERBOSE) message('- Optimising via pairwise stochastic approximation...')
    stoFit <- list()
    tmp$method <- "SA"


    if(is.null(CONTROL$STEP0)) CONTROL$STEP0 <- sqrt(1/dims$d)
    tmp$control <- CONTROL

    sa_args <-  c(list(FREQ        = DATA$freq,
                       VALFREQ     = if(is.null(DATA$valfreq)) DATA$freq else DATA$valfreq,
                       N           = dims$n,
                       C_VEC       = dims$cat,
                       CONSTRMAT   = MODEL$CONSTRMAT,
                       CONSTRLOGSD = MODEL$CONSTRLOGSD,
                       LLC         = MODEL$LLC,
                       THETA_INIT  = tmp$init,
                       NTHR        = dims$nthr,
                       NLOAD       = dims$nload,
                       NCORR       = dims$ncorr,
                       NVAR        = dims$nvar
    ), unclass(CONTROL), VERBOSE = VERBOSE, UFREQ = list(DATA$ufreq),
    VALUFREQ = list(if(is.null(DATA$valfreq)) DATA$ufreq else numeric(0)))

    opt <- do.call(cpp_plSA2, sa_args)
    end_time <- Sys.time()
    tmp$RTime <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

    converge <- list()
    converge$convergence_full <- opt$convergence_full
    converge$convergence_burn <- opt$convergence_burn
    converge$proj_after_burn  <- opt$proj_after_burn
    converge$finite_gr        <- opt$finite_gr
    converge$neg_pdiff        <- opt$neg_pdiff


    if(VERBOSE & !opt$convergence_burn)    message(paste0("Burn-in period did not reach tolerance level ", sa_args$TOL_BURN, ". Stopped at ", sa_args$BURNE, " epochs Try increasing BURNE or STEP0."))
    if(VERBOSE & opt$proj_after_burn)      message("Projections performed after the burn-in. Trajectories might be unstable. Check trajectories visually and evaluate decreasing STEP0.")
    if(VERBOSE & opt$convergence_full==0)  message(paste0("Algorithm ended before reaching tolerance level ", sa_args$TOL_END, ". Stopped at ", sa_args$MAXE, ". epochs Try increasing MAXE or STEP0."))
    if(VERBOSE & opt$convergence_full==-1) message(paste0("Stopped because of divergent trajectories. Decrease STEP0."))
    if(VERBOSE & !opt$finite_gr)           message(paste0("The algorithm automatically dealt with divergent trajectories. Check trajectories visually and evaluate decreasing STEP0."))
    if(VERBOSE & opt$neg_pdiff)            message("Possible divergent trajectories detected monitoring the full negative pairwise likelihood. Check trajectories visually and evaluate decreasing STEP0.")


    stoFit$convergence  <- converge
    stoFit$path$iters   <- opt$path_iters
    stoFit$path$theta   <- opt$path_theta
    stoFit$path$avtheta <- opt$path_avtheta
    stoFit$path$nll     <- opt$path_nll
    stoFit$nll          <- opt$nll
    stoFit$last_iter    <- opt$last_iter
    stoFit$burnt        <- opt$burnt
    tmp$stoFit          <- stoFit
    tmp$theta           <- opt$avtheta


    if(VERBOSE) message('Done! (', round(tmp$RTime,2),' secs)')
    RcppParallel::setThreadOptions(numThreads = 1)
    return(tmp)







  }
}
