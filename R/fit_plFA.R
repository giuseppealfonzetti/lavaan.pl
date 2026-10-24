#'Fit factor models for ordinal data with pairwise likelihood methods
#'
#'@description
#'
#'Fit latent variable models for ordinal variable using pairwise likelihood
#'methods running via quasi-Newton BFGS or stochastic approximations.
#'
#'@param DATA Integer data matrix of dimension \eqn{n*p}. Categories must be
#'  coded starting from zero.
#'@param CONSTR_LIST List of constraints. It must contain \tabular{ll}{
#'  \code{CONSTRMAT} \tab \eqn{p*q}-dimensional matrix. Elements set to `NA`
#'  refers to free loading parameters. Elements set to numerical values denote
#'  fixed values constraints. \cr \tab \cr \code{CONSTRVAR} \tab
#'  \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent
#'  variance parameters. Elements set to numerical values denote fixed values
#'  constraints. \cr \tab \cr \code{CORRFLAG} \tab Logical indicator. Set it to
#'  `FALSE` if the latent variables are assumed to be independent. Set it `TRUE`
#'  otherwise. \cr \tab \cr \code{STDLV} \tab Logical indicator. Set it to
#'  `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.
#'  \cr    \tab \cr \code{LLC} \tab Linear loadings constraints. Expects a list
#'  of constraints. See DETAILS. \cr }
#'@param METHOD Label for the method chosen. Possible values are: \tabular{ll}{
#'  \code{'ucminf'} \tab for estimation via the quasi-Newton BFGS optimiser from
#'  the \code{ucminf} package. Used as default method. \cr \tab \cr \code{'SA'}
#'  \tab for estimation via Stochastic Approximations. \cr }
#'@param INIT_METHOD \tabular{ll}{ \code{'custom'} \tab Uses the vector provided
#'  via `INIT` as starting point. \cr \tab \cr \code{'standard'} \tab Uses cold
#'  initialisation where loadings are initialised at `0.5` and the latent
#'  covariance matrix as an Identity matrix. \cr \tab \cr \code{'SA'} \tab
#'  Computes the starting point using a short chain of stochastic updates.
#'  Updates start from the `standard` initialisation point. Set as default
#'  method. \cr }
#'@param NCORES Integer value setting the number of threads to carry out the
#'  estimation.
#'@param CONTROL List of control options to pass to `ucminf`. See
#'  [ucminf::ucminf] documentation.
#'@param CPP_CONTROL_MAIN List of control options to pass to the `SA` optimiser
#'  when `METHOD="SA"`: \tabular{ll}{ \code{PAIRS_PER_ITERATION} \tab  Number of
#'  pairs to draw at each iteration.   \cr \tab \cr \code{MAXT} \tab  Maximum
#'  number of iterations.\cr \tab \cr \code{BURN} \tab  Maximum number of
#'  iterations to burn before trajectory averaging.\cr \tab \cr \code{STEP0}
#'  \tab  Initial step length.\cr }
#'@param CPP_CONTROL_INIT List of control options to pass to the `SA` optimiser
#'  when `INIT_METHOD="SA"`. See `CPP_CONTROL_MAIN` for details.
#'@param VALDATA Validation data used to monitor convergence of stochastic
#'  approximations. If `NULL`, data passed via `DATA` is used for monitoring
#'  purposes.
#'@param INIT Initialising vector. If not provided, the starting point is
#'  computed according to `INIT_METHOD`
#'@param VERBOSE `TRUE` for verbose output
#'
#'@details The argument CONSTR_LIST$LLC is expected to be a list of constraints,
#'  e.g. `CONSTR_LIST$LLC <- list(constraint_1, constraint2, ...)`, where each
#'  constraint is defined by a list itself. For example, to impose the
#'  constraint "L_(2,1)=0.5L_(5,2)+0.25L_(9,3)" you have to write `constraint_1
#'  <- list(c(2,1), c(0.5,5,2), c(0.25, 9,3))`. That is: the first vector of the
#'  list is 2-dimensional and stores the coordinates of the constrained loading.
#'  Each of the successive triplets represent a linear coefficient followed by
#'  the coordinates of the loading of reference. You can set up an arbitrary
#'  number of triplets in each constraint.
#'
#'@return Object of class `plFaFit`.
#'@export
fit_plFA <- function(
    DATA,
    CONSTR_LIST,
    METHOD = c("ucminf", "SA"),
    INIT_METHOD = c("standard", "custom",  "SA"),
    NCORES = 1,
    VALDATA = NULL,
    CONTROL = list(),
    CPP_CONTROL_MAIN = NULL,
    CPP_CONTROL_INIT = NULL,
    INIT = NULL,
    VERBOSE = FALSE
){

  start_time <- Sys.time()

  dat <- check_data(DATA)
  constr_list <- check_cnstr(CONSTR_LIST)
  METHOD <- match.arg(METHOD)
  INIT_METHOD <- match.arg(INIT_METHOD)

  # Identify model dimensions
  dims <- check_dims(dat, constr_list)
  if(METHOD=="SA" && dims$p<20) message("METHOD='SA' not recommended with less than 20 items.")
  INIT_METHOD <- check_init_method(dims, INIT_METHOD, INIT)

  tmp <- structure(
    list(cnstr = list(loadings = constr_list$CONSTRMAT,
                      corrflag = constr_list$CORRFLAG,
                      stdlv    = constr_list$STDLV,
                      loglatsd = constr_list$CONSTRLOGSD,
                      llc      = constr_list$LLC),
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


  if(VERBOSE) message('- Computing frequencies...')
  freq_start_time <- Sys.time()
  tmp$freq <- pairs_freq(dat, dims$cat)
  freq_end_time <- Sys.time()
  tmp$freqTime <- as.numeric(difftime(freq_end_time, freq_start_time, units = 'secs')[1])

  # Check Initialisation
  tmp$init <- check_init(INIT, tmp$freq, dims, constr_list, INIT_METHOD, VERBOSE, CPP_ARGS = CPP_CONTROL_INIT)

  diag_inv_H <- rep(1,length(tmp$init))

  # Numerical optimisation
  if(METHOD == 'ucminf'){

    if(VERBOSE) message('- Optimising with ucminf...')
    tmp$method <- "ucminf"

    # Compute frequency table bivariate patterns

    Rwr_ncl <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = constr_list$CONSTRMAT,
        CONSTRLOGSD= constr_list$CONSTRLOGSD,
        LLC        = constr_list$LLC,
        FREQ       = tmp$freq,
        THETA      = par_vec,
        CORRFLAG   = constr_list$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        SILENTFLAG = 1
      )
      out <- mod$iter_nll/dims$n
      return(out)
    }

    # function for gradient
    Rwr_ngr <- function(par_vec){
      mod <- cpp_multiThread_completePairwise(
        N          = dims$n,
        C_VEC      = dims$cat,
        CONSTRMAT  = constr_list$CONSTRMAT,
        CONSTRLOGSD= constr_list$CONSTRLOGSD,
        LLC        = constr_list$LLC,
        FREQ       = tmp$freq,
        THETA      = par_vec,
        CORRFLAG   = constr_list$CORRFLAG,
        NTHR       = dims$nthr,
        NLOAD      = dims$nload,
        NCORR      = dims$ncorr,
        NVAR       = dims$nvar,
        GRFLAG     = 1,
        SILENTFLAG = 1
      )

      out <- mod$iter_ngradient/dims$n
      return(out)
    }

    # list of ucminf args
    args <- list(
      'par' = tmp$init,
      'fn' = Rwr_ncl,
      'gr' = Rwr_ngr,
      'control' = if (is.null(CONTROL$ctrl)) list() else CONTROL$ctrl,
      'hessian' = if (is.null(CONTROL$hessian)) 0 else CONTROL$hessian)

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
    if(!is.null(VALDATA)){
      if(VERBOSE) message('- Computing frequencies on validation data...')
      tmp$valfreq <- pairs_freq(VALDATA, dims$cat)
    }

    if(VERBOSE) message('- Optimising via pairwise stochastic approximation...')
    stoFit <- list()
    tmp$method <- "SA"


    sa_args <- check_plSA_args(FREQ        = tmp$freq,
                               VALFREQ     = if(is.null(VALDATA)) NULL else tmp$valfreq,
                               DIMS        = dims,
                               CONSTR_LIST = constr_list,
                               LIST        = CPP_CONTROL_MAIN,
                               SETTING     = "main")

    stoFit$control <- sa_args

    sa_args <-  c(list(N           = dims$n,
                       C_VEC       = dims$cat,
                       CONSTRMAT   = constr_list$CONSTRMAT,
                       CONSTRLOGSD = constr_list$CONSTRLOGSD,
                       LLC         = constr_list$LLC,
                       THETA_INIT  = tmp$init,
                       DIH         = diag_inv_H,
                       NTHR        = dims$nthr,
                       NLOAD       = dims$nload,
                       NCORR       = dims$ncorr,
                       NVAR        = dims$nvar
    ), sa_args)
    sa_args$VERBOSE <- VERBOSE

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
