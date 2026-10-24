#' Fit Confirmatory Factor Analysis Models
#'
#' Fit a Confirmatory Factor Analysis (CFA) model using the `lavaan.pl()` engine
#' using the [`lavaan`] framework. By default, a quasi-Newton BFGS optimiser
#' from the [`ucminf`] package is used. For large data sets, the stochastic
#' approximation algorithm can be used.
#'
#' @details Not all [`lavaan`] options can be used at present. Some options of
#'   interest are:
#' \describe{
#'   \item{\code{information}}{The information matrix to use. Only `"observed"` is currently supported.}
#'   \item{\code{se}}{The standard error method to use. Only `"robust.huber.white"` is currently supported.}
#'   \item{\code{test}}{No GOF tests are available as of now, so this is set to `"none"`.}
#' }
#'
#' @param model A description of the user-specified model. Typically, the model
#'   is described using the lavaan model syntax. See [`lavaan::model.syntax`]
#'   for more information. Alternatively, a parameter table (eg. the output of
#'   the [`lavaan::lavaanify()`] function) is also accepted.
#' @param data A data frame containing the observed variables used in the model.
#'   Variables must be declared as ordered factors.
#' @param std.lv If `TRUE`, the metric of each latent variable is determined by
#'   fixing their (residual) variances to 1.0. If `FALSE`, the metric of each
#'   latent variable is determined by fixing the factor loading of the first
#'   indicator to 1.0.
#' @param start A vector of starting values to use (in the order of free
#'   loadings, thresholds, and then factor correlations). If not provided, the
#'   starting point is computed according to [`fit_plFA()`]'s `INIT_METHOD`.
#' @param estimator The estimator is `PML`. If any other estimator is provided,
#'   then [`lavaan::cfa()`] is used.
#' @param estimator.args A list of arguments for [fit_plFA()]--see the help file
#'   for more details. Possible options are:
#'   \describe{
#'     \item{\code{method}}{One of \code{"ucminf"} (default) for the quasi-Newton BFGS optimiser, or \code{"SA"} for stochastic approximation.}
#'     \item{\code{init_method}}{One of \code{"SA"} (default) for stochastic approximation, \code{"custom"} for user-defined starting values, or \code{"standard"} for standard starting values.}
#'     \item{\code{cpp_control_init}}{A list of control parameters for the initialisation algorithm.}
#'     \item{\code{ncores}}{The number of cores to use for parallel computation.}
#'     \item{\code{valdata}}{Validation data.}
#'     \item{\code{computevar_numderiv}}{If `TRUE`, the asymptotic variance-covariance matrix is computed using numerical derivatives.}
#'   }
#' @param control A list of control parameters for the estimation algorithm. See
#'   [fit_plFA()] for more information.
#' @param verbose If `TRUE`, print additional information during the estimation
#'   process.
#' @param ... Additional arguments to be passed to [lavaan()].
#'
#' @returns A `plFAlavaan` object, which is a subclass of the `lavaan` class.
#'   Therefore, all methods available for `lavaan` objects are expected to be
#'   compatible with `plFAlavaan` objects.
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
#' fit <- cfa("eta =~ y1 + y2 + y3 + y4 + y5", LSAT, std.lv = TRUE)
#' summary(fit)
#'
#' @export
cfa <- function(
    model,
    data,
    std.lv = FALSE,
    estimator = "PML",
    estimator.args = list(
      method = c("ucminf", "SA"),
      init_method = c("SA", "custom", "standard"),
      cpp_control_init = NULL,
      ncores = 1,
      valdata = NULL,
      computevar_numderiv = FALSE
    ),
    start = NULL,
    control = list(),
    verbose = FALSE,
    ...
) {

  # Validate arguments ---------------------------------------------------------
  method <- estimator.args$method
  if (is.null(method)) method <- "ucminf"
  method <- rlang::arg_match(method, c("ucminf", "SA"))
  init_method <- estimator.args$init_method
  if (is.null(init_method)) init_method <- "SA"
  ncores <- estimator.args$ncores
  if (is.null(ncores)) ncores <- 1
  valdata <- estimator.args$valdata
  if (method == "SA") {
    cpp_control_main <- control
    control <- list()
  }
  cpp_control_init <- estimator.args$cpp_control_init
  computevar_numderiv <- estimator.args$computevar_numderiv
  if (is.null(computevar_numderiv)) computevar_numderiv <- FALSE

  # Initialise {lavaan} model object -------------------------------------------
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$std.lv <- std.lv
  lavargs$estimator <- estimator
  lavargs$do.fit <- FALSE

  # Cheeky way to see if vars is passed, then no need to recompute it
  vars <- lavargs$vars
  lavargs$vars <- NULL  # empty it!

  if (estimator == "PML") {
    if ("information" %in% names(lavargs)) {
      if (lavargs$information != "observed") {
        cli::cli_abort("Only 'observed' information is currently supported.")
      }
    }
    lavargs$information <- "observed"
    if ("se" %in% names(lavargs)) {
      if (!(lavargs$se %in% c("robust.huber.white", "none"))) {
        cli::cli_alert_info("Only 'robust.huber.white' is currently supported.")
        lavargs$se <- "robust.huber.white"
      }
    } else {
      lavargs$se <- "robust.huber.white"
    }
    if (!("test" %in% names(lavargs))) {
      lavargs$test <- "none"
    }
    fit0 <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)
  } else{
    # In case user selects "DWLS" or "WLSMV", revert to lavaan::cfa()
    lavargs$do.fit <- TRUE
    fit <- do.call(get("cfa", envir = asNamespace("lavaan")), lavargs)
    return(fit)
  }

  # Check any ordinal data or not?
  if (!all(fit0@Data@ov$type == "ordered")) {
    cli::cli_abort("All measurement items must be declared as ordered factors.")
  }
  if (fit0@Data@ngroups > 1L) {
    cli::cli_abort("Multigroup analysis is not currently supported.")
  }

  # Fit lavaan.pl -------------------------------------------------------------------
  D <- fit0@Data@X[[1]] - 1

  # Build the p x q loading constraints matrix
  FREE <- lavaan::inspect(fit0, what = "free", add.class = FALSE)
  lambda <- lavaan::inspect(fit0, what = "est")$lambda
  p <- nrow(lambda)
  q <- ncol(lambda)
  A <- lambda
  class(A) <- "matrix"
  A[] <- NA
  A[FREE$lambda == 0] <- lambda[FREE$lambda == 0]

  # CORRFLAG
  psi <- lavaan::inspect(fit0, "est")$psi
  LTpsi <- psi[lower.tri(psi)]
  if (length(LTpsi) > 0) {
    if (any(LTpsi == 0 & FREE$psi[lower.tri(psi)] == 0))
      corrflag <- FALSE
    else
      corrflag <- TRUE
  } else {
    corrflag <- FALSE
  }
  constrvar <- diag(psi)
  constrvar[diag(FREE$psi) > 0] <- NA

  # Linear loadings constraints
  llc <- NULL
  pt <- lavaan::partable(fit0)
  ptlc <- pt[pt$op == "==", ]
  Lambdaid <- lavaan::inspect(fit0, "free")$lambda

  if (nrow(ptlc) > 0L) {
    llc <- list()
    coord1 <- ptlc$lhs
    coord2 <- strsplit(ptlc$rhs, "\\+")[[1]]
    for (i in seq_along(coord1)) {
      idx <- pt$id[pt$label == coord1[i]]
      llc[[i]] <- list(which(Lambdaid == idx, arr.ind = TRUE))

      for (j in seq_along(coord2)) {
        cz <- strsplit(coord2[j], "\\*")[[1]]
        if (length(cz) == 1L) {
          cz <- c("1", cz)
        }
        idx <- pt$id[pt$label == cz[2]]
        tmp <- which(Lambdaid == idx, arr.ind = TRUE)
        llc[[i]] <- c(llc[[i]], list(as.numeric(c(cz[1], tmp))))
      }
    }
  }

  # Build constraint list
  constr_list <- list(
    CONSTRMAT = A,
    CONSTRVAR = constrvar,  # FIXME: Correct way to specify variance constraints?
    CORRFLAG = corrflag,
    STDLV = std.lv,
    LLC = llc
  )

  # Starting values are always computed (for lav<->lavaan.pl indices). Note: the
  # values in lavaan.pl are always in the order of thresholds, loadings, latent
  # correlations, then latent variances.
  startlist <- lavaan::inspect(fit0, what = "start", add.class = FALSE)
  lambda <- startlist$lambda[FREE$lambda > 0]
  tau <- startlist$tau[FREE$tau > 0]
  rho <- startlist$psi[FREE$psi > 0 & lower.tri(FREE$psi)]
  logsd <- log(sqrt(startlist$psi[FREE$psi > 0 & diag(TRUE, nrow(FREE$psi))]))
  init <- c(tau, lambda, rho, logsd)

  # init <- get_theta(
  #   THRESHOLDS = as.numeric(startlist$tau),
  #   LOADINGS = startlist$lambda,
  #   LATENT_COV = startlist$psi,
  #   CAT = fit0@Data@ov$nlev,
  #   CONSTRMAT = A,
  #   CONSTRVAR = constrvar,
  #   CORRFLAG = corrflag,
  #   STDLV = std.lv
  # )
  idx_lav2plFA <- c(
    FREE$tau[FREE$tau > 0],
    FREE$lambda[FREE$lambda > 0],
    FREE$psi[FREE$psi > 0 & lower.tri(FREE$psi)],
    FREE$psi[FREE$psi > 0 & diag(TRUE, nrow(FREE$psi))]
  )
  idx_plFA2lav <- order(idx_lav2plFA)

  if (!is.null(start)) {
    init <- as.numeric(start)
    init <- init[idx_lav2plFA]
  }

  # Send to lavaan.pl
  fit1 <- fit_plFA(
    DATA = D,
    CONSTR_LIST = constr_list,
    METHOD = method,
    INIT_METHOD = init_method,
    NCORES = ncores,
    VALDATA = valdata,
    CONTROL = control,
    CPP_CONTROL_MAIN = cpp_control_main,
    CPP_CONTROL_INIT = cpp_control_init,
    INIT = init,
    VERBOSE = verbose
  )

  # Compute standard errors ----------------------------------------------------
  if (is.null(vars)) {
    if (lavargs$se != "none") {
      vars <- computeVar(
        OBJ = fit1,
        DATA = D,
        NUMDERIV = computevar_numderiv,
        OPTION = "transformed",
        VERBOSE = verbose
      )
    }
  }

  out <- create_lav_from_fitplFA(fit0, fit1, vars, D, idx_plFA2lav)
  new("plFAlavaan", out)
}

create_lav_from_fitplFA <- function(fit0, fit1, vars, D, idx_plFA2lav) {

  # Get coefficients and standard errors ---------------------------------------
  FREE <- lavaan::inspect(fit0, what = "free")
  n <- fit0@Data@nobs[[1]]  # FIXME: Group 1 only

  parlist <- extract_par(
    THETA = fit1@theta,
    OPTION = "list",
    CONSTRMAT = fit1@cnstr@loadings,
    CONSTRLOGSD = fit1@cnstr@loglatsd,
    LLC = fit1@cnstr@llc,
    NTHR = fit1@dims@nthr,
    NLOAD = fit1@dims@nload,
    NCORR = fit1@dims@ncorr,
    NVAR = fit1@dims@nvar
  )
  lambda <- parlist$loadings[FREE$lambda > 0]
  tau <- parlist$thresholds[FREE$tau > 0]
  psi <- c(
    parlist$latent_correlations[FREE$psi > 0 & diag(TRUE, nrow(FREE$psi))],
    parlist$latent_correlations[FREE$psi > 0 & lower.tri(FREE$psi)]
  )

  x <- c(lambda, tau, psi)
  # x <- fit1@theta[idx_plFA2lav]

  if (is.null(vars)) {
    SE <- NULL
    vcov <- NULL
  } else {
    SE <- sqrt(vars$asymptotic_variance / n)[idx_plFA2lav]
    vcov <- vars$vcov / n
    vcov <- vcov[idx_plFA2lav, idx_plFA2lav]
  }

  # Change Model and implied slots ---------------------------------------------
  fit0@Model <- lavaan::lav_model_set_parameters(fit0@Model, x)
  fit0@implied <- lavaan::lav_model_implied(fit0@Model)

  # Change ParTable slot -------------------------------------------------------

  # Find Theta matrix (residuals) and Sigmay (implied covariance matrix)
  thetadiag <- diag(fit0@Model@GLIST$theta)
  Sigmay <- fit0@implied$cov[[1]]  # FIXME: Group 1 only

  # Change ParTable and pta slots
  pt <- lavaan::partable(fit0)
  pt$est[pt$free > 0] <- x
  pt$se <- 0
  if (length(SE) > 0) pt$se[pt$free > 0] <- SE
  pt$start[pt$free > 0] <- fit1@init[idx_plFA2lav]

  # Put the diag theta values in the pt
  ov_names <- fit0@Data@ov.names[[1]]  # FIXME: Group 1 only
  pt$est[pt$lhs %in% ov_names &
           pt$rhs %in% ov_names &
           pt$lhs == pt$rhs &
           pt$op == "~~"] <- thetadiag

  # Manually change the slack column
  slack_values <- as.vector(fit0@Model@con.jac %*% x - fit0@Model@ceq.rhs)
  pt$est[pt$op == "=="] <- slack_values
  fit0@ParTable <- as.list(pt)
  fit0@pta$names <- names(pt)

  # Change Options slot --------------------------------------------------------
  fit0@Options$optim.method <- fit1@method
  fit0@Options$do.fit <- TRUE

  # Change Fit slot (depends whether it is numFit or stoFit) -------------------
  fit0@Fit@x <- x
  fit0@Fit@TH[[1]] <- tau  # FIXME: Group 1
  fit0@Fit@est <- pt$est
  fit0@Fit@se <- pt$se
  fit0@Fit@start <- pt$start
  if (fit1@method == "ucminf") {
    fit0@Fit@iterations <- as.integer(fit1@numFit$info["neval"])
    fit0@Fit@converged <- fit1@numFit$convergence == 1L  # FIXME: Check!!
    fit0@Fit@fx <- fit1@numFit$value
  } else if (fit1@method == "SA") {
    fit0@Fit@iterations <- as.integer(fit1@stoFit@last_iter)
    fit0@Fit@converged <- fit1@stoFit@convergence@convergence_full == 1L  # FIXME: Check!!
    fit0@Fit@fx <- fit1@stoFit@nll
  }
  fit0@Fit@Sigma.hat[[1]] <- Sigmay  # implied variance-covariance matrix for group 1!!

  # Change optim slot ----------------------------------------------------------
  fit0@optim$x <- x
  fit0@optim$dx <- 0
  fit0@optim$npar <- length(x)
  fit0@optim$fx <- fx <- fit0@Fit@fx
  fit0@optim$fx.group <- fit0@Fit@fx.group
  fit0@optim$iterations <- fit0@Fit@iterations
  fit0@optim$converged <- fit0@Fit@converged

  # Change loglik slot ---------------------------------------------------------
  if (fit1@method == "ucminf") {
    fit0@loglik$loglik <- fit1@numFit$value
  } else if (fit1@method == "SA") {
    fit0@loglik$loglik <- fit1@stoFit@nll
  }
  fit0@loglik$estimator <- "ML"  # FIXME: to turn off warning for now
  if (!is.null(vars)) {
    fit0@loglik$AIC <- get_AIC(NLL = fit0@loglik$loglik, INVH = vars$invH,
                               J = vars$J)
    fit0@loglik$BIC <- get_BIC(NLL = fit0@loglik$loglik, INVH = vars$invH,
                               J = vars$J, N = n)
  }

  # Change vcov slot -----------------------------------------------------------
  fit0@vcov$vcov <- vcov

  # Change test slot -----------------------------------------------------------
  Options <- fit0@Options
  Options$optim.method <- "nlminb"  # hack to get no warnings from lavaan (ucminf not recognised.)
  # fxval <- fx
  # attr(fx, "fx.pml") <- fxval
  # attr(fx, "fx.group") <- fxval
  # attr(x, "fx") <- fx
  # VCOV <- vcov
  # attr(VCOV, "B0.group")[[1]] <- vars$J[idx_plFA2lav, idx_plFA2lav]
  # attr(VCOV, "E.inv") <- vars$invH[idx_plFA2lav, idx_plFA2lav]

  start_time <- Sys.time()
  fit0@test <- lavaan___lav_model_test(
    lavoptions     = Options,
    lavmodel       = fit0@Model,
    lavsamplestats = fit0@SampleStats,
    lavdata        = fit0@Data,
    lavpartable    = fit0@ParTable,
    lavcache       = fit0@Cache,
    lavimplied     = fit0@implied,
    lavh1          = fit0@h1,
    # x              = x,  # by not including these, lavaan will recompute it
    # VCOV           = VCOV,
    lavloglik      = fit0@loglik
  )
  end_time <- Sys.time()
  timing_test <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

  # Change baseline slot -------------------------------------------------------
  fit0@baseline <- lavaan___lav_step15_baseline(
    lavoptions = fit0@Options,
    lavsamplestats = fit0@SampleStats,
    lavdata = fit0@Data,
    lavcache = fit0@Cache,
    lavh1 = fit0@h1,
    lavpartable = fit0@ParTable
  )

  # Change timing slot ---------------------------------------------------------
  fit0@timing$optim <- fit0@timing$optim + fit1@RTime
  fit0@timing$vcov <- vars$RTime
  fit0@timing$test <- timing_test
  fit0@timing$total <- sum(unlist(fit0@timing))

  # Include the entire output of lavaan.pl fit --------------------------------------
  fit0@external <- list(
    lavaan.pl = fit1,
    vars = vars,
    D = D,
    idx_plFA2lav = idx_plFA2lav
  )

  # Output ---------------------------------------------------------------------
  fit0
}

# Export lavaan functions
lavaan___lav_step15_baseline <-
  utils::getFromNamespace("lav_step15_baseline", "lavaan")
lavaan___lav_model_test <-
  utils::getFromNamespace("lav_model_test", "lavaan")
