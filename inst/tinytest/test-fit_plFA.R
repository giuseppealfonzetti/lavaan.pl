Sys.setenv("R_TESTS" = "")

# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

###########################
## Check METHOD="SA" init #
###########################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q)
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions

  d <- nthr + nload + ncorr + nvar
  d


  #
  theta <- get_theta(
    THRESHOLDS = rep(thr, p),
    LOADINGS = Load,
    LATENT_COV = S,
    CAT = cat,
    CONSTRMAT = A,
    CONSTRVAR = exp(constr_lsd)^2,
    CORRFLAG = corrflag,
    STDLV = stdlv
  )


  dat <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

  plFit <- suppressWarnings(fit_plFA(
    DATA = set_plfa_data(dat),
    MODEL = do.call(set_plfa_model, constr_list),
    INIT = set_init(METHOD = "standard"),
    CONTROL = set_sa_control(BURNE=0, MAXE=1, PAIRS_PER_ITERATION=1, EPOCHS_PER_CHECK=1),
    VERBOSE = FALSE))

  plFit$stoFit$path$nll


  fit <- lavaan.pl:::cpp_plSA2(
    FREQ = plFit$data$freq,
    VALFREQ = plFit$data$freq,
    N = plFit$dims$n,
    C_VEC = plFit$dims$cat,
    CONSTRMAT=plFit$model$CONSTRMAT,
    CONSTRLOGSD=plFit$model$CONSTRLOGSD,
    THETA_INIT = plFit$init,
    NTHR = plFit$dims$nthr,
    NLOAD= plFit$dims$nload,
    NCORR= plFit$dims$ncorr,
    NVAR= plFit$dims$nvar,
    PAIRS_PER_ITERATION=plFit$control$PAIRS_PER_ITERATION,
    STEP0=plFit$control$STEP0,
    STEP1=plFit$control$STEP1,
    STEP2=plFit$control$STEP2,
    STEP3=plFit$control$STEP3,
    BURNE=plFit$control$BURNE,
    MAXE=plFit$control$MAXE,
    EPOCHS_PER_CHECK = plFit$control$EPOCHS_PER_CHECK,
    TOL_END=plFit$control$TOL_END,
    TOL_BURN = plFit$control$TOL_BURN,
    CHECK_TOL=plFit$control$CHECK_TOL,
    MAX_TOL_COUNTER = plFit$control$MAX_TOL_COUNTER,
    SEED=plFit$control$SEED,
    VERBOSE=FALSE,
    VERBOSE_ITER=FALSE
  )

  # check initial nll from SA
  local({expect_equal(plFit$stoFit$path$nll[1],
                                                       fit$path_nll[1])})


}


###################################
## Check METHOD="SA" vs cpp_plSA ##
###################################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q)
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions

  d <- nthr + nload + ncorr + nvar
  d


  #
  theta <- get_theta(
    THRESHOLDS = rep(thr, p),
    LOADINGS = Load,
    LATENT_COV = S,
    CAT = cat,
    CONSTRMAT = A,
    CONSTRVAR = exp(constr_lsd)^2,
    CORRFLAG = corrflag,
    STDLV = stdlv
  )


  dat <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

  plFit <- suppressWarnings(fit_plFA(
    DATA = set_plfa_data(dat),
    MODEL = do.call(set_plfa_model, constr_list),
    INIT = set_init(METHOD = "standard"),
    CONTROL = set_sa_control(BURNE=1, MAXE=3, PAIRS_PER_ITERATION=8, EPOCHS_PER_CHECK=1),
    VERBOSE = FALSE))

  plFit$stoFit$path$nll


  fit <- lavaan.pl:::cpp_plSA2(
    FREQ = plFit$data$freq,
    VALFREQ = plFit$data$freq,
    N = plFit$dims$n,
    C_VEC = plFit$dims$cat,
    CONSTRMAT=plFit$model$CONSTRMAT,
    CONSTRLOGSD=plFit$model$CONSTRLOGSD,
    THETA_INIT = plFit$init,
    NTHR = plFit$dims$nthr,
    NLOAD= plFit$dims$nload,
    NCORR= plFit$dims$ncorr,
    NVAR= plFit$dims$nvar,
    PAIRS_PER_ITERATION=plFit$control$PAIRS_PER_ITERATION,
    STEP0=plFit$control$STEP0,
    STEP1=plFit$control$STEP1,
    STEP2=plFit$control$STEP2,
    STEP3=plFit$control$STEP3,
    BURNE=plFit$control$BURNE,
    MAXE=plFit$control$MAXE,
    EPOCHS_PER_CHECK = plFit$control$EPOCHS_PER_CHECK,
    TOL_END=plFit$control$TOL_END,
    TOL_BURN = plFit$control$TOL_BURN,
    CHECK_TOL=plFit$control$CHECK_TOL,
    MAX_TOL_COUNTER = plFit$control$MAX_TOL_COUNTER,
    SEED=plFit$control$SEED,
    VERBOSE=FALSE,
    VERBOSE_ITER=FALSE
  )

  # check nll against direct cpp call
  local({
    expect_true(all(plFit$stoFit$path$nll== fit$path_nll))})
}

#######################
## Check METHOD="SA" ##
#######################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q)
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions

  d <- nthr + nload + ncorr + nvar
  d


  #
  theta <- get_theta(
    THRESHOLDS = rep(thr, p),
    LOADINGS = Load,
    LATENT_COV = S,
    CAT = cat,
    CONSTRMAT = A,
    CONSTRVAR = exp(constr_lsd)^2,
    CORRFLAG = corrflag,
    STDLV = stdlv
  )


  dat <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

  plFit <- suppressWarnings(fit_plFA(
    DATA = set_plfa_data(dat),
    MODEL = do.call(set_plfa_model, constr_list),
    INIT = set_init(METHOD = "standard"),
    CONTROL = set_sa_control(),
    VERBOSE = FALSE))

  # finite nll
  local({expect_true(all(is.finite(plFit$stoFit$path$nll)))})
  # finite par
  local({lapply(plFit$stoFit$path$avtheta, function(x){expect_true(all(is.finite(x)))})})
}

############################
## Check INIT_METHOD="SA" ##
############################
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(321)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv)
  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q)
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.25, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions

  d <- nthr + nload + ncorr + nvar
  d


  #
  theta <- get_theta(
    THRESHOLDS = rep(thr, p),
    LOADINGS = Load,
    LATENT_COV = S,
    CAT = cat,
    CONSTRMAT = A,
    CONSTRVAR = exp(constr_lsd)^2,
    CORRFLAG = corrflag,
    STDLV = stdlv
  )


  dat <- sim_data(
    SAMPLE_SIZE = n,
    LOADINGS = Load,
    THRESHOLDS = thr,
    LATENT_COV = S)
  constr_list <- list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv)

  plFit <- suppressWarnings(fit_plFA(
    DATA = set_plfa_data(dat),
    MODEL = do.call(set_plfa_model, constr_list),
    INIT = set_init(METHOD = "SA"),
    CONTROL = set_sa_control(),
    VERBOSE = FALSE))

  # finite nll
  local({expect_true(all(is.finite(plFit$stoFit$path$nll)))})
  # finite par
  local({lapply(plFit$stoFit$path$avtheta, function(x){expect_true(all(is.finite(x)))})})
}
