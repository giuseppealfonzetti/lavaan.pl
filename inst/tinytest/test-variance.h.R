# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 1000

# Thresholds vector for each item
thr <- c(-1.5, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

full_ngr <- function(par_vec){
  n <- dims$n
  mod <- lavaan.pl:::cpp_multiThread_completePairwise(
    N          = dims$n,
    C_VEC      = dims$cat,
    CONSTRMAT  = constr_list$CONSTRMAT,
    CONSTRLOGSD= constr_list$CONSTRLOGSD,
    FREQ       = f,
    THETA      = par_vec,
    CORRFLAG   = constr_list$CORRFLAG,
    NTHR       = dims$nthr,
    NLOAD      = dims$nload,
    NCORR      = dims$ncorr,
    NVAR       = dims$nvar,
    GRFLAG     = 1,
    SILENTFLAG = 1
  )
  out <- mod$iter_ngr/dims$n
  return(out)
}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv)


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
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -.5, .1)
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
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  dat <- lavaan.pl:::check_data(dat)
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  RcppParallel::setThreadOptions(numThreads = 1)
  HJ <- lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = theta,
    FREQ = f,
    DATA = dat,
    C_VEC = cat,
    CONSTRMAT = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    N = dims$n,
    CORRFLAG = constr_list$CORRFLAG,
    NTHR = dims$nthr,
    NLOAD = dims$nload,
    NCORR= dims$ncorr,
    NVAR = dims$nvar
  )

  # check total gradient
  local({
    expect_equal(
      -HJ$gradH/n,
      full_ngr(theta)
    )

    expect_equal(
      -HJ$gradJ/n,
      full_ngr(theta)
    )
  })

}

hessian_setups <- list(
  "free latent variances"         = list(stdlv = FALSE, corrflag = TRUE),
  "standardised latent variables" = list(stdlv = TRUE,  corrflag = TRUE),
  "orthogonal latent variables"   = list(stdlv = TRUE,  corrflag = FALSE)
)
for (setup_name in names(hessian_setups)) {
  setup <- hessian_setups[[setup_name]]
  set.seed(1)
  p <- 9; q <- 3; n <- 1000
  thr <- c(-1, 0.2, 1.4)
  cat <- rep(length(thr) + 1, p)
  constrvar <- if (setup$stdlv) NULL else rep(NA, q)

  A <- lavaan.pl:::check_cnstr_loadings(build_constrMat(P = p, Q = q, STRUCT = "crossed", CROSS = 1),
                            STDLV = setup$stdlv)
  Load <- gen_loadings(CONSTRMAT = A, STDLV = setup$stdlv, LB = .3, UB = .6)
  ncorr <- if (setup$corrflag) q * (q - 1) / 2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC = rnorm(q * (q - 1) / 2) * setup$corrflag, NCORR = ncorr, Q = q)
  lsd <- lavaan.pl:::check_cnstr_latvar(constrvar, q, setup$stdlv)
  lsd[is.na(lsd)] <- rnorm(sum(is.na(lsd)), -.3, .1)
  S <- diag(exp(lsd), q) %*% R %*% diag(exp(lsd), q)

  dat <- sim_data(SAMPLE_SIZE = n, LOADINGS = Load, THRESHOLDS = thr, LATENT_COV = S)
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT = A, CORRFLAG = setup$corrflag, CONSTRVAR = constrvar,
                                  STDLV = setup$stdlv))
  dims <- lavaan.pl:::check_dims(dat, constr_list)
  f <- compute_frequencies(Y = dat, C_VEC = cat)

  theta <- get_theta(THRESHOLDS = rep(thr, p), LOADINGS = Load, LATENT_COV = S, CAT = cat, CONSTRMAT = A,
                     CONSTRVAR = constrvar, CORRFLAG = setup$corrflag, STDLV = setup$stdlv)
  theta <- theta + rnorm(length(theta), 0, .05)

  RcppParallel::setThreadOptions(numThreads = 2)
  HJ <- lavaan.pl:::cpp_sample_estimators_HJ(
    THETA = theta, FREQ = f, DATA = dat, C_VEC = cat,
    CONSTRMAT = constr_list$CONSTRMAT, CONSTRLOGSD = constr_list$CONSTRLOGSD,
    N = n, CORRFLAG = constr_list$CORRFLAG,
    NTHR = dims$nthr, NLOAD = dims$nload, NCORR = dims$ncorr, NVAR = dims$nvar
  )
  Hnum <- numDeriv::jacobian(full_ngr, theta)

  local({
    expect_equal(-HJ$gradH / n, full_ngr(theta))
    expect_equal(HJ$H, t(HJ$H))
    expect_equal(HJ$H, (Hnum + t(Hnum)) / 2, tolerance = 1e-4)
  })
}
RcppParallel::setThreadOptions(numThreads = 1)
