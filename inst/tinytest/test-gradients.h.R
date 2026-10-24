###################################################
# Test file for gradients.h based on wrappers for #
# pair contributions in pairs.h as exported by    #
# exportedFuns.h                                  #
###################################################

# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3; n <- 100

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p
pair_nll <- function(PAR, K, L){
  pair <- lavaan.pl:::cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    LLC         = constr_list$LLC,
    C_VEC       = cat,
    THETA       = PAR,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    GRADFLAG    = 0
  )
  out <- pair$nll/n
  return(out)
}
pair_ngr <- function(PAR, K, L){
  pair <- lavaan.pl:::cpp_compute_pair_ext(
    CONSTRMAT   = constr_list$CONSTRMAT,
    CONSTRLOGSD = constr_list$CONSTRLOGSD,
    LLC         = constr_list$LLC,
    C_VEC       = cat,
    THETA       = PAR,
    NTHR        = nthr,
    NLOAD       = nload,
    NCORR       = ncorr,
    NVAR        = nvar,
    K           = K,
    L           = L,
    PAIRS_TABLE = f,
    GRADFLAG    = 1
  )
  out <- pair$ngr/n
  return(out)
}

full_ngr <- function(par_vec){
  n <- dims$n
  mod <- lavaan.pl:::cpp_multiThread_completePairwise(
    N          = dims$n,
    C_VEC      = dims$cat,
    CONSTRMAT  = constr_list$CONSTRMAT,
    CONSTRLOGSD= constr_list$CONSTRLOGSD,
    LLC        = constr_list$LLC,
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
  llc <- NULL
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)


  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })

}
#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
#### with linear constraints on loadings ####
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(.5,6,2)),
              list(c(7,2), c(1,11,3))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)


  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })
}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
#### with linear constraints on loadings ####
if(1){
  set.seed(123)
  stdlv <- TRUE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(.5,6,2), c(.2,11,3)),
              list(c(7,2), c(1,11,3))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })
}

### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
#### some variances specified by the user ####
#### with linear constraints on loadings ####
if(1){
  set.seed(123)
  stdlv <- TRUE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(.5,6,2), c(.2,11,3)),
              list(c(7,2), c(1,11,3))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q); constr_var[1:2] <- .5
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })
}

#### (STDLV=FALSE, CORRFLAG=FALSE) ####
#### no correlation matrix but some free latent variances ######
#### with linear constraints on loadings ####
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- FALSE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(.5,6,2), c(.2,11,3)),
              list(c(7,2), c(1,11,3))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q); constr_var[2] <- .5
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -1, .1)
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })
}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### no correlation matrix but all free latent variances ######
#### with linear constraints on loadings ####
if(1){
  set.seed(123)
  stdlv <- FALSE
  corrflag <- FALSE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(.5,6,2), c(.2,11,3)),
              list(c(7,2), c(1,11,3))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, STDLV=stdlv, LLC = llc)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q);
  constr_lsd <- lavaan.pl:::check_cnstr_latvar(constr_var, q, stdlv)
  nvar  <- sum(is.na(constr_lsd))
  tsdvec <- constr_lsd
  tsdvec[is.na(constr_lsd)] <- rnorm(sum(is.na(constr_lsd)), -1, .1)
  Dmat <- diag(exp(tsdvec),q,q)

  S <- Dmat %*% R %*% Dmat

  # Dimensions
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv))
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
  constr_list <- lavaan.pl:::check_cnstr(list(CONSTRMAT=A, CORRFLAG=corrflag, CONSTRVAR=exp(constr_lsd)^2, STDLV=stdlv, LLC=llc))
  dims <- lavaan.pl:::check_dims(dat, constr_list)

  tot_gr <- rep(0, dims$d)
  for (k in 2:p) {
    for (l in 1:(k-1)) {
      local({
        expect_equal(pair_ngr(PAR=theta, K=k-1, L=l-1), numDeriv::grad(func=pair_nll, x=theta, K=k-1,L=l-1), tolerance = 1e-4)
      })
      tot_gr <- tot_gr+pair_ngr(PAR=theta, K=k-1, L=l-1)
    }
  }

  # check total gradient
  local({
    expect_equal(
      tot_gr,
      full_ngr(theta)
    )
  })
}
