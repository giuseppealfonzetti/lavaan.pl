############################################
# Test file for utils.R and generalUtils.h #
# functions as exported by exportedFuns.h  #                                 #
############################################


set.seed(123)

# p = number of items, q = number of latent variables, n = number of observations
p <- 12; q <- 3

# Thresholds vector for each item
thr <- c(-1.5, 0, 1.5)
cat <- rep(length(thr)+1, p)
nthr <- sum(cat)-p

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### free correlation matrix and latent variances ######
{
  stdlv <- FALSE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1.5,6,2), c(2,10, 3)),
              list(c(3,1), c(3,7,2))
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
  local({
    expect_equal(length(theta), d)
  })


  local({
    expect_equal(
      get_lambda(
        THETA = theta,
        NTHR = nthr,
        NLOAD = nload
      ),
      Load[is.na(A)]

    )
  })

  local({
    expect_equal(
      lavaan.pl:::cpp_loadings_theta2mat(
        THETA = theta,
        CONSTRMAT = A,
        LLC = llc,
        NTHR = nthr,
        NLOAD = nload
      ),
      Load

    )
  })

  local({
    expect_equal(
      lavaan.pl:::cpp_latvar_theta2mat(
        THETA = theta,
        CONSTRLOGSD = constr_lsd,
        NTHR = nthr,
        NLOAD=nload,
        NCORR = ncorr,
        NVAR = nvar,
        Q=q
      ),
      S
    )
  })

  local({
    expect_equal(
      get_S(
        THETA = theta,
        CONSTRLOGSD = constr_lsd,
        NTHR = nthr,
        NLOAD=nload,
        NCORR = ncorr,
        NVAR = nvar,
        Q=q
      ),
      S
    )
  })

  local({
    expect_equal(
      get_R(
        THETA = theta,
        NTHR = nthr,
        NLOAD=nload,
        NCORR = ncorr,
        NVAR = nvar,
        Q=q
      ),
      R
    )
  })

  local({
    expect_equal(
      get_corr(
        THETA = theta,
        NTHR = nthr,
        NLOAD=nload,
        NCORR = ncorr,
        NVAR = nvar,
        Q=q
      ),
      R[upper.tri(R)]
    )
  })

}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
{
  stdlv <- TRUE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1,6,2), c(2,10, 3)),
              list(c(3,1), c(3,7,2))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv, llc)

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
  local({
                     expect_equal(length(theta), d)
                   })


  local({
                     expect_equal(
                       get_lambda(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load[is.na(A)]

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_loadings_theta2mat(
                         THETA = theta,
                         CONSTRMAT = A,
                         LLC = llc,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2mat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_S(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_R(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R
                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2dmat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       Dmat
                     )
                   })

  local({
                     expect_equal(
                       get_corr(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R[upper.tri(R)]
                     )
                   })

}

#### (STDLV=TRUE, CORRFLAG=TRUE) ####
#### free correlation matrix but all latent variances fixed ######
#### some variances specified by the user ####
{
  stdlv <- TRUE
  corrflag <- TRUE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1,6,2), c(2,10, 3)),
              list(c(3,1), c(3,7,2))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv, llc)

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
  local({
                     expect_equal(length(theta), d)
                   })


  local({
                     expect_equal(
                       get_lambda(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load[is.na(A)]

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_loadings_theta2mat(
                         THETA = theta,
                         CONSTRMAT = A,
                         LLC = llc,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2mat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_S(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_R(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R
                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2dmat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       Dmat
                     )
                   })

  local({
                     expect_equal(
                       get_corr(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R[upper.tri(R)]
                     )
                   })

}

#### (STDLV=FALSE, CORRFLAG=FALSE) ####
#### no correlation matrix but some free latent variances ######
{
  stdlv <- FALSE
  corrflag <- FALSE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1,6,2), c(2,10, 3)),
              list(c(3,1), c(3,7,2))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv, llc)

  # Draw some random loadings according to constraints
  Load <- gen_loadings(CONSTRMAT = A, STDLV = stdlv, LLC=llc)
  nload <- sum(is.na(A))

  # Generate random latent correlation matrix
  tcorrvec <- rep(0, q*(q-1)/2); if(corrflag) tcorrvec <- rnorm(q*(q-1)/2)
  ncorr <- if(corrflag)q*(q-1)/2 else 0
  R <- lavaan.pl:::cpp_latvar_vec2cmat(VEC=tcorrvec, NCORR=ncorr, Q=q)


  # Generate random latent variances
  constr_var <- rep(NA, q); constr_var[2] <- 1.5
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
  local({
                     expect_equal(length(theta), d)
                   })


  local({
                     expect_equal(
                       get_lambda(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load[is.na(A)]

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_loadings_theta2mat(
                         THETA = theta,
                         CONSTRMAT = A,
                         LLC = llc,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2mat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_S(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_R(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R
                     )
                   })

  local({
                     expect_equal(
                       get_corr(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R[upper.tri(R)]
                     )
                   })

}

#### (STDLV=FALSE, CORRFLAG=TRUE) ####
#### no correlation matrix but all free latent variances ######
{
  stdlv <- FALSE
  corrflag <- FALSE

  # Simple loading matrix constraints
  A <- build_constrMat(P = p, Q = q, STRUCT = 'simple')
  llc <- list(list(c(2,1), c(1,6,2), c(2,10, 3)),
              list(c(3,1), c(3,7,2))
  )
  A <- lavaan.pl:::check_cnstr_loadings(A, stdlv, llc)

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
  local({
                     expect_equal(length(theta), d)
                   })


  local({
                     expect_equal(
                       get_lambda(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load[is.na(A)]

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_loadings_theta2mat(
                         THETA = theta,
                         CONSTRMAT = A,
                         LLC = llc,
                         NTHR = nthr,
                         NLOAD = nload
                       ),
                       Load

                     )
                   })

  local({
                     expect_equal(
                       lavaan.pl:::cpp_latvar_theta2mat(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_S(
                         THETA = theta,
                         CONSTRLOGSD = constr_lsd,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       S
                     )
                   })

  local({
                     expect_equal(
                       get_R(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R
                     )
                   })

  local({
                     expect_equal(
                       get_corr(
                         THETA = theta,
                         NTHR = nthr,
                         NLOAD=nload,
                         NCORR = ncorr,
                         NVAR = nvar,
                         Q=q
                       ),
                       R[upper.tri(R)]
                     )
                   })

}

