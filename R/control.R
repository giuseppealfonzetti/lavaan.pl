check_data <- function(DATA){
  stopifnot(is.matrix(DATA))
  stopifnot(is.finite(DATA) | is.na(DATA))
  stopifnot(is.numeric(DATA))
  stopifnot(min(DATA, na.rm = TRUE)==0)
  return(round(DATA,0))
}
check_cnstr_corrflag <- function(CORRFLAG){
  if(is.null(CORRFLAG)) CORRFLAG <- TRUE
  stopifnot(CORRFLAG%in%c(0,1))
  return(as.logical(CORRFLAG))
}
check_cnstr_stdlv <- function(STDLV){
  if(is.null(STDLV)) STDLV <- TRUE
  stopifnot(is.logical(STDLV))
  return(STDLV)
}
check_cnstr_loadings <- function(MAT, STDLV, LLC=NULL){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(!is.infinite(sum(MAT, na.rm = TRUE)))

  if(!STDLV){
    for(j in 1:ncol(MAT)){
      if(sum(is.na(MAT[,j]))>0){
        wna <- which(is.na(MAT[,j]))
        fna <- wna[1]
        if(is.null(fna)) break
        if(fna==1){MAT[fna,j] <- 1}else{
          if(all(MAT[1:(fna-1),j]==0))MAT[fna,j] <- 1
        }
      }
    }
  }

  if(!is.null(LLC)){
    for (idx_lc in 1:length(LLC)) {
      MAT[LLC[[idx_lc]][[1]][1], LLC[[idx_lc]][[1]][2]] <- 1
    }
  }

  return(MAT)
}
check_cnstr_latvar <- function(VEC, Q, STDLV){
  if(is.null(VEC) & STDLV) VEC <- rep(1, Q)
  if(is.null(VEC) & !STDLV) VEC <- rep(NA, Q)

  stopifnot(is.vector(VEC))
  stopifnot(!is.infinite(sum(VEC, na.rm = TRUE)))
  stopifnot(VEC[!is.na(VEC)]>0)

  if(STDLV){
    VEC[is.na(VEC)] <- 1
  }

  return(log(sqrt(VEC)))
}
check_cnstr <- function(LIST){
  ls <- list()
  ls$CORRFLAG <- check_cnstr_corrflag(LIST$CORRFLAG)
  ls$STDLV    <- check_cnstr_stdlv(LIST$STDLV)
  ls$LLC <- LIST$LLC
  ls$CONSTRMAT <- check_cnstr_loadings(LIST$CONSTRMAT, STDLV = ls$STDLV, LLC = ls$LLC)
  q <- ncol(ls$CONSTRMAT)
  ls$CONSTRLOGSD <- check_cnstr_latvar(LIST$CONSTRVAR, Q=q, STDLV = ls$STDLV)

  stopifnot(q==length(ls$CONSTRLOGSD))
  return(ls)
}
check_p <- function(P){
  stopifnot(is.numeric(P))
  stopifnot(is.finite(P))
  stopifnot(P>1)
  return(as.integer(P))
}
check_q <- function(Q){
  stopifnot(is.numeric(Q))
  stopifnot(is.finite(Q))
  stopifnot(Q>0)
  return(as.integer(Q))
}
check_n <- function(N){
  stopifnot(is.numeric(N))
  stopifnot(is.finite(N))
  stopifnot(N>0)
  return(as.integer(N))
}
check_loadings <- function(MAT){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(is.finite(MAT))
  stopifnot(nrow(MAT)>ncol(MAT))
  return(MAT)
}
check_thresholds_j <- function(TH){
  stopifnot(is.vector(TH))
  stopifnot(is.numeric(TH))
  stopifnot(!is.unsorted(TH))
  return(TH)
}
check_thresholds <- function(THF, CAT){
  s <- 0
  for(j in 1:length(CAT)){
    check_thresholds_j(THF[(s + 1):(s + CAT[j] - 1)])
    s <- s + CAT[j] - 1
  }
  return(THF)
}
check_latcov <- function(MAT){
  stopifnot(is.matrix(MAT))
  stopifnot(is.numeric(MAT))
  stopifnot(is.finite(MAT))
  stopifnot(nrow(MAT)==ncol(MAT))
  stopifnot(isSymmetric.matrix(MAT))
  #check also positive definiteness
  return(MAT)
}
check_dims <- function(DATA, CONSTR_LIST){
  p <- ncol(DATA)
  n <- nrow(DATA)
  q <- ncol(CONSTR_LIST$CONSTRMAT)
  categories <- apply(DATA, 2, max, na.rm = T) + 1

  nthr <- sum(categories)-p
  nload <- sum(is.na(CONSTR_LIST$CONSTRMAT))
  ncorr <- 0; if(CONSTR_LIST$CORRFLAG==1) ncorr <- q*(q-1)/2
  nvar  <- sum(is.na(CONSTR_LIST$CONSTRLOGSD))
  d <- nthr + nload + ncorr + nvar
  stopifnot(q<p)
  stopifnot(categories>1)




  return(list(
    n=n,
    p=p,
    q=q,
    nthr=nthr,
    nload=nload,
    ncorr=ncorr,
    nvar=nvar,
    d=d,
    cat=categories,
    pairs=p*(p-1)/2
  ))
}
check_theta <- function(THETA){

  stopifnot(is.vector(THETA))
  stopifnot(is.numeric(THETA))
  stopifnot(is.finite(THETA))

  return(THETA)
}
check_init_par <- function(PAR, DIMS, CONSTR_LIST){
  stopifnot(is.vector(PAR))
  stopifnot(is.numeric(PAR))
  stopifnot(is.finite(PAR))
  stopifnot(length(PAR)==DIMS$d)

  L <- cpp_loadings_theta2mat(THETA     = PAR,
                              CONSTRMAT = CONSTR_LIST$CONSTRMAT,
                              LLC       = CONSTR_LIST$LLC,
                              NTHR      = DIMS$nthr,
                              NLOAD     = DIMS$nload)
  check_loadings(L)
  S <- cpp_latvar_theta2mat(THETA       = PAR,
                            CONSTRLOGSD = CONSTR_LIST$CONSTRLOGSD,
                            NTHR        = DIMS$nthr,
                            NLOAD       = DIMS$nload,
                            NCORR       = DIMS$ncorr,
                            NVAR        = DIMS$nvar,
                            Q           = DIMS$q)
  check_latcov(S)
  Sy <- L%*%S%*%t(L)
  if(!all(diag(Sy)<=1)) warning("Negative error variances implied by init parameters")
  PAR <- cpp_sa_proj( THETA=PAR,
                      CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                      CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
                      LLC=CONSTR_LIST$LLC,
                      C_VEC=DIMS$cat,
                      CORRFLAG=CONSTR_LIST$CORRFLAG,
                      NTHR=DIMS$nthr,
                      NLOAD=DIMS$nload,
                      NCORR=DIMS$ncorr,
                      NVAR=DIMS$nvar)

  return(PAR)
}
check_init <- function(INIT, FREQ, DIMS, CONSTR_LIST, INIT_METHOD=c("SA", "custom", "standard"), VERBOSE=FALSE, CPP_ARGS=NULL){

  INIT_METHOD <- check_init_method(DIMS=DIMS, INIT_METHOD = INIT_METHOD, INIT = INIT)

  # Setup fixed values initialisation
  if(!is.null(INIT)) {
    if(VERBOSE) message('- Initialised at user-supplied values.')
  }else{
    lambda0_init <- init_thresholds(DIMS, CONSTR_LIST)
    lambda_init  <- init_loadings(DIMS, CONSTR_LIST)
    transformed_rhos_init  <- init_transformed_latcorr(DIMS, CONSTR_LIST)
    transformed_latsd_init <- init_transformed_latsd(DIMS, CONSTR_LIST)
    INIT <- c(lambda0_init, lambda_init, transformed_rhos_init, transformed_latsd_init)
    INIT <- cpp_sa_proj(THETA=INIT,
                        CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                        CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
                        LLC=CONSTR_LIST$LLC,
                        C_VEC=DIMS$cat,
                        CORRFLAG=CONSTR_LIST$CORRFLAG,
                        NTHR=DIMS$nthr,
                        NLOAD=DIMS$nload,
                        NCORR=DIMS$ncorr,
                        NVAR=DIMS$nvar)
    if(VERBOSE) message(' Done.')
  }

  if(INIT_METHOD=="SA"){
    if(VERBOSE) message('- Initialising with SA ...', appendLF = FALSE)
    CPP_ARGS <- check_plSA_args(FREQ=FREQ, DIMS=DIMS, CONSTR_LIST=CONSTR_LIST, LIST=CPP_ARGS, SETTING = "init")
    INIT <- init_par_with_plSA(DIMS, CONSTR_LIST, CPP_ARGS, INIT)
    if(VERBOSE) message(' Done.')
  }

  INIT <- check_init_par(INIT, DIMS, CONSTR_LIST)
  return(INIT)
}
check_init_method <- function(DIMS, INIT_METHOD = c("SA", "custom", "standard"), INIT =NULL){
  INIT_METHOD <- match.arg(INIT_METHOD)
  if(INIT_METHOD == "custom") INIT_METHOD <- "standard"

  if(INIT_METHOD=="SA"){
    if(DIMS$p<15)INIT_METHOD <- "standard"
  }
  return(INIT_METHOD)
}
check_plSA_args <- function(FREQ, VALFREQ=NULL, DIMS, CONSTR_LIST, LIST=NULL, SETTING=c("init", "main")){
  SETTING <- match.arg(SETTING)
  out <- list()
  out$FREQ <-  FREQ

  if(is.null(VALFREQ)) {out$VALFREQ <- FREQ} else {out$VALFREQ <- VALFREQ}

  if(is.null(LIST$PAIRS_PER_ITERATION) & SETTING=="init") LIST$PAIRS_PER_ITERATION <- 16
  if(is.null(LIST$PAIRS_PER_ITERATION) & SETTING=="main") LIST$PAIRS_PER_ITERATION <- 16

  stopifnot(is.numeric(LIST$PAIRS_PER_ITERATION))
  stopifnot(LIST$PAIRS_PER_ITERATION>0)
  out$PAIRS_PER_ITERATION <- as.integer(LIST$PAIRS_PER_ITERATION)

  if(is.null(LIST$STEP0) & SETTING=="init") LIST$STEP0 <- sqrt(1/DIMS$d)
  if(is.null(LIST$STEP0) & SETTING=="main") LIST$STEP0 <- 1

  stopifnot(is.numeric(LIST$STEP0))
  stopifnot(LIST$STEP0>0)
  out$STEP0 <- LIST$STEP0

  if(is.null(LIST$STEP1)) LIST$STEP1 <- 1
  stopifnot(is.numeric(LIST$STEP1))
  stopifnot(LIST$STEP1>0)
  out$STEP1 <- LIST$STEP1

  if(is.null(LIST$STEP2)) LIST$STEP2 <- 1e-8
  stopifnot(is.numeric(LIST$STEP2))
  stopifnot(LIST$STEP2>0)
  out$STEP2 <- LIST$STEP2

  if(is.null(LIST$STEP3)) LIST$STEP3 <- .75
  stopifnot(is.numeric(LIST$STEP3))
  stopifnot(LIST$STEP3>0)
  out$STEP3 <- LIST$STEP3

  if(is.null(LIST$BURNE) & (SETTING=="init")) LIST$BURNE <- 3
  if(is.null(LIST$BURNE) & (SETTING=="main")) LIST$BURNE <- 20
  stopifnot(is.numeric(LIST$BURNE))
  stopifnot(LIST$BURNE>=0)
  out$BURNE <- as.integer(LIST$BURNE)

  if(is.null(LIST$MAXE) & SETTING=="init") LIST$MAXE <- 5
  if(is.null(LIST$MAXE) & SETTING=="main") LIST$MAXE <- 30
  stopifnot(is.numeric(LIST$MAXE))
  stopifnot(LIST$MAXE>=0)
  out$MAXE <- as.integer(LIST$MAXE)

  if(is.null(LIST$MAX_TOL_COUNTER)) LIST$MAX_TOL_COUNTER <- 1
  stopifnot(is.numeric(LIST$MAX_TOL_COUNTER))
  stopifnot(LIST$MAX_TOL_COUNTER>0)
  out$MAX_TOL_COUNTER <- as.integer(LIST$MAX_TOL_COUNTER)

  if(is.null(LIST$TOL_END)) LIST$TOL_END <- 1e-7
  stopifnot(is.numeric(LIST$TOL_END))
  stopifnot(LIST$TOL_END>0)
  out$TOL_END <- LIST$TOL_END

  if(is.null(LIST$TOL_BURN)) LIST$TOL_BURN <- 5e-7
  stopifnot(is.numeric(LIST$TOL_BURN))
  stopifnot(LIST$TOL_BURN>0)
  out$TOL_BURN <- LIST$TOL_BURN

  if(is.null(LIST$CHECK_TOL) & SETTING=="init") LIST$CHECK_TOL <- FALSE
  if(is.null(LIST$CHECK_TOL) & SETTING=="main") LIST$CHECK_TOL <- TRUE

  stopifnot(is.logical(LIST$CHECK_TOL))
  out$CHECK_TOL <- LIST$CHECK_TOL


  if(is.null(LIST$EPOCHS_PER_CHECK) & SETTING=="init") LIST$EPOCHS_PER_CHECK <- 1
  if(is.null(LIST$EPOCHS_PER_CHECK) & SETTING=="main") LIST$EPOCHS_PER_CHECK <- 5
  stopifnot(is.numeric(LIST$EPOCHS_PER_CHECK))
  stopifnot(LIST$EPOCHS_PER_CHECK>0)
  out$EPOCHS_PER_CHECK <- as.integer(LIST$EPOCHS_PER_CHECK)

  if(is.null(LIST$SEED)) LIST$SEED <- 123
  stopifnot(is.numeric(LIST$SEED))
  stopifnot(LIST$SEED>0)
  out$SEED <- as.integer(LIST$SEED)

  if(is.null(LIST$VERBOSE) & SETTING=="init") LIST$VERBOSE <- FALSE
  if(is.null(LIST$VERBOSE) & SETTING=="main") LIST$VERBOSE <- TRUE
  stopifnot(is.logical(LIST$VERBOSE))
  out$VERBOSE <- LIST$VERBOSE


  return(out)

}
