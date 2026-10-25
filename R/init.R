init_thresholds <- function(DIMS, CONSTR_LIST){
  out <- c()
  s <- 0
  for (i in 1:length(DIMS$cat)) {
    if(DIMS$cat[i]==2){
      out[(s + 1):(s + DIMS$cat[i] - 1)] <- 0
    }else{
      vec <- 1:(DIMS$cat[i]-1)
      vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
      out[(s + 1):(s + DIMS$cat[i] - 1)] <- vec
    }
    s <- s + DIMS$cat[i] - 1
  }
  return(out)
}

init_loadings <- function(DIMS, CONSTR_LIST, FXD=0.25){
  out <-  rep(FXD, sum(is.na(CONSTR_LIST$CONSTRMAT)))
  return(out)
}

init_transformed_latcorr <- function(DIMS, CONSTR_LIST){

  out <- NULL
  if(CONSTR_LIST$CORRFLAG==1) out <- rep(0, DIMS$q*(DIMS$q-1)/2)
  return(out)
}

init_transformed_latsd <- function(DIMS, CONSTR_LIST){
  rep(0, sum(is.na(CONSTR_LIST$CONSTRLOGSD)))
}

init_par_with_plSA <- function(DIMS, CONSTR_LIST, FREQ, INIT_CONTROL, INIT){

  start_par <- INIT
  if(is.null(INIT_CONTROL$STEP0)) INIT_CONTROL$STEP0 <- sqrt(1/DIMS$d)

  cpp_args <- c(list(FREQ        = FREQ,
                     VALFREQ     = FREQ,
                     N           = DIMS$n,
                     C_VEC       = DIMS$cat,
                     CONSTRMAT   = CONSTR_LIST$CONSTRMAT,
                     CONSTRLOGSD = CONSTR_LIST$CONSTRLOGSD,
                     LLC         = CONSTR_LIST$LLC,
                     THETA_INIT  = start_par,
                     NTHR        = DIMS$nthr,
                     NLOAD       = DIMS$nload,
                     NCORR       = DIMS$ncorr,
                     NVAR        = DIMS$nvar
  ), unclass(INIT_CONTROL), VERBOSE = FALSE)

  sto_init <- do.call(cpp_plSA2, cpp_args)


  return(sto_init$avtheta)
}
