build_constr_mat <- function(P, Q, STRUCT = c('simple', 'triangular', 'crossed'), CROSS = NULL){
  stopifnot(P > 1, Q > 0)
  p <- as.integer(P)
  q <- as.integer(Q)

  STRUCT <- match.arg(STRUCT)
  if(STRUCT == 'triangular'){
    # Build lower trinagular matrix
    out <- matrix(NA, nrow = p, ncol = q)
    for(j in 1:p){
      for (h in 1:q) {
        if(h > j) out[j,h] <- 0
      }
    }

  } else if ( STRUCT == 'simple'){
    # Each item to one factor
    loadings_per_factor <- p %/% q
    remaining_loadings <- p %% q

    out <- matrix(0, nrow = p, ncol = q)
    for (h in 1:(q-1)) {
      out[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor), h] <- rep(NA,loadings_per_factor)
    }
    h <- q
    out[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor + remaining_loadings), h] <- rep(NA,loadings_per_factor + remaining_loadings)
  } else if( STRUCT == 'crossed'){
    if(is.null(CROSS)) CROSS <- 1

    # Each item to one factor
    loadings_per_factor <- p %/% q
    remaining_loadings <- p %% q

    out <- matrix(0, nrow = p, ncol = q)
    for (h in 1:(q-1)) {
      out[ max(1, (loadings_per_factor*(h-1) + 1 - CROSS)) : (loadings_per_factor*(h-1) + loadings_per_factor), h] <- NA
    }
    h <- q
    out[ max(1, (loadings_per_factor*(h-1) + 1 - CROSS)) : (loadings_per_factor*(h-1) + loadings_per_factor + remaining_loadings), h] <- NA

  }

  return(out)
}

gen_loadings <- function(CONSTRMAT, FIXED = NULL, SEED = 123, LB = 0, UB = 1){
  out <- t(CONSTRMAT)
  set.seed(SEED)
  out[is.na(out)] <- if (is.null(FIXED)) runif(sum(is.na(out)), LB, UB) else FIXED
  t(out)
}

sim_data <- function(SAMPLE_SIZE, LOADINGS, THRESHOLDS, LATENT_COV, SEED = 123){

  stopifnot(SAMPLE_SIZE > 0, is.matrix(LOADINGS), all(is.finite(LOADINGS)), nrow(LOADINGS) > ncol(LOADINGS),
            !is.unsorted(THRESHOLDS), isSymmetric(LATENT_COV), ncol(LOADINGS) == ncol(LATENT_COV))
  n <- as.integer(SAMPLE_SIZE)
  loadings <- LOADINGS
  thresholds <- THRESHOLDS
  lat_cov <- LATENT_COV

  set.seed(SEED)

  p <- nrow(loadings)
  error_variance <- diag(1, p, p) - diag(diag(loadings%*%lat_cov%*%t(loadings)),p,p)

  ev <- eigen(error_variance, symmetric = TRUE)
  errors <- matrix(rnorm(n * p), nrow = n, byrow = TRUE) %*% t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  ev <- eigen(lat_cov, symmetric = TRUE)
  factors <- matrix(rnorm(n * ncol(lat_cov)), nrow = n, byrow = TRUE) %*% t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))

  URV <- factors %*% t(loadings) + errors

  out <- matrix(0, n, p)
  for(i in 1:length(thresholds)){
    out[URV>thresholds[i]] <- i
  }

  return(out)
}
