#' Automatic construction of simple constraints matrices
#'
#' build_constr_mat() creates a binary matrix with entries equal to 1
#' when referring to free loadings, and 0 entries otherwise.
#'
#' @param P Number of items.
#' @param Q Number of latent vaiables.
#' @param STRUCT Label for the constraint structure chosen. Only two possible values:
#' 'triangular', 'simple' and 'crossed'.
#' @param CROSS Integer referring to the number of items latent i shares with latent i+1
#'
#' @export
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


#' Construct loading matrix
#'
#' @param FIXED Fixed value to assign to all free loadings. If 'NULL' it draws
#' them randomly from Unif(0,1)
#' @param CONSTRMAT Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
#' that the corresponding element in the loading matrix is free to be estimated.
#' A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
#' to \eqn{0}.
#' @param SEED Random seed.
#' @param LB Lower bound for uniform random generator. Default set to 0.
#' @param UB Upper bound for uniform random generator. Default set to 1.
#' @param STDLV  Logical indicator. Set it to `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.
#' @export
gen_loadings <- function(CONSTRMAT, FIXED = NULL,  SEED = 123, LB = 0, UB = 1, STDLV=TRUE){
  out <- t(set_latpl_model(CONSTRMAT, STDLV = STDLV)$CONSTRMAT)
  set.seed(SEED)
  out[is.na(out)] <- if (is.null(FIXED)) runif(sum(is.na(out)), LB, UB) else FIXED
  t(out)
}

#' Simulate new data
#'
#' Simulate data given the values of the parameters for thresholds, loadings and latent correlations
#'
#' @param SAMPLE_SIZE Number of observations to simulate.
#' @param LOADINGS Matrix of loadings, with dimensions \eqn{p*q}.
#' @param THRESHOLDS Vector of thresholds fo a single item.
#' It is assumed to be the same for all the items.
#' @param LATENT_COV Latent correlation matrix. Dimension \eqn{q*q}.
#' @param SEED Numeric seed for random generation.
#'
#'@export
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
