#' Compute pairwise frequencies
#'
#' It returns a 5-rows matrix, with one column for each combination of items and
#' categories.
#' Row0: item k, Row1: item l, Row2; category item k, Row3: category item l, Row4: freq
#'
#' @param Y Integer matrix of dimension \eqn{n*p}, where \eqn{n} is the sample size
#' and \eqn{p} is the number of items considered. Categories must be coded starting from zero.
#' For example, an item with three categories can only accept values contained in
#' \eqn{{0, 1, 2}}.
#' @param C_VEC Integer vector indicating how many possible categories are associated to
#' each item in 'Y'.
#'
#' @export
compute_frequencies <- function(Y, C_VEC){
  if(sum(!is.finite(C_VEC))>0 | is.matrix(C_VEC))stop('C_VEC is not a numeric vector.')
  if(length(C_VEC) != ncol(Y))stop('Y and C_VEC dimensions do not match.')

  pairs_freq(Y = Y, C_VEC = C_VEC)
}

#' Automatic construction of simple constraints matrices
#'
#' build_constrMat() creates a binary matrix with entries equal to 1
#' when referring to free loadings, and 0 entries otherwise.
#'
#' @param P Number of items.
#' @param Q Number of latent vaiables.
#' @param STRUCT Label for the constraint structure chosen. Only two possible values:
#' 'triangular', 'simple' and 'crossed'.
#' @param CROSS Integer referring to the number of items latent i shares with latent i+1
#'
#' @export
build_constrMat <- function(P, Q, STRUCT = c('simple', 'triangular', 'crossed'), CROSS = NULL){
  p <- check_p(P)
  q <- check_q(Q)

  STRUCT <- match.arg(STRUCT)
  if(STRUCT == 'triangular'){
    # Build lower trinagular matrix
    constrMat <- matrix(NA, nrow = p, ncol = q)
    for(j in 1:p){
      for (h in 1:q) {
        if(h > j) constrMat[j,h] <- 0
      }
    }

  } else if ( STRUCT == 'simple'){
    # Each item to one factor
    loadings_per_factor <- p %/% q
    remaining_loadings <- p %% q

    constrMat <- matrix(0, nrow = p, ncol = q)
    for (h in 1:(q-1)) {
      constrMat[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor), h] <- rep(NA,loadings_per_factor)
    }
    h <- q
    constrMat[ (loadings_per_factor*(h-1)+1) : (loadings_per_factor*(h-1) + loadings_per_factor + remaining_loadings), h] <- rep(NA,loadings_per_factor + remaining_loadings)
  } else if( STRUCT == 'crossed'){
    if(is.null(CROSS)) CROSS <- 1

    # Each item to one factor
    loadings_per_factor <- p %/% q
    remaining_loadings <- p %% q

    constrMat <- matrix(0, nrow = p, ncol = q)
    for (h in 1:(q-1)) {
      constrMat[ max(1, (loadings_per_factor*(h-1) + 1 - CROSS)) : (loadings_per_factor*(h-1) + loadings_per_factor), h] <- NA
    }
    h <- q
    constrMat[ max(1, (loadings_per_factor*(h-1) + 1 - CROSS)) : (loadings_per_factor*(h-1) + loadings_per_factor + remaining_loadings), h] <- NA

  }



  return(constrMat)
}


#' Construct loading matrix
#'
#' @param FIXED Fixed value to assign to all free loadings. If 'NULL' it draws
#' them randomly from Unif(0,1)
#' @param CONSTRMAT Binary matrix of dimension \eqn{p*q}. A cell equal to \eqn{1} indicates
#' that the corresponding element in the loading matrix is free to be estimated.
#' A cell equal to \eqn{0} fixes the corresponding element in the loading matrix
#' to \eqn{0}.
#' @param LLC Linear loadings constraints. Expects a list of constraints. See [fit_plFA] documentation.
#' @param SEED Random seed.
#' @param LB Lower bound for uniform random generator. Default set to 0.
#' @param UB Upper bound for uniform random generator. Default set to 1.
#' @param STDLV  Logical indicator. Set it to `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.
#' @export
gen_loadings <- function(CONSTRMAT, FIXED = NULL,  SEED = 123, LB = 0, UB = 1, STDLV=TRUE, LLC=NULL){
  CONSTRMAT <- check_cnstr_loadings(CONSTRMAT, STDLV = STDLV, LLC=LLC)

  set.seed(SEED)
  p <- check_p(nrow(CONSTRMAT))
  q <- check_q(ncol(CONSTRMAT))

  out <-  CONSTRMAT
  for (j in 1:p) {
    for (h in 1:q) {
      if(is.na(out[j,h])) out[j,h] <- runif(1, LB, UB)
    }
  }

  if(!is.null(FIXED)){
    for (j in 1:p) {
      for (h in 1:q) {
        if(is.na(out[j,h])) out[j,h] <- FIXED
      }
    }
  }

  if(!is.null(LLC)){
    for (idx_lc in 1:length(LLC)) {
      val <- 0
      for (idx_var in 2:length(LLC[[idx_lc]])) {
        val <- val + LLC[[idx_lc]][[idx_var]][1]*out[LLC[[idx_lc]][[idx_var]][2], LLC[[idx_lc]][[idx_var]][3]]
      }
      out[LLC[[idx_lc]][[1]][1], LLC[[idx_lc]][[1]][2]] <- val
    }
  }


  return(out)
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

  n          <- check_n(SAMPLE_SIZE)
  loadings   <- check_loadings(LOADINGS)
  thresholds <- check_thresholds_j(THRESHOLDS)
  lat_cov    <- check_latcov(LATENT_COV)
  stopifnot(ncol(loadings)==ncol(lat_cov))

  set.seed(SEED)

  p <- check_p(nrow(loadings))
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
