
#' Extract latent covariance matrix from unconstrained parameter vector
#'
#' @param THETA Unconstrained parameter vector
#' @param CONSTRLOGSD  \eqn{Q}-dimensional vector. Elements set to `NA` refers to free latent log standard deviations parameters. Elements set to numerical values denote fixed values constraints.
#' @param NTHR Number of thresholds parameters.
#' @param NLOAD Number of free loadings parameters
#' @param NCORR Number of free latent correlations parameters.
#' @param NVAR Number of free latent variance parameters.
#' @param Q Number of latent variables
#'
#' @return
#' Latent covariance matrix of dimension  \eqn{Q*Q}.
#'
#' @export
get_S <- function(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q){
  THETA <- check_theta(THETA)
  return(cpp_latvar_theta2mat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, Q))
}


#' Extract latent correlation matrix from unconstrained parameter vector
#'
#' @param THETA Unconstrained parameter vector
#' @param NTHR Number of thresholds parameters.
#' @param NLOAD Number of free loadings parameters
#' @param NCORR Number of free latent correlations parameters.
#' @param NVAR Number of free latent variance parameters.
#' @param Q Number of latent variables
#'
#' @return Latent correlation matrix
#'
#' @export
get_R <- function(THETA, NTHR, NLOAD, NCORR, NVAR, Q){
  Q <- check_q(Q)
  THETA <- check_theta(THETA)
  R <- cpp_latvar_theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, Q)
  return(R)
}

#' Extract vector of latent correlations from unconstrained parameter vector
#'
#' @param THETA Unconstrained parameter vector
#' @param NTHR Number of thresholds parameters.
#' @param NLOAD Number of free loadings parameters
#' @param NCORR Number of free latent correlations parameters.
#' @param NVAR Number of free latent variance parameters.
#' @param Q Number of latent variables
#'
#' @return Vector of latent correlations
#'
#' @export
#' @export
get_corr <- function(THETA, NTHR, NLOAD, NCORR, NVAR, Q){
  Q <- check_q(Q)
  THETA <- check_theta(THETA)

  R <- cpp_latvar_theta2cmat(THETA, NTHR, NLOAD, NCORR, NVAR, Q)
  out <- R[upper.tri(R)]

  return(out)
}

#' Extract free loadings from parameter vector
#'
#' @param THETA Unconstrained parameter vector
#' @param NTHR Number of thresholds parameters.
#' @param NLOAD Number of free loadings parameters
#'
#'@export
get_lambda <- function(THETA, NTHR, NLOAD){
  THETA <- check_theta(THETA)

  lambda <- cpp_loadings_theta2vec(THETA, NTHR, NLOAD)

  return(lambda)
}

#' Get full parameter vector
#'
#' @description
#' Assemble the parameter vector starting from the threshold vector,
#' the latent covariance matrix and the loading matrix
#'
#' @param THRESHOLDS Numerical vector of thresholds parameters.
#' If all the \eqn{p} items have \eqn{c} possible categories, it contains
#' \eqn{p(c-1)} elements.
#' @param LOADINGS Numerical matrix of dimension \eqn{p*q},
#' where \eqn{q} is the number of latent variables.
#' @param LATENT_COV Latent covariance matrix of dimension \eqn{q*q}.
#' @param CAT Integer vector storing the possible number of categories
#' for each item. Values must be ordered following items order in the dataset.
#' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA` refers to free loading parameters. Elements set to numerical values denote fixed values constraints.
#' @param CONSTRVAR \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent variance parameters. Elements set to numerical values denote fixed values constraints.
#' @param CORRFLAG Logical indicator. Set it to `FALSE` if the latent variables are independent. Set it `TRUE` otherwise.
#' @param STDLV Logical indicator. Set it to `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.

#' @export
get_theta <- function(THRESHOLDS, LOADINGS, LATENT_COV, CAT, CONSTRMAT, CONSTRVAR, CORRFLAG, STDLV){

  CONSTRMAT <- check_cnstr_loadings(CONSTRMAT, STDLV=STDLV)
  CONSTRLOGSD <- check_cnstr_latvar(CONSTRVAR, Q=ncol(CONSTRMAT), STDLV=STDLV)
  THRESHOLDS <- check_thresholds(THRESHOLDS, CAT)
  LOADINGS <- check_loadings(LOADINGS)
  LATENT_COV <- check_latcov(LATENT_COV)


  stopifnot(dim(LOADINGS)==dim(CONSTRMAT))
  stopifnot(ncol(LOADINGS)==ncol(LATENT_COV))
  stopifnot(length(THRESHOLDS)==(sum(CAT)-nrow(CONSTRMAT)))



  load_vec <- cpp_loadings_mat2vec(LOADINGS, CONSTRMAT, sum(is.na(CONSTRMAT)))

  q <- ncol(CONSTRMAT)
  ncorr <- if(CORRFLAG) q*(q-1)/2 else 0

  lat_vec <- cpp_latvar_mat2vec(S=LATENT_COV,
                                CONSTRLOGSD = CONSTRLOGSD,
                                NCORR = ncorr,
                                NVAR = sum(is.na(CONSTRLOGSD)))

  theta <- c(THRESHOLDS, load_vec, lat_vec)
  return(theta)
}
