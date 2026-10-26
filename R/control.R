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
check_cnstr_loadings <- function(MAT, STDLV){
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
  ls$CONSTRAINTS <- LIST$CONSTRAINTS
  ls$CONSTRMAT <- check_cnstr_loadings(LIST$CONSTRMAT, STDLV = ls$STDLV)
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
                              NTHR      = DIMS$nthr,
                              NLOAD     = DIMS$nload)
  check_loadings(L)
  S <- cpp_latvar_theta2mat(THETA       = PAR,
                            CONSTRLOGSD = CONSTR_LIST$CONSTRLOGSD,
                            NTHR        = DIMS$nthr,
                            NLOAD       = DIMS$nload,
                            NCORR       = DIMS$ncorr,
                            NVAR        = DIMS$nvar,
                            Q           = DIMS$q,
                            CORRFLAG    = corrflag_code(CONSTR_LIST))
  check_latcov(S)
  Sy <- L%*%S%*%t(L)
  if(!all(diag(Sy)<=1)) warning("Negative error variances implied by init parameters")
  if (corrflag_code(CONSTR_LIST) == 2) return(PAR)
  PAR <- cpp_sa_proj( THETA=PAR,
                      CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                      CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
                      C_VEC=DIMS$cat,
                      CORRFLAG=CONSTR_LIST$CORRFLAG,
                      NTHR=DIMS$nthr,
                      NLOAD=DIMS$nload,
                      NCORR=DIMS$ncorr,
                      NVAR=DIMS$nvar)

  return(PAR)
}
check_init <- function(INIT, FREQ, DIMS, CONSTR_LIST, INIT_METHOD=c("SA", "custom", "standard"), VERBOSE=FALSE, INIT_CONTROL=set_sa_init_control(), UFREQ=numeric(0)){

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
    if (corrflag_code(CONSTR_LIST) != 2) INIT <- cpp_sa_proj(THETA=INIT,
                        CONSTRMAT=CONSTR_LIST$CONSTRMAT,
                        CONSTRLOGSD=CONSTR_LIST$CONSTRLOGSD,
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
    INIT <- init_par_with_plSA(DIMS, CONSTR_LIST, FREQ, INIT_CONTROL, INIT, UFREQ)
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
#' Ordinal data for pairwise likelihood estimation
#'
#'@param DATA Integer data matrix of dimension \eqn{n*p}. Categories must be
#'  coded starting from zero.
#'@param VALDATA Validation data used to monitor convergence of stochastic
#'  approximations. If `NULL`, data passed via `DATA` is used for monitoring
#'  purposes.
#'@param MISSING Treatment of missing values (`NA`): `"listwise"` drops incomplete observations,
#'  `"pairwise"` uses the complete pairs and `"available.cases"` the available cases of
#'  Katsikatsou et al. (2022).
#'@param WEIGHTS Sampling weights, one per observation, normalised to sum to the number of
#'  observations. If `NULL`, unit weights.
#'@param WEIGHTS_TYPE Type of sampling weights, as in lavaan: `"design"` or `"frequency"`.
#'@param CLUSTER Cluster membership, one per observation, for cluster-robust standard errors. If
#'  `NULL`, observations are independent.
#'
#' @return An object of class `plfa_data` with the data matrix `Y`, its dimensions `n` and `p`,
#'   the number of categories of each item `cat`, the pairwise frequency tables `freq` and
#'   `valfreq`, the treatment of missing values `missing`, the unit weights `uweights` and
#'   univariate frequencies `ufreq` of the available-cases term, the normalised sampling weights
#'   `weights` and their factors on the scores `jweights`, the cluster indices `cluster`, and the
#'   time spent computing the frequencies `freqTime`.
#'
#' @export
set_plfa_data <- function(DATA, VALDATA = NULL, MISSING = c("listwise", "pairwise", "available.cases"),
                          WEIGHTS = NULL, WEIGHTS_TYPE = c("design", "frequency"), CLUSTER = NULL) {
  MISSING <- match.arg(MISSING)
  WEIGHTS_TYPE <- match.arg(WEIGHTS_TYPE)
  Y <- check_data(DATA)
  stopifnot(
    is.null(WEIGHTS) || (is.numeric(WEIGHTS) && length(WEIGHTS) == nrow(Y) && all(is.finite(WEIGHTS)) && all(WEIGHTS >= 0)),
    is.null(CLUSTER) || (length(CLUSTER) == nrow(Y) && !anyNA(CLUSTER))
  )
  if (MISSING == "listwise" && anyNA(Y)) {
    keep <- stats::complete.cases(Y)
    message(sum(!keep), " incomplete observations removed (MISSING = \"listwise\").")
    Y <- Y[keep, , drop = FALSE]
    if (!is.null(WEIGHTS)) WEIGHTS <- WEIGHTS[keep]
    if (!is.null(CLUSTER)) CLUSTER <- CLUSTER[keep]
  }
  stopifnot(nrow(Y) > 0, colSums(!is.na(Y)) > 0)
  cat <- apply(Y, 2, max, na.rm = TRUE) + 1
  stopifnot(cat > 1)
  weights <- jweights <- cluster <- numeric(0)
  if (!is.null(WEIGHTS)) {
    weights <- WEIGHTS / sum(WEIGHTS) * nrow(Y)
    jweights <- if (WEIGHTS_TYPE == "design") weights else sqrt(weights)
  }
  if (!is.null(CLUSTER)) {
    cluster <- as.numeric(match(CLUSTER, unique(CLUSTER)) - 1)
    stopifnot(length(unique(cluster)) > 1)
  }
  freq_start_time <- Sys.time()
  freq <- pairs_freq(Y, cat, weights)
  uweights <- if (MISSING == "available.cases") rowSums(is.na(Y)) else rep(0, nrow(Y))
  uw <- if (length(weights)) uweights * weights else uweights
  ufreq <- unlist(lapply(seq_len(ncol(Y)), function(j) {
    vapply(seq_len(cat[j]) - 1, function(a) sum(uw[Y[, j] %in% a]), numeric(1))
  }))
  valfreq <- NULL
  if (!is.null(VALDATA)) {
    stopifnot(
      is.matrix(VALDATA),
      is.numeric(VALDATA),
      all(is.finite(VALDATA)),
      ncol(VALDATA) == ncol(Y),
      all(VALDATA >= 0),
      all(t(VALDATA) <= cat - 1)
    )
    valfreq <- pairs_freq(round(VALDATA), cat)
  }
  out <- list(
    Y = Y,
    n = nrow(Y),
    p = ncol(Y),
    cat = cat,
    freq = freq,
    valfreq = valfreq,
    missing = MISSING,
    uweights = uweights,
    ufreq = ufreq,
    weights = weights,
    jweights = jweights,
    cluster = cluster,
    freqTime = as.numeric(difftime(Sys.time(), freq_start_time, units = "secs"))
  )
  class(out) <- "plfa_data"
  out
}

#' Factor model for pairwise likelihood estimation
#'
#' @param CONSTRMAT \eqn{p*q}-dimensional matrix. Elements set to `NA`
#'  refers to free loading parameters. Elements set to numerical values denote
#'  fixed values constraints.
#' @param CONSTRVAR \eqn{q}-dimensional vector. Elements set to `NA` refers to free latent
#'  variance parameters. Elements set to numerical values denote fixed values
#'  constraints.
#' @param CORRFLAG Logical indicator. Set it to
#'  `FALSE` if the latent variables are assumed to be independent. Set it `TRUE`
#'  otherwise.
#' @param STDLV Logical indicator. Set it to
#'  `TRUE` to fix latent variables scale. Set it `FALSE` to fix loadings scale.
#' @param CONSTRAINTS Linear equality constraints `A %*% theta = b` on the parameter vector, as
#'  `list(A = , b = )`, with one column of `A` per parameter in the order of [theta_names()]. Each
#'  constraint involves either thresholds, loadings, latent covariances (with `LATCOV = "direct"`) or
#'  latent log standard deviations.
#' @param LATCOV Parametrisation of the latent correlations: `"lkj"` (Cholesky factor with Fisher-z
#'  angles, always positive definite) or `"direct"` (the latent covariances themselves, which allows
#'  linear constraints on them).
#'
#' @return An object of class `plfa_model` with the identified loading constraints `CONSTRMAT`,
#'   the log standard deviations of the latent variables `CONSTRLOGSD` (`NA` when free),
#'   `CORRFLAG`, `STDLV`, `CONSTRAINTS`, `LATCOV`, and the numbers of items `p` and latent variables `q`.
#'
#' @export
set_plfa_model <- function(CONSTRMAT, CONSTRVAR = NULL, CORRFLAG = TRUE, STDLV = TRUE, CONSTRAINTS = NULL,
                           LATCOV = c("lkj", "direct")) {
  LATCOV <- match.arg(LATCOV)
  if (is.matrix(CONSTRMAT) && all(is.na(CONSTRMAT))) storage.mode(CONSTRMAT) <- "double"
  stopifnot(is.matrix(CONSTRMAT))
  if (!is.null(CONSTRAINTS)) {
    stopifnot(
      is.list(CONSTRAINTS), is.matrix(CONSTRAINTS$A), is.numeric(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$A)),
      is.numeric(CONSTRAINTS$b), length(CONSTRAINTS$b) == nrow(CONSTRAINTS$A), all(is.finite(CONSTRAINTS$b))
    )
  }
  out <- check_cnstr(list(CONSTRMAT = CONSTRMAT, CONSTRVAR = CONSTRVAR, CORRFLAG = CORRFLAG, STDLV = STDLV,
                          CONSTRAINTS = CONSTRAINTS))
  if (LATCOV == "direct" && !out$CORRFLAG) stop("LATCOV = \"direct\" requires CORRFLAG = TRUE.")
  out$LATCOV <- LATCOV
  out$p <- nrow(out$CONSTRMAT)
  out$q <- ncol(out$CONSTRMAT)
  class(out) <- "plfa_model"
  out
}

#' Names of the model parameters
#'
#' @param DATA Data object from [set_plfa_data()].
#' @param MODEL Model object from [set_plfa_model()].
#' @return The names of the parameters in the order used by [fit_plFA()]: thresholds `tau[i,k]`,
#'   free loadings `lambda[j,h]`, latent correlation parameters `eta[r,c]` and free latent log
#'   standard deviations `logsd[h]`.
#' @export
theta_names <- function(DATA, MODEL) {
  stopifnot(inherits(DATA, "plfa_data"), inherits(MODEL, "plfa_model"))
  idx <- which(is.na(MODEL$CONSTRMAT), arr.ind = TRUE)
  q <- MODEL$q
  c(
    unlist(lapply(seq_len(DATA$p), function(i) sprintf("tau[%d,%d]", i, seq_len(DATA$cat[i] - 1)))),
    sprintf("lambda[%d,%d]", idx[, 1], idx[, 2]),
    if (MODEL$CORRFLAG && q > 1) unlist(lapply(2:q, function(h) {
      sprintf(if (corrflag_code(MODEL) == 2) "psi[%d,%d]" else "eta[%d,%d]", seq_len(h - 1), h)
    })),
    sprintf("logsd[%d]", which(is.na(MODEL$CONSTRLOGSD)))
  )
}

corrflag_code <- function(MODEL) {
  if (!MODEL$CORRFLAG) 0L else if (identical(MODEL$LATCOV, "direct")) 2L else 1L
}

lincon <- function(CONSTRAINTS, DIMS, DIRECT = FALSE) {
  d <- DIMS$d
  if (is.null(CONSTRAINTS) || nrow(CONSTRAINTS$A) == 0) return(list(K = diag(d), c = rep(0, d), free_idx = seq_len(d)))
  A <- CONSTRAINTS$A
  b <- CONSTRAINTS$b
  if (ncol(A) != d) stop("CONSTRAINTS$A must have one column per parameter (", d, "), in the order of theta_names().")
  block <- rep(1:4, c(DIMS$nthr, DIMS$nload, DIMS$ncorr, DIMS$nvar))
  for (r in seq_len(nrow(A))) {
    blk <- unique(block[A[r, ] != 0])
    if (any(blk == 3) && !DIRECT) stop("Constraints on the latent correlations require LATCOV = \"direct\".")
    if (length(blk) > 1) stop("Each constraint must involve either thresholds, loadings, latent covariances or latent log standard deviations.")
  }
  piv <- integer(0)
  for (j in seq_len(d)) {
    r <- length(piv) + 1
    if (r > nrow(A)) break
    i <- r - 1 + which.max(abs(A[r:nrow(A), j]))
    if (abs(A[i, j]) < 1e-10) next
    A[c(r, i), ] <- A[c(i, r), ]
    b[c(r, i)] <- b[c(i, r)]
    b[r] <- b[r] / A[r, j]
    A[r, ] <- A[r, ] / A[r, j]
    b[-r] <- b[-r] - A[-r, j] * b[r]
    A[-r, ] <- A[-r, , drop = FALSE] - outer(A[-r, j], A[r, ])
    piv <- c(piv, j)
  }
  if (any(abs(b[-seq_along(piv)]) > 1e-8)) stop("CONSTRAINTS are inconsistent.")
  free_idx <- setdiff(seq_len(d), piv)
  K <- matrix(0, d, length(free_idx))
  K[cbind(free_idx, seq_along(free_idx))] <- 1
  cvec <- rep(0, d)
  for (r in seq_along(piv)) {
    K[piv[r], ] <- -A[r, free_idx]
    cvec[piv[r]] <- b[r]
  }
  list(K = K, c = cvec, free_idx = free_idx)
}

#' Starting values for pairwise likelihood estimation
#'
#' @param METHOD \tabular{ll}{ \code{'standard'} \tab Uses cold
#'  initialisation where loadings are initialised at `0.25` and the latent
#'  covariance matrix as an Identity matrix. \cr \tab \cr \code{'SA'} \tab
#'  Computes the starting point using a short chain of stochastic updates.
#'  Updates start from the `standard` initialisation point. \cr }
#'  With fewer than 15 items, `'SA'` falls back to `'standard'`.
#' @param THETA Initialising vector. If `NULL`, the `standard` initialisation point is used.
#' @param CONTROL Control object from [set_sa_init_control()], used when `METHOD = "SA"`.
#'
#' @return An object of class `plfa_init` to be passed to [fit_plFA()] via `INIT`, or to [cfa()]
#'   via `estimator.args$init`.
#'
#' @export
set_init <- function(METHOD = c("standard", "SA"), THETA = NULL, CONTROL = set_sa_init_control()) {
  METHOD <- match.arg(METHOD)
  stopifnot(
    is.null(THETA) || (is.numeric(THETA) && is.vector(THETA) && all(is.finite(THETA))),
    inherits(CONTROL, "plfa_sa_control")
  )
  out <- list(METHOD = METHOD, THETA = THETA, CONTROL = CONTROL)
  class(out) <- "plfa_init"
  out
}

#' Control options for the ucminf optimiser
#'
#' @param GRTOL Tolerance on the largest absolute gradient element.
#' @param XTOL Tolerance on the relative step length.
#' @param STEPMAX Initial maximal step length of the trust region.
#' @param MAXEVAL Maximum number of function evaluations.
#'
#' @return An object of class `plfa_ucminf_control` to be passed to [fit_plFA()] via `CONTROL`
#'   when `METHOD = "ucminf"`, or to [cfa()] via `control`.
#'
#' @export
set_ucminf_control <- function(
  GRTOL = 1e-6,
  XTOL = 1e-12,
  STEPMAX = 1,
  MAXEVAL = 500
) {
  stopifnot(
    is.numeric(GRTOL), length(GRTOL) == 1, GRTOL > 0,
    is.numeric(XTOL), length(XTOL) == 1, XTOL > 0,
    is.numeric(STEPMAX), length(STEPMAX) == 1, STEPMAX > 0,
    is.numeric(MAXEVAL), length(MAXEVAL) == 1, MAXEVAL >= 1
  )
  out <- list(GRTOL = GRTOL, XTOL = XTOL, STEPMAX = STEPMAX, MAXEVAL = as.integer(MAXEVAL))
  class(out) <- "plfa_ucminf_control"
  out
}

#' Control options for the pairwise stochastic approximation
#'
#' @param PAIRS_PER_ITERATION Number of item pairs drawn at each iteration.
#' @param STEP0 Initial step length.
#' @param STEP1 Scaling of the step length.
#' @param STEP2 Rate of the step length decay.
#' @param STEP3 Exponent of the step length decay.
#' @param BURNE Maximum number of burn-in epochs before trajectory averaging.
#' @param MAXE Maximum number of epochs.
#' @param EPOCHS_PER_CHECK Number of epochs between evaluations of the full negative pairwise
#'   log-likelihood.
#' @param CHECK_TOL If `TRUE`, stop the burn-in and the optimisation when the relative change in
#'   the full negative pairwise log-likelihood falls below `TOL_BURN` and `TOL_END`.
#' @param MAX_TOL_COUNTER Number of consecutive checks below tolerance required to stop.
#' @param TOL_BURN Tolerance ending the burn-in.
#' @param TOL_END Tolerance ending the optimisation.
#' @param SEED Seed of the pair sampling.
#'
#' @return An object of class `plfa_sa_control` to be passed to [fit_plFA()] via `CONTROL` when
#'   `METHOD = "SA"`, or to [cfa()] via `control`.
#'
#' @export
set_sa_control <- function(
  PAIRS_PER_ITERATION = 16,
  STEP0 = 1,
  STEP1 = 1,
  STEP2 = 1e-8,
  STEP3 = 0.75,
  BURNE = 20,
  MAXE = 30,
  EPOCHS_PER_CHECK = 5,
  CHECK_TOL = TRUE,
  MAX_TOL_COUNTER = 1,
  TOL_BURN = 5e-7,
  TOL_END = 1e-7,
  SEED = 123
) {
  stopifnot(
    is.numeric(PAIRS_PER_ITERATION), length(PAIRS_PER_ITERATION) == 1, PAIRS_PER_ITERATION >= 1,
    is.null(STEP0) || (is.numeric(STEP0) && length(STEP0) == 1 && STEP0 > 0),
    is.numeric(STEP1), length(STEP1) == 1, STEP1 > 0,
    is.numeric(STEP2), length(STEP2) == 1, STEP2 > 0,
    is.numeric(STEP3), length(STEP3) == 1, STEP3 > 0,
    is.numeric(BURNE), length(BURNE) == 1, BURNE >= 0,
    is.numeric(MAXE), length(MAXE) == 1, MAXE >= 0,
    BURNE < MAXE,
    is.numeric(EPOCHS_PER_CHECK), length(EPOCHS_PER_CHECK) == 1, EPOCHS_PER_CHECK >= 1,
    is.logical(CHECK_TOL), length(CHECK_TOL) == 1,
    is.numeric(MAX_TOL_COUNTER), length(MAX_TOL_COUNTER) == 1, MAX_TOL_COUNTER >= 1,
    is.numeric(TOL_BURN), length(TOL_BURN) == 1, TOL_BURN > 0,
    is.numeric(TOL_END), length(TOL_END) == 1, TOL_END > 0,
    is.numeric(SEED), length(SEED) == 1, SEED > 0
  )
  out <- list(
    PAIRS_PER_ITERATION = as.integer(PAIRS_PER_ITERATION),
    STEP0 = STEP0,
    STEP1 = STEP1,
    STEP2 = STEP2,
    STEP3 = STEP3,
    BURNE = as.integer(BURNE),
    MAXE = as.integer(MAXE),
    EPOCHS_PER_CHECK = as.integer(EPOCHS_PER_CHECK),
    CHECK_TOL = CHECK_TOL,
    MAX_TOL_COUNTER = as.integer(MAX_TOL_COUNTER),
    TOL_BURN = TOL_BURN,
    TOL_END = TOL_END,
    SEED = as.integer(SEED)
  )
  class(out) <- "plfa_sa_control"
  out
}

#' Control options for the stochastic approximation initialisation
#'
#' @inheritParams set_sa_control
#' @param STEP0 Initial step length. If `NULL`, it is set to \eqn{\sqrt{1/d}} at fit time, where
#'   \eqn{d} is the number of free parameters.
#'
#' @return An object of class `plfa_sa_control` to be passed to [fit_plFA()] via `INIT_CONTROL`,
#'   or to [cfa()] via `estimator.args$init_control`, when the initialisation method is `"SA"`.
#'
#' @export
set_sa_init_control <- function(
  PAIRS_PER_ITERATION = 16,
  STEP0 = NULL,
  STEP1 = 1,
  STEP2 = 1e-8,
  STEP3 = 0.75,
  BURNE = 3,
  MAXE = 5,
  EPOCHS_PER_CHECK = 1,
  CHECK_TOL = FALSE,
  MAX_TOL_COUNTER = 1,
  TOL_BURN = 5e-7,
  TOL_END = 1e-7,
  SEED = 123
) {
  set_sa_control(
    PAIRS_PER_ITERATION = PAIRS_PER_ITERATION,
    STEP0 = STEP0,
    STEP1 = STEP1,
    STEP2 = STEP2,
    STEP3 = STEP3,
    BURNE = BURNE,
    MAXE = MAXE,
    EPOCHS_PER_CHECK = EPOCHS_PER_CHECK,
    CHECK_TOL = CHECK_TOL,
    MAX_TOL_COUNTER = MAX_TOL_COUNTER,
    TOL_BURN = TOL_BURN,
    TOL_END = TOL_END,
    SEED = SEED
  )
}
