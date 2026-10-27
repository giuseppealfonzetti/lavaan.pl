theta_block <- function(MODEL, NTHR) {
  q <- MODEL$q
  rep(c("tau", "lambda", "corr", "logsd", "beta", "rescov"),
      c(NTHR, sum(is.na(MODEL$CONSTRMAT)), if (MODEL$CORRFLAG) q * (q - 1) / 2 else 0, sum(is.na(MODEL$CONSTRLOGSD)),
        sum(is.na(MODEL$REGMAT)), sum(is.na(MODEL$RESCOV[upper.tri(MODEL$RESCOV)]))))
}

latpl_dims <- function(DATA, MODEL) {
  stopifnot(MODEL$q < DATA$p)
  nb <- table(factor(theta_block(MODEL, sum(DATA$cat - 1)), c("tau", "lambda", "corr", "logsd", "beta", "rescov")))
  list(n = DATA$n, p = DATA$p, q = MODEL$q, cat = DATA$cat, pairs = DATA$p * (DATA$p - 1) / 2,
       nthr = nb[[1]], nload = nb[[2]], ncorr = nb[[3]], nvar = nb[[4]], nreg = nb[[5]], nres = nb[[6]], d = sum(nb))
}

corrflag_code <- function(MODEL) {
  if (!MODEL$CORRFLAG) 0L else if (identical(MODEL$LATCOV, "direct")) 2L else 1L
}


#' Names of the model parameters
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()].
#' @return The names of the parameters in the order used by [fit_latpl()]: thresholds `tau[i,k]`,
#'   free loadings `lambda[j,h]`, latent correlation parameters `eta[r,c]` (or covariances `psi[r,c]`),
#'   free latent log standard deviations `logsd[h]`, free regression coefficients `beta[i,j]` and free
#'   residual covariances `theta[k,l]`.
#' @export
theta_names <- function(DATA, MODEL) {
  stopifnot(inherits(DATA, "latpl_data"), inherits(MODEL, "latpl_model"))
  idx <- which(is.na(MODEL$CONSTRMAT), arr.ind = TRUE)
  q <- MODEL$q
  c(
    unlist(lapply(seq_len(DATA$p), function(i) sprintf("tau[%d,%d]", i, seq_len(DATA$cat[i] - 1)))),
    sprintf("lambda[%d,%d]", idx[, 1], idx[, 2]),
    if (MODEL$CORRFLAG && q > 1) unlist(lapply(2:q, function(h) {
      sprintf(if (corrflag_code(MODEL) == 2) "psi[%d,%d]" else "eta[%d,%d]", seq_len(h - 1), h)
    })),
    sprintf("logsd[%d]", which(is.na(MODEL$CONSTRLOGSD))),
    sprintf("beta[%d,%d]", row(MODEL$REGMAT)[is.na(MODEL$REGMAT)], col(MODEL$REGMAT)[is.na(MODEL$REGMAT)]),
    sprintf("theta[%d,%d]", row(MODEL$RESCOV)[is.na(MODEL$RESCOV) & upper.tri(MODEL$RESCOV)],
            col(MODEL$RESCOV)[is.na(MODEL$RESCOV) & upper.tri(MODEL$RESCOV)])
  )
}

lkj <- function(H, Q) {
  U <- diag(Q)
  z <- matrix(0, Q, Q)
  z[upper.tri(z)] <- tanh(H)
  for (s in seq_len(Q)[-1]) {
    P <- cumprod(c(1, sqrt(1 - z[seq_len(s - 1), s]^2)))
    U[seq_len(s), s] <- c(z[seq_len(s - 1), s], 1) * P
  }
  U
}

lkj_inv <- function(R) {
  U <- chol(R)
  z <- matrix(0, ncol(R), ncol(R))
  for (s in seq_len(ncol(R))[-1]) {
    for (r in seq_len(s - 1)) z[r, s] <- U[r, s] / prod(sqrt(1 - z[seq_len(r - 1), s]^2))
  }
  atanh(z[upper.tri(z)])
}

#' Convert between the parameter vector and the model matrices
#'
#' @param MODEL Model object from [set_latpl_model()].
#' @param THETA Parameter vector, in the order of [theta_names()].
#' @param PAR List with the thresholds `thresholds`, the loading matrix `loadings`, the latent
#'   (disturbance) covariance matrix `latent_correlations`, and, when the model has them, the
#'   regression matrix `regressions` and the residual covariance matrix `residual_covariances`.
#' @return `theta2par()` returns the list `PAR`; `par2theta()` returns the parameter vector.
#' @export
theta2par <- function(MODEL, THETA) {
  stopifnot(inherits(MODEL, "latpl_model"), is.numeric(THETA))
  b <- theta_block(MODEL, length(THETA) - length(theta_block(MODEL, 0)))
  q <- MODEL$q
  p <- MODEL$p
  L <- MODEL$CONSTRMAT
  L[is.na(L)] <- THETA[b == "lambda"]
  sd <- MODEL$CONSTRLOGSD
  sd[is.na(sd)] <- THETA[b == "logsd"]
  sd <- exp(sd)
  Psi <- diag(sd^2, q)
  if (corrflag_code(MODEL) == 2) {
    Psi[upper.tri(Psi)] <- THETA[b == "corr"]
    Psi[lower.tri(Psi)] <- t(Psi)[lower.tri(Psi)]
  } else if (any(b == "corr")) {
    Psi <- sd * crossprod(lkj(THETA[b == "corr"], q)) * rep(sd, each = q)
  }
  B <- matrix(0, q, q)
  if (length(MODEL$REGMAT)) {
    B <- MODEL$REGMAT
    B[is.na(B)] <- THETA[b == "beta"]
  }
  R <- matrix(0, p, p)
  if (length(MODEL$RESCOV)) {
    R <- MODEL$RESCOV
    R[is.na(R) & upper.tri(R)] <- THETA[b == "rescov"]
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
  }
  list(thresholds = THETA[b == "tau"], loadings = L, latent_correlations = Psi, regressions = B, residual_covariances = R)
}

#' @rdname theta2par
#' @export
par2theta <- function(MODEL, PAR) {
  stopifnot(inherits(MODEL, "latpl_model"))
  Psi <- PAR$latent_correlations
  c(PAR$thresholds,
    PAR$loadings[is.na(MODEL$CONSTRMAT)],
    if (corrflag_code(MODEL) == 2) Psi[upper.tri(Psi)] else if (MODEL$CORRFLAG && MODEL$q > 1) lkj_inv(stats::cov2cor(Psi)),
    log(sqrt(diag(Psi)))[is.na(MODEL$CONSTRLOGSD)],
    if (length(MODEL$REGMAT)) PAR$regressions[is.na(MODEL$REGMAT)],
    if (length(MODEL$RESCOV)) PAR$residual_covariances[is.na(MODEL$RESCOV) & upper.tri(MODEL$RESCOV)])
}

lincon <- function(MODEL, NTHR) {
  block <- theta_block(MODEL, NTHR)
  d <- length(block)
  if (is.null(MODEL$CONSTRAINTS) || nrow(MODEL$CONSTRAINTS$A) == 0) return(list(K = diag(d), c = rep(0, d), free_idx = seq_len(d)))
  A <- MODEL$CONSTRAINTS$A
  b <- MODEL$CONSTRAINTS$b
  if (ncol(A) != d) stop("CONSTRAINTS$A must have one column per parameter (", d, "), in the order of theta_names().")
  for (r in seq_len(nrow(A))) {
    blk <- unique(block[A[r, ] != 0])
    if (any(blk == "corr") && corrflag_code(MODEL) != 2) stop("Constraints on the latent correlations require LATCOV = \"direct\".")
    if (length(blk) > 1) stop("Each constraint must involve either thresholds, loadings, latent covariances, latent log standard deviations, regressions or residual covariances.")
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
