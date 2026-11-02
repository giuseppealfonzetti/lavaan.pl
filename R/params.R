#' Parameter blocks
#'
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param NTHR Number of thresholds in each group.
#' @return A character vector with the block of each parameter, in the order used by [fit_latpl()]: `"tau"`,
#'   `"lambda"`, `"corr"`, `"logsd"`, `"beta"`, `"rescov"`, `"nu"`, `"alpha"` or `"scale"`, group after group.
#' @seealso [theta_names()], [latpl_dims()]
#' @keywords internal
theta_block <- function(MODEL, NTHR) {
  if (inherits(MODEL, "latpl_groups")) return(unlist(lapply(MODEL$MODELS, theta_block, NTHR)))
  q <- MODEL$q
  rep(c("tau", "lambda", "corr", "logsd", "beta", "rescov", "nu", "alpha", "scale"),
      c(NTHR, sum(is.na(MODEL$CONSTRMAT)), if (MODEL$CORRFLAG) q * (q - 1) / 2 else 0, sum(is.na(MODEL$CONSTRLOGSD)),
        sum(is.na(MODEL$REGMAT)), sum(is.na(MODEL$RESCOV[upper.tri(MODEL$RESCOV)])), sum(is.na(MODEL$CONSTRINT)),
        sum(is.na(MODEL$CONSTRMEAN)), sum(is.na(MODEL$CONSTRLOGSCALE))))
}

#' Dimensions of a fit
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @return A list with the sample size `n`, the numbers of items `p` and latent variables `q`, the numbers of
#'   categories `cat`, the number of groups `ngroups`, the number of pairs `pairs` over all groups, the size of each
#'   parameter block summed over groups (`nthr`, `nload`, `ncorr`, `nvar`, `nreg`, `nres`, `nint`, `nmean`,
#'   `nscale`) and the number of parameters `d`.
#' @seealso [theta_block()], [fit_latpl()]
#' @keywords internal
latpl_dims <- function(DATA, MODEL) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
  stopifnot(all(vapply(models, function(x) x$q, numeric(1)) < DATA$p))
  nb <- table(factor(theta_block(MODEL, sum(DATA$cat - 1)), c("tau", "lambda", "corr", "logsd", "beta", "rescov", "nu", "alpha", "scale")))
  list(n = DATA$n, p = DATA$p, q = models[[1]]$q, cat = DATA$cat, ngroups = DATA$ngroups,
       pairs = DATA$ngroups * DATA$p * (DATA$p - 1) / 2,
       nthr = nb[[1]], nload = nb[[2]], ncorr = nb[[3]], nvar = nb[[4]], nreg = nb[[5]], nres = nb[[6]],
       nint = nb[[7]], nmean = nb[[8]], nscale = nb[[9]], d = sum(nb))
}

#' Code of the latent covariance parameterisation
#'
#' @param MODEL Model object from [set_latpl_model()].
#' @return `0` for orthogonal latent variables, `1` for LKJ correlations, `2` for direct covariances.
#' @seealso [engine()], [lincon()]
#' @keywords internal
corrflag_code <- function(MODEL) {
  if (!MODEL$CORRFLAG) 0L else if (identical(MODEL$LATCOV, "direct")) 2L else 1L
}


#' Names of the model parameters
#'
#' @param DATA Data object from [set_latpl_data()].
#' @param MODEL Model object from [set_latpl_model()].
#' @return The names of the parameters in the order used by [fit_latpl()]: thresholds `tau[i,k]`,
#'   free loadings `lambda[j,h]`, latent correlation parameters `eta[r,c]` (or covariances `psi[r,c]`),
#'   free latent log standard deviations `logsd[h]`, free regression coefficients `beta[i,j]`, free
#'   residual covariances `theta[k,l]`, free intercepts `nu[j]`, free latent means `alpha[h]`, and free
#'   log scale factors `logdelta[j]` or, with `PARAMETERIZATION = "theta"`, free log residual standard
#'   deviations `logrsd[j]`. With [set_latpl_groups()], the names of group 2 onwards end in `.g2`,
#'   `.g3`, and so on.
#' @export
theta_names <- function(DATA, MODEL) {
  stopifnot(inherits(DATA, "latpl_data"), inherits(MODEL, c("latpl_model", "latpl_groups")))
  if (inherits(MODEL, "latpl_groups")) {
    return(unlist(lapply(seq_len(MODEL$ngroups), function(g) {
      paste0(theta_names(DATA, MODEL$MODELS[[g]]), if (g > 1) paste0(".g", g))
    })))
  }
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
            col(MODEL$RESCOV)[is.na(MODEL$RESCOV) & upper.tri(MODEL$RESCOV)]),
    sprintf("nu[%d]", which(is.na(MODEL$CONSTRINT))),
    sprintf("alpha[%d]", which(is.na(MODEL$CONSTRMEAN))),
    sprintf(if (MODEL$PARAMETERIZATION == "theta") "logrsd[%d]" else "logdelta[%d]", which(is.na(MODEL$CONSTRLOGSCALE)))
  )
}

#' Cholesky factor from LKJ parameters
#'
#' @param H Unconstrained latent correlation parameters, one per pair of latent variables.
#' @param Q Number of latent variables.
#' @return The upper triangular Cholesky factor of the latent correlation matrix.
#' @seealso [lkj_inv()], [theta2par()]
#' @keywords internal
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

#' LKJ parameters from a correlation matrix
#'
#' @param R Latent correlation matrix.
#' @return The unconstrained latent correlation parameters, one per pair of latent variables.
#' @seealso [lkj()], [par2theta()]
#' @keywords internal
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
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param THETA Parameter vector, in the order of [theta_names()].
#' @param PAR List with the thresholds `thresholds`, the loading matrix `loadings`, the latent
#'   (disturbance) covariance matrix `latent_correlations`, the regression matrix `regressions`, the
#'   residual covariance matrix `residual_covariances` (with the residual variances on the diagonal
#'   when `PARAMETERIZATION = "theta"`), the intercepts `intercepts`, the latent means `means` and,
#'   with `PARAMETERIZATION = "delta"`, the scale factors `scales`.
#' @return `theta2par()` returns the list `PAR`, or with [set_latpl_groups()] a list of them, one
#'   per group; `par2theta()` returns the parameter vector.
#' @export
theta2par <- function(MODEL, THETA) {
  stopifnot(inherits(MODEL, c("latpl_model", "latpl_groups")), is.numeric(THETA))
  if (inherits(MODEL, "latpl_groups")) {
    d <- vapply(MODEL$MODELS, function(x) length(theta_block(x, 0)), numeric(1))
    d <- d + (length(THETA) - sum(d)) / length(d)
    return(Map(theta2par, MODEL$MODELS, split(THETA, rep(seq_along(d), d))))
  }
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
  nu <- MODEL$CONSTRINT
  nu[is.na(nu)] <- THETA[b == "nu"]
  al <- MODEL$CONSTRMEAN
  al[is.na(al)] <- THETA[b == "alpha"]
  ls <- MODEL$CONSTRLOGSCALE
  ls[is.na(ls)] <- THETA[b == "scale"]
  theta <- MODEL$PARAMETERIZATION == "theta"
  if (theta) diag(R) <- exp(2 * ls)
  c(list(thresholds = THETA[b == "tau"], loadings = L, latent_correlations = Psi, regressions = B, residual_covariances = R,
         intercepts = nu, means = al), if (!theta) list(scales = exp(ls)))
}

#' @rdname theta2par
#' @export
par2theta <- function(MODEL, PAR) {
  stopifnot(inherits(MODEL, c("latpl_model", "latpl_groups")))
  if (inherits(MODEL, "latpl_groups")) return(unlist(Map(par2theta, MODEL$MODELS, PAR), use.names = FALSE))
  Psi <- PAR$latent_correlations
  c(PAR$thresholds,
    PAR$loadings[is.na(MODEL$CONSTRMAT)],
    if (corrflag_code(MODEL) == 2) Psi[upper.tri(Psi)] else if (MODEL$CORRFLAG && MODEL$q > 1) lkj_inv(stats::cov2cor(Psi)),
    log(sqrt(diag(Psi)))[is.na(MODEL$CONSTRLOGSD)],
    if (length(MODEL$REGMAT)) PAR$regressions[is.na(MODEL$REGMAT)],
    if (length(MODEL$RESCOV)) PAR$residual_covariances[is.na(MODEL$RESCOV) & upper.tri(MODEL$RESCOV)],
    PAR$intercepts[is.na(MODEL$CONSTRINT)],
    PAR$means[is.na(MODEL$CONSTRMEAN)],
    if (anyNA(MODEL$CONSTRLOGSCALE)) {
      log(if (MODEL$PARAMETERIZATION == "theta") sqrt(diag(PAR$residual_covariances)) else PAR$scales)[is.na(MODEL$CONSTRLOGSCALE)]
    })
}

#' Linear constraints on the parameters
#'
#' @param MODEL Model object from [set_latpl_model()] or [set_latpl_groups()].
#' @param NTHR Number of thresholds in each group.
#' @return A list with the matrix `K` and the vector `c` that map the free parameters to the full parameter
#'   vector, and the positions `free_idx` of the free parameters in it.
#' @seealso [fit_latpl()], [theta_names()]
#' @keywords internal
lincon <- function(MODEL, NTHR) {
  models <- if (inherits(MODEL, "latpl_groups")) MODEL$MODELS else list(MODEL)
  block <- theta_block(MODEL, NTHR)
  d <- length(block)
  dg <- vapply(models, function(x) length(theta_block(x, NTHR)), numeric(1))
  direct <- rep(vapply(models, function(x) corrflag_code(x) == 2, logical(1)), dg)
  cons <- c(lapply(models, function(x) x$CONSTRAINTS), list(if (length(models) > 1) MODEL$CONSTRAINTS))
  cols <- c(split(seq_len(d), rep(seq_along(dg), dg)), list(seq_len(d)))
  A <- matrix(0, 0, d)
  b <- numeric(0)
  for (g in seq_along(cons)) {
    if (is.null(cons[[g]]) || nrow(cons[[g]]$A) == 0) next
    if (ncol(cons[[g]]$A) != length(cols[[g]])) {
      stop("CONSTRAINTS$A must have one column per parameter (", length(cols[[g]]), "), in the order of theta_names().")
    }
    Ag <- matrix(0, nrow(cons[[g]]$A), d)
    Ag[, cols[[g]]] <- cons[[g]]$A
    A <- rbind(A, Ag)
    b <- c(b, cons[[g]]$b)
  }
  if (nrow(A) == 0) return(list(K = diag(d), c = rep(0, d), free_idx = seq_len(d)))
  for (r in seq_len(nrow(A))) {
    blk <- unique(block[A[r, ] != 0])
    if (any(blk == "corr" & !direct[A[r, ] != 0])) stop("Constraints on the latent correlations require LATCOV = \"direct\".")
    if (length(blk) > 1) stop("Each constraint must involve either thresholds, loadings, latent covariances, latent log standard deviations, regressions, residual covariances, intercepts, latent means or log scales.")
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
