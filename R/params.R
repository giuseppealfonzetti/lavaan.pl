#' Parameter blocks
#'
#' @param MODEL Model object from [lavaanpl_model()] or [lavaanpl_groups()].
#' @param NTHR Number of thresholds in each group.
#' @return A character vector with the block of each parameter, in the order used by [fit_lavaanpl()]: `"tau"`,
#'   `"lambda"`, `"corr"`, `"logsd"`, `"beta"`, `"rescov"`, `"nu"`, `"alpha"`, `"scale"`, `"gamma"` or `"kappa"`, group
#'   after group.
#' @seealso [parameter_names()], [model_dims()]
#' @keywords internal
parameter_blocks <- function(MODEL, NTHR) {
  if (inherits(MODEL, "lavaanpl_groups")) return(unlist(lapply(MODEL$MODELS, parameter_blocks, NTHR)))
  q <- MODEL$q
  rep(c("tau", "lambda", "corr", "logsd", "beta", "rescov", "nu", "alpha", "scale", "gamma", "kappa"),
      c(NTHR, sum(is.na(MODEL$LAMBDA)), if (MODEL$LATCOV != "none") q * (q - 1) / 2 else 0, sum(is.na(diag(MODEL$PSI))),
        sum(is.na(MODEL$BETA)), sum(is.na(MODEL$THETA[upper.tri(MODEL$THETA)])), sum(is.na(MODEL$NU)),
        sum(is.na(MODEL$ALPHA)), sum(is.na(item_scales(MODEL))), sum(is.na(MODEL$GAMMA)), sum(is.na(MODEL$KAPPA))))
}

#' Scales of the items
#'
#' @param MODEL Model object from [lavaanpl_model()].
#' @return The residual standard deviation of each item under `PARAMETERIZATION = "theta"` and of the numeric items,
#'   the scale factor `DELTA` of the other items, `NA` when free.
#' @seealso [parameter_blocks()], [vector_to_matrices()]
#' @keywords internal
item_scales <- function(MODEL) {
  ifelse(MODEL$PARAMETERIZATION == "theta" | MODEL$NUMERIC, sqrt(diag(MODEL$THETA)), MODEL$DELTA)
}

#' Dimensions of a fit
#'
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_groups()].
#' @return A list with the number of latent variables `q`, the size of each parameter block summed over groups
#'   (`nthr`, `nload`, `ncorr`, `nvar`, `nreg`, `nres`, `nint`, `nmean`, `nscale`, `ngam`, `nkap`) and the number of
#'   parameters `d`.
#' @seealso [parameter_blocks()], [fit_lavaanpl()]
#' @keywords internal
model_dims <- function(DATA, MODEL) {
  models <- MODEL$MODELS
  stopifnot(all(vapply(models, function(x) x$q, numeric(1)) < DATA$p))
  for (m in models) {
    if (ncol(m$GAMMA) > 0 && ncol(m$GAMMA) != ncol(DATA$X)) stop("GAMMA and KAPPA must have one column per covariate of DATA.")
  }
  nb <- table(factor(parameter_blocks(MODEL, sum(DATA$cat - 1)), c("tau", "lambda", "corr", "logsd", "beta", "rescov", "nu", "alpha", "scale", "gamma", "kappa")))
  list(q = models[[1]]$q, nthr = nb[[1]], nload = nb[[2]], ncorr = nb[[3]], nvar = nb[[4]], nreg = nb[[5]], nres = nb[[6]],
       nint = nb[[7]], nmean = nb[[8]], nscale = nb[[9]], ngam = nb[[10]], nkap = nb[[11]], d = sum(nb))
}

#' Names of the model parameters
#'
#' @param DATA Data object from [lavaanpl_data()].
#' @param MODEL Model object from [lavaanpl_model()].
#' @return The names of the parameters in the order used by [fit_lavaanpl()]: thresholds `tau[i,k]`, then the free
#'   loadings `lambda[j,h]`, latent correlation parameters `eta[r,c]` (covariances `psi[r,c]` with
#'   `LATCOV = "direct"`), latent log standard deviations `logsd[h]`, regressions `beta[i,j]`, residual covariances
#'   `theta[k,l]`, intercepts `nu[j]`, latent means `alpha[h]`, log scale factors `logdelta[j]` or log residual
#'   standard deviations `logrsd[j]`, effects of the covariates `gamma[h,c]` and `kappa[j,c]`. With
#'   [lavaanpl_groups()], the names of group 2 onwards end in `.g2`, `.g3`, and so on.
#' @seealso [parameter_blocks()]
#' @keywords internal
parameter_names <- function(DATA, MODEL) {
  stopifnot(inherits(DATA, "lavaanpl_data"), inherits(MODEL, c("lavaanpl_model", "lavaanpl_groups")))
  if (inherits(MODEL, "lavaanpl_groups")) {
    return(unlist(lapply(seq_len(MODEL$ngroups), function(g) {
      paste0(parameter_names(DATA, MODEL$MODELS[[g]]), if (g > 1) paste0(".g", g))
    })))
  }
  idx <- which(is.na(MODEL$LAMBDA), arr.ind = TRUE)
  q <- MODEL$q
  scales <- item_scales(MODEL)
  c(
    unlist(lapply(seq_len(DATA$p), function(i) sprintf("tau[%d,%d]", i, seq_len(DATA$cat[i] - 1)))),
    sprintf("lambda[%d,%d]", idx[, 1], idx[, 2]),
    if (MODEL$LATCOV != "none" && q > 1) unlist(lapply(2:q, function(h) {
      sprintf(if (MODEL$LATCOV == "direct") "psi[%d,%d]" else "eta[%d,%d]", seq_len(h - 1), h)
    })),
    sprintf("logsd[%d]", which(is.na(diag(MODEL$PSI)))),
    sprintf("beta[%d,%d]", row(MODEL$BETA)[is.na(MODEL$BETA)], col(MODEL$BETA)[is.na(MODEL$BETA)]),
    sprintf("theta[%d,%d]", row(MODEL$THETA)[is.na(MODEL$THETA) & upper.tri(MODEL$THETA)],
            col(MODEL$THETA)[is.na(MODEL$THETA) & upper.tri(MODEL$THETA)]),
    sprintf("nu[%d]", which(is.na(MODEL$NU))),
    sprintf("alpha[%d]", which(is.na(MODEL$ALPHA))),
    sprintf(ifelse(MODEL$PARAMETERIZATION == "theta" | MODEL$NUMERIC, "logrsd[%d]", "logdelta[%d]")[is.na(scales)], which(is.na(scales))),
    sprintf("gamma[%d,%d]", row(MODEL$GAMMA)[is.na(MODEL$GAMMA)], col(MODEL$GAMMA)[is.na(MODEL$GAMMA)]),
    sprintf("kappa[%d,%d]", row(MODEL$KAPPA)[is.na(MODEL$KAPPA)], col(MODEL$KAPPA)[is.na(MODEL$KAPPA)])
  )
}

#' Cholesky factor from LKJ parameters
#'
#' @param H Unconstrained latent correlation parameters, one per pair of latent variables.
#' @param Q Number of latent variables.
#' @return The upper triangular Cholesky factor of the latent correlation matrix.
#' @seealso [lkj_inv()], [vector_to_matrices()]
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
#' @seealso [lkj()], [matrices_to_vector()]
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
#' @param MODEL Model object from [lavaanpl_model()] or [lavaanpl_groups()].
#' @param PARAMETER_VECTOR Parameter vector, in the order of [parameter_names()].
#' @param MATRICES List of model matrices named as in lavaan: `tau`, `lambda`, `psi`, `beta`, `theta`, `nu`, `alpha`,
#'   `delta` and, with covariates, `gamma` and `kappa`.
#' @return `vector_to_matrices()` returns the list `MATRICES`, or with [lavaanpl_groups()] one per group;
#'   `matrices_to_vector()` returns the parameter vector.
#' @seealso [parameter_names()]
#' @keywords internal
vector_to_matrices <- function(MODEL, PARAMETER_VECTOR) {
  stopifnot(inherits(MODEL, c("lavaanpl_model", "lavaanpl_groups")), is.numeric(PARAMETER_VECTOR))
  if (inherits(MODEL, "lavaanpl_groups")) {
    d <- vapply(MODEL$MODELS, function(x) length(parameter_blocks(x, 0)), numeric(1))
    d <- d + (length(PARAMETER_VECTOR) - sum(d)) / length(d)
    return(Map(vector_to_matrices, MODEL$MODELS, split(PARAMETER_VECTOR, rep(seq_along(d), d))))
  }
  b <- parameter_blocks(MODEL, length(PARAMETER_VECTOR) - length(parameter_blocks(MODEL, 0)))
  q <- MODEL$q
  p <- MODEL$p
  L <- MODEL$LAMBDA
  L[is.na(L)] <- PARAMETER_VECTOR[b == "lambda"]
  sd <- sqrt(diag(MODEL$PSI))
  sd[is.na(sd)] <- exp(PARAMETER_VECTOR[b == "logsd"])
  Psi <- diag(sd^2, q)
  if (MODEL$LATCOV == "direct") {
    Psi[upper.tri(Psi)] <- PARAMETER_VECTOR[b == "corr"]
    Psi[lower.tri(Psi)] <- t(Psi)[lower.tri(Psi)]
  } else if (any(b == "corr")) {
    Psi <- sd * crossprod(lkj(PARAMETER_VECTOR[b == "corr"], q)) * rep(sd, each = q)
  }
  B <- MODEL$BETA
  B[is.na(B)] <- PARAMETER_VECTOR[b == "beta"]
  R <- MODEL$THETA
  R[is.na(R) & upper.tri(R)] <- PARAMETER_VECTOR[b == "rescov"]
  R[lower.tri(R)] <- t(R)[lower.tri(R)]
  nu <- MODEL$NU
  nu[is.na(nu)] <- PARAMETER_VECTOR[b == "nu"]
  al <- MODEL$ALPHA
  al[is.na(al)] <- PARAMETER_VECTOR[b == "alpha"]
  ls <- log(item_scales(MODEL))
  ls[is.na(ls)] <- PARAMETER_VECTOR[b == "scale"]
  thetapar <- MODEL$PARAMETERIZATION == "theta"
  rsd <- thetapar | MODEL$NUMERIC
  diag(R)[rsd] <- exp(2 * ls[rsd])
  scales <- exp(ls)
  scales[MODEL$NUMERIC] <- 1
  Ga <- MODEL$GAMMA
  Ga[is.na(Ga)] <- PARAMETER_VECTOR[b == "gamma"]
  Ka <- MODEL$KAPPA
  Ka[is.na(Ka)] <- PARAMETER_VECTOR[b == "kappa"]
  c(list(tau = PARAMETER_VECTOR[b == "tau"], lambda = L, psi = Psi, beta = B, theta = R, nu = nu, alpha = al),
    if (!thetapar) list(delta = scales), if (ncol(Ga) > 0) list(gamma = Ga, kappa = Ka))
}

#' @rdname vector_to_matrices
matrices_to_vector <- function(MODEL, MATRICES) {
  stopifnot(inherits(MODEL, c("lavaanpl_model", "lavaanpl_groups")))
  if (inherits(MODEL, "lavaanpl_groups")) return(unlist(Map(matrices_to_vector, MODEL$MODELS, MATRICES), use.names = FALSE))
  Psi <- MATRICES$psi
  scales <- item_scales(MODEL)
  c(MATRICES$tau,
    MATRICES$lambda[is.na(MODEL$LAMBDA)],
    if (MODEL$LATCOV == "direct") Psi[upper.tri(Psi)] else if (MODEL$LATCOV == "lkj" && MODEL$q > 1) lkj_inv(stats::cov2cor(Psi)),
    log(sqrt(diag(Psi)))[is.na(diag(MODEL$PSI))],
    MATRICES$beta[is.na(MODEL$BETA)],
    MATRICES$theta[is.na(MODEL$THETA) & upper.tri(MODEL$THETA)],
    MATRICES$nu[is.na(MODEL$NU)],
    MATRICES$alpha[is.na(MODEL$ALPHA)],
    if (anyNA(scales)) {
      rsd <- MODEL$PARAMETERIZATION == "theta" | MODEL$NUMERIC
      log(ifelse(rsd, sqrt(diag(MATRICES$theta)), if (is.null(MATRICES$delta)) 1 else MATRICES$delta))[is.na(scales)]
    },
    MATRICES$gamma[is.na(MODEL$GAMMA)],
    MATRICES$kappa[is.na(MODEL$KAPPA)])
}

#' Linear constraints on the parameters
#'
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param NTHR Number of thresholds in each group.
#' @return A list with the matrix `K` and the vector `c` that map the free parameters to the full parameter
#'   vector, and the positions `free_idx` of the free parameters in it.
#' @seealso [fit_lavaanpl()], [parameter_names()]
#' @keywords internal
lincon <- function(MODEL, NTHR) {
  models <- MODEL$MODELS
  block <- parameter_blocks(MODEL, NTHR)
  d <- length(block)
  dg <- vapply(models, function(x) length(parameter_blocks(x, NTHR)), numeric(1))
  direct <- rep(vapply(models, function(x) x$LATCOV == "direct", logical(1)), dg)
  cons <- c(lapply(models, function(x) x$CONSTRAINTS), list(MODEL$CONSTRAINTS))
  cols <- c(split(seq_len(d), rep(seq_along(dg), dg)), list(seq_len(d)))
  A <- matrix(0, 0, d)
  b <- numeric(0)
  for (g in seq_along(cons)) {
    if (is.null(cons[[g]]) || nrow(cons[[g]]$A) == 0) next
    if (ncol(cons[[g]]$A) != length(cols[[g]])) {
      stop("CONSTRAINTS$A must have one column per parameter (", length(cols[[g]]), "), in the order of parameter_names().")
    }
    Ag <- matrix(0, nrow(cons[[g]]$A), d)
    Ag[, cols[[g]]] <- cons[[g]]$A
    A <- rbind(A, Ag)
    b <- c(b, cons[[g]]$b)
  }
  fixed <- do.call(rbind, Map(function(m, pos) {
    tau <- if (is.null(m$TAU)) rep(NA, NTHR) else m$TAU
    if (length(tau) != NTHR) stop("TAU must have one element per threshold (", NTHR, ").")
    value <- c(tau, if (m$LATCOV == "direct") m$PSI[upper.tri(m$PSI)])
    at <- pos[block[pos] == "tau" | (m$LATCOV == "direct" & block[pos] == "corr")]
    cbind(at, value)[!is.na(value), , drop = FALSE]
  }, models, cols[seq_along(models)]))
  Af <- matrix(0, nrow(fixed), d)
  Af[cbind(seq_len(nrow(fixed)), fixed[, 1])] <- 1
  A <- rbind(A, Af)
  b <- c(b, fixed[, 2])
  if (nrow(A) == 0) return(list(K = diag(d), c = rep(0, d), free_idx = seq_len(d)))
  for (r in seq_len(nrow(A))) {
    blk <- unique(block[A[r, ] != 0])
    if (any(blk == "corr" & !direct[A[r, ] != 0])) stop("Constraints on the latent correlations require LATCOV = \"direct\".")
    if (length(blk) > 1) stop("Each constraint must involve either thresholds, loadings, latent covariances, latent log standard deviations, regressions, residual covariances, intercepts, latent means, log scales, or effects of the covariates.")
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

#' Interpretable parameters
#'
#' @param MODEL Model object from [lavaanpl_groups()].
#' @param PARAMETER_VECTOR Parameter vector, in the order of [parameter_names()].
#' @return The parameter vector with latent correlations instead of their parameters, latent standard deviations, and
#'   scale factors or residual variances instead of their logarithms.
#' @keywords internal
transform_parameters <- function(MODEL, PARAMETER_VECTOR) {
  unlist(Map(function(m, mat) {
    R <- stats::cov2cor(mat$psi)
    rsd <- m$PARAMETERIZATION == "theta" | m$NUMERIC
    c(mat$tau,
      mat$lambda[is.na(m$LAMBDA)],
      if (m$LATCOV != "none" && m$q > 1) R[upper.tri(R)],
      sqrt(diag(mat$psi))[is.na(diag(m$PSI))],
      mat$beta[is.na(m$BETA)],
      mat$theta[is.na(m$THETA) & upper.tri(m$THETA)],
      mat$nu[is.na(m$NU)],
      mat$alpha[is.na(m$ALPHA)],
      ifelse(rsd, diag(mat$theta), if (is.null(mat$delta)) 1 else mat$delta)[is.na(item_scales(m))],
      mat$gamma[is.na(m$GAMMA)],
      mat$kappa[is.na(m$KAPPA)])
  }, MODEL$MODELS, vector_to_matrices(MODEL, PARAMETER_VECTOR)), use.names = FALSE)
}
