#' Pairwise likelihood ratio test
#' @param OBJ Object of class latpl_fit
#' @param VARS Output of [compute_var()] for `OBJ`
#' @param H1 Object of class latpl_fit in which `OBJ` is nested, or `NULL` to test `OBJ` against the
#'   saturated model
#' @return A list with the test of the model against `H1` or the saturated model `model` and, without
#'   `H1`, of the independence model against the saturated model `baseline`. Each holds the
#'   statistic `stat`, its degrees of freedom `df`, and the mean and variance adjusted `stat.scaled`
#'   (`stat` divided by `scaling.factor`), `df.scaled` and `pvalue.scaled`. Also, without `H1`, the
#'   saturated and independence fits `saturated` and `independence`, and the timing `RTime`.
#' @export
compute_plrt <- function(OBJ, VARS = compute_var(OBJ), H1 = NULL) UseMethod('compute_plrt')
#' Pairwise likelihood ratio test
#' @param OBJ Object of class latpl_fit
#' @param VARS Output of [compute_var()] for `OBJ`
#' @param H1 Object of class latpl_fit in which `OBJ` is nested, or `NULL` to test `OBJ` against the
#'   saturated model
#' @return A list with the test of the model against `H1` or the saturated model `model` and, without
#'   `H1`, of the independence model against the saturated model `baseline`. Each holds the
#'   statistic `stat`, its degrees of freedom `df`, and the mean and variance adjusted `stat.scaled`
#'   (`stat` divided by `scaling.factor`), `df.scaled` and `pvalue.scaled`. Also, without `H1`, the
#'   saturated and independence fits `saturated` and `independence`, and the timing `RTime`.
#' @export
compute_plrt.latpl_fit <- function(OBJ, VARS = compute_var(OBJ), H1 = NULL){
  start_time <- Sys.time()
  RcppParallel::setThreadOptions(numThreads = OBJ$cores)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  out <- if (is.null(H1)) plrt_fit(OBJ, VARS) else list(model = plrt_nested(OBJ, H1))
  out$RTime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  out
}

plrt_adjust <- function(STAT, DF, E, V){
  c1 <- E / (V / 2)
  list(stat = STAT, df = DF, stat.scaled = c1 * STAT, df.scaled = E^2 / (V / 2), scaling.factor = 1 / c1,
       pvalue.scaled = if(isTRUE(c1 > 0)) stats::pchisq(c1 * STAT, E^2 / (V / 2), lower.tail = FALSE) else NA)
}

# Nested models: H1 evaluated at the H0 estimates, with the exact restriction given
# by the Jacobian of the embedding of OBJ into H1. The embedding goes through the model
# matrices; without regressions in H1, through the implied latent covariance, rescaled
# to the fixed latent variances of H1
plrt_nested <- function(OBJ, H1){
  stopifnot(inherits(H1, "latpl_fit"))
  if(!identical(OBJ$data$freq, H1$data$freq) || !identical(OBJ$data$ufreq, H1$data$ufreq) ||
     !identical(OBJ$data$weights, H1$data$weights) || !identical(OBJ$data$cluster, H1$data$cluster) ||
     H1$dims$npar <= OBJ$dims$npar) {
    stop('OBJ is not nested in H1.')
  }
  pl <- function(FIT, THETA) cpp_pl(THETA, FIT$engine$model, FIT$engine$data)$pl
  lc0 <- OBJ$lincon
  lc1 <- H1$lincon
  emb <- function(x){
    par <- theta2par(OBJ$model, as.vector(lc0$K %*% x + lc0$c))
    if(!length(H1$model$REGMAT)){
      A <- solve(diag(OBJ$dims$q) - par$regressions)
      S <- A %*% par$latent_correlations %*% t(A)
      sc <- ifelse(is.finite(H1$model$CONSTRLOGSD), sqrt(diag(S)) / exp(H1$model$CONSTRLOGSD), 1)
      par$loadings <- par$loadings %*% diag(sc, OBJ$dims$q)
      par$latent_correlations <- S / tcrossprod(sc)
    }
    par2theta(H1$model, par)[lc1$free_idx]
  }
  theta0f <- OBJ$theta[lc0$free_idx]
  theta1 <- as.vector(lc1$K %*% emb(theta0f) + lc1$c)
  pl0 <- pl(OBJ, OBJ$theta)
  if(abs(pl(H1, theta1) - pl0) > 1e-8 * abs(pl0)) stop('OBJ is not nested in H1.')
  N <- numDeriv::jacobian(emb, theta0f)
  HJ <- cpp_HJ(theta1, H1$engine$model, H1$engine$data)
  Hf <- crossprod(lc1$K, HJ$H %*% lc1$K)
  M <- crossprod(lc1$K, HJ$J %*% lc1$K) %*% (solve(Hf) - N %*% solve(crossprod(N, Hf %*% N), t(N)))
  plrt_adjust(STAT = 2 * (pl(H1, H1$theta) - pl0), DF = H1$dims$npar - OBJ$dims$npar, E = sum(diag(M)), V = 2 * sum(M * t(M)))
}

# Model and independence model against the saturated model
plrt_fit <- function(OBJ, VARS){
  dims <- OBJ$dims
  n <- dims$n
  nthr <- dims$nthr
  eng <- OBJ$engine$data
  sat <- function(tau, rho, order, scores = FALSE) cpp_saturated(tau, rho, eng, order, scores)

  # Independence model: thresholds from the univariate margins, with unit weights
  # given by the number of pairs and univariate terms each observed item enters
  Y <- OBJ$data$Y
  v <- (rowSums(!is.na(Y)) - 1 + OBJ$data$uweights) * (if(length(OBJ$data$weights)) OBJ$data$weights else 1)
  tau_ind <- unlist(lapply(seq_len(dims$p), function(j){
    pj <- vapply(seq_len(dims$cat[j]) - 1, function(a) sum(v[Y[, j] %in% a]), numeric(1))
    stats::qnorm(cumsum(pj)[-dims$cat[j]] / sum(pj))
  }))

  # Saturated model: Newton steps through the Schur complement of the diagonal
  # correlations block
  tau <- tau_ind
  rho <- rep(0, dims$pairs)
  o <- sat(tau, rho, 2)
  for(iter in 1:100){
    if(max(abs(o$gradient)) / n < 1e-8) break
    gt <- o$gradient[1:nthr]
    gr <- o$gradient[-(1:nthr)]
    dt <- as.vector(solve(o$hschur, gt - o$htr %*% (gr / o$hrr)))
    dr <- as.vector(gr - crossprod(o$htr, dt)) / o$hrr
    step <- 1
    repeat{
      if(all(abs(rho + step * dr) < 1)){
        o_new <- sat(tau + step * dt, rho + step * dr, 2)
        if(isTRUE(o_new$ll >= o$ll - 1e-12 * abs(o$ll))) break
      }
      step <- step / 2
      if(step < 1e-10) stop('Saturated model fit failed.')
    }
    tau <- tau + step * dt
    rho <- rho + step * dr
    o <- o_new
  }
  if(max(abs(o$gradient)) / n >= 1e-8) warning('Saturated model fit did not converge.')

  # Traces of the saturated model at (TAU, RHO), from per unit information blocks
  # and the weighted scores of the distinct response patterns
  sat_terms <- function(TAU, RHO, DELTA = NULL, FULL = FALSE){
    s <- sat(TAU, RHO, 2, TRUE)
    htt <- s$htt / n
    htr <- s$htr / n
    d <- s$hrr / n
    w <- function(S) if(length(eng$WM)) S %*% eng$WM else S * rep(eng$WV, each = nrow(S))
    tr <- cpp_saturated_traces(htt, htr, d, s$hschur / n, w(s$st), w(s$sr), n, FULL)
    out <- list(ll = s$ll, e = tr[1], v = 2 * tr[2])
    if(!is.null(DELTA) && FULL){
      ht <- crossprod(DELTA$tau, htr %*% DELTA$rho)
      out$dad <- crossprod(DELTA$tau, htt %*% DELTA$tau) + ht + t(ht) + crossprod(DELTA$rho, DELTA$rho * d)
    }else if(!is.null(DELTA)){
      hd <- htr %*% DELTA
      out$dad <- crossprod(DELTA, DELTA * d) - crossprod(hd, solve(htt, hd))
    }
    out
  }

  # Model against saturated, on the free parameters. With constrained thresholds,
  # the traces run over the full parameter vectors (lavaan's lav_pml_plrt2)
  lc <- OBJ$lincon
  r <- cpp_rho(OBJ$theta, OBJ$engine$model, eng)
  jac <- if (ncol(lc$K) < nrow(lc$K)) r$jacobian %*% lc$K else r$jacobian
  invH <- solve(VARS$H)
  if(all(seq_len(nthr) %in% lc$free_idx)){
    psi <- which(lc$free_idx > nthr)
    G <- (invH %*% VARS$J %*% invH)[psi, psi, drop = FALSE]
    P0 <- solve(invH[psi, psi, drop = FALSE], G)
    tz <- sat_terms(OBJ$theta[1:nthr], r$rho, jac[, psi, drop = FALSE])
    cv <- sum(diag(tz$dad %*% G %*% P0))
  }else{
    P0 <- invH %*% VARS$J
    tz <- sat_terms(OBJ$theta[1:nthr], r$rho, list(tau = lc$K[seq_len(nthr), , drop = FALSE], rho = jac), FULL = TRUE)
    cv <- sum(diag(tz$dad %*% P0 %*% P0 %*% invH))
  }
  model <- plrt_adjust(STAT = 2 * (o$ll - tz$ll),
                       DF   = nthr + dims$pairs - dims$npar,
                       E    = tz$e - sum(diag(P0)),
                       V    = tz$v + 2 * sum(P0 * t(P0)) - 4 * cv)

  # Independence model against saturated
  tb <- sat_terms(tau_ind, rep(0, dims$pairs))
  baseline <- plrt_adjust(STAT = 2 * (o$ll - tb$ll), DF = dims$pairs, E = tb$e, V = tb$v)

  list(
    model = model,
    baseline = baseline,
    saturated = list(tau = tau, rho = rho, logl = o$ll, iterations = iter),
    independence = list(tau = tau_ind, logl = tb$ll)
  )
}
