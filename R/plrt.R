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
compute_plrt <- function(OBJ, VARS = compute_var(OBJ), H1 = NULL) {
  UseMethod('compute_plrt')
}
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
compute_plrt.latpl_fit <- function(OBJ, VARS = compute_var(OBJ), H1 = NULL) {
  start_time <- Sys.time()
  RcppParallel::setThreadOptions(numThreads = OBJ$cores)
  on.exit(RcppParallel::setThreadOptions(numThreads = 1))
  out <- if (is.null(H1)) {
    plrt_fit(OBJ, VARS)
  } else {
    list(model = plrt_nested(OBJ, H1))
  }
  out$RTime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  out
}

#' Mean and variance adjusted test
#'
#' @param STAT Test statistic, one per group.
#' @param DF Degrees of freedom.
#' @param E Mean of the statistic.
#' @param V Variance of the statistic.
#' @return A list with `stat`, `stat.group`, `df`, `stat.scaled`, `df.scaled`, `scaling.factor` and
#'   `pvalue.scaled`.
#' @seealso [compute_plrt()]
#' @keywords internal
plrt_adjust <- function(STAT, DF, E, V) {
  c1 <- E / (V / 2)
  list(
    stat = sum(STAT),
    stat.group = STAT,
    df = DF,
    stat.scaled = c1 * sum(STAT),
    df.scaled = E^2 / (V / 2),
    scaling.factor = 1 / c1,
    pvalue.scaled = if (isTRUE(c1 > 0)) {
      stats::pchisq(c1 * sum(STAT), E^2 / (V / 2), lower.tail = FALSE)
    } else {
      NA
    }
  )
}


#' Test of nested models
#'
#' @param OBJ Object of class latpl_fit, nested in `H1`.
#' @param H1 Object of class latpl_fit.
#' @return The output of [plrt_adjust()] for `OBJ` against `H1`.
#' @seealso [compute_plrt()], [lavTestLRT()]
#' @keywords internal
plrt_nested <- function(OBJ, H1) {
  stopifnot(inherits(H1, "latpl_fit"))
  if (
    !identical(OBJ$data$freq, H1$data$freq) ||
      !identical(OBJ$data$ufreq, H1$data$ufreq) ||
      !identical(OBJ$data$weights, H1$data$weights) ||
      !identical(OBJ$data$cluster, H1$data$cluster) ||
      H1$dims$npar <= OBJ$dims$npar
  ) {
    stop('OBJ is not nested in H1.')
  }
  pl <- function(FIT, THETA) cpp_pl(THETA, FIT$engine$model, FIT$engine$data)$pl
  lc0 <- OBJ$lincon
  lc1 <- H1$lincon
  m1 <- if (inherits(H1$model, "latpl_groups")) {
    H1$model$MODELS
  } else {
    list(H1$model)
  }
  emb <- function(x) {
    par <- theta2par(OBJ$model, as.vector(lc0$K %*% x + lc0$c))
    if (!inherits(OBJ$model, "latpl_groups")) {
      par <- list(par)
    }
    par <- Map(
      function(par, m) {
        q <- ncol(par$loadings)
        Lt <- par$loadings %*% solve(diag(q) - par$regressions)
        theta <- m$PARAMETERIZATION == "theta"
        if (
          all(m$CONSTRINT %in% 0) &&
            all(m$CONSTRMEAN %in% 0) &&
            all(m$CONSTRLOGSCALE %in% 0)
        ) {
          R <- par$residual_covariances
          s <- if (theta) {
            1 / sqrt(diag(R))
          } else if (is.null(par$scales)) {
            1 / sqrt(rowSums((Lt %*% par$latent_correlations) * Lt) + diag(R))
          } else {
            par$scales
          }
          k <- OBJ$dims$cat - 1
          par$thresholds <- rep(s, k) *
            (par$thresholds -
              rep(as.vector(par$intercepts + Lt %*% par$means), k))
          par$loadings <- s * par$loadings
          par$residual_covariances <- s * R * rep(s, each = length(s))
          par$intercepts[] <- 0
          par$means[] <- 0
          if (!theta) {
            diag(par$residual_covariances) <- 0
            par$scales <- rep(1, length(s))
          }
        } else if (!anyNA(m$CONSTRMEAN)) {
          par$intercepts <- as.vector(
            par$intercepts + Lt %*% (par$means - m$CONSTRMEAN)
          )
          par$means <- m$CONSTRMEAN
        }
        if (!length(m$REGMAT)) {
          A <- solve(diag(q) - par$regressions)
          S <- A %*% par$latent_correlations %*% t(A)
          mk <- cbind(
            apply(!is.na(m$CONSTRMAT) & m$CONSTRMAT != 0, 2, match, x = TRUE),
            seq_len(q)
          )
          sc <- ifelse(
            is.finite(m$CONSTRLOGSD),
            sqrt(diag(S)) / exp(m$CONSTRLOGSD),
            ifelse(is.na(mk[, 1]), 1, m$CONSTRMAT[mk] / par$loadings[mk])
          )
          par$loadings <- par$loadings %*% diag(sc, q)
          par$latent_correlations <- S / tcrossprod(sc)
          par$means <- par$means / sc
        }
        par
      },
      par,
      m1
    )
    par2theta(
      H1$model,
      if (inherits(H1$model, "latpl_groups")) par else par[[1]]
    )[lc1$free_idx]
  }
  theta0f <- OBJ$theta[lc0$free_idx]
  theta1 <- as.vector(lc1$K %*% emb(theta0f) + lc1$c)
  pl0 <- pl(OBJ, OBJ$theta)
  if (abs(pl(H1, theta1) - pl0) > 1e-8 * abs(pl0)) {
    stop('OBJ is not nested in H1.')
  }
  N <- numDeriv::jacobian(emb, theta0f)
  HJ <- cpp_HJ(theta1, H1$engine$model, H1$engine$data)
  Hf <- crossprod(lc1$K, HJ$H %*% lc1$K)
  M <- crossprod(lc1$K, HJ$J %*% lc1$K) %*%
    (solve(Hf) - N %*% solve(crossprod(N, Hf %*% N), t(N)))
  plrt_adjust(
    STAT = 2 * (pl(H1, H1$theta) - pl0),
    DF = H1$dims$npar - OBJ$dims$npar,
    E = sum(diag(M)),
    V = 2 * sum(M * t(M))
  )
}

#' Test against the saturated model
#'
#' @param OBJ Object of class latpl_fit.
#' @param VARS Output of [compute_var()] for `OBJ`.
#' @return A list with the tests of the model `model` and of the independence model `baseline` against the
#'   saturated model, and the saturated and independence fits `saturated` and `independence`.
#' @seealso [compute_plrt()], [plrt_adjust()]
#' @keywords internal
plrt_fit <- function(OBJ, VARS) {
  dims <- OBJ$dims
  n <- dims$n
  nthr <- sum(dims$cat - 1)
  P <- dims$p * (dims$p - 1) / 2
  eng <- OBJ$engine$data
  sat <- function(g, tau, rho, order, scores = FALSE) {
    cpp_saturated(tau, rho, eng[[g]], order, scores)
  }

  Y <- OBJ$data$Y
  v <- (rowSums(!is.na(Y)) - 1 + OBJ$data$uweights) *
    (if (length(OBJ$data$weights)) OBJ$data$weights else 1)
  fits <- lapply(seq_len(dims$ngroups), function(g) {
    rows <- OBJ$data$group == g - 1
    tau_ind <- unlist(lapply(seq_len(dims$p), function(j) {
      pj <- vapply(
        seq_len(dims$cat[j]) - 1,
        function(a) sum(v[rows & Y[, j] %in% a]),
        numeric(1)
      )
      stats::qnorm(cumsum(pj)[-dims$cat[j]] / sum(pj))
    }))
    tau <- tau_ind
    rho <- rep(0, P)
    o <- sat(g, tau, rho, 2)
    for (iter in 1:100) {
      if (max(abs(o$gradient)) / n < 1e-8) {
        break
      }
      gt <- o$gradient[1:nthr]
      gr <- o$gradient[-(1:nthr)]
      dt <- as.vector(solve(o$hschur, gt - o$htr %*% (gr / o$hrr)))
      dr <- as.vector(gr - crossprod(o$htr, dt)) / o$hrr
      step <- 1
      repeat {
        if (all(abs(rho + step * dr) < 1)) {
          o_new <- sat(g, tau + step * dt, rho + step * dr, 2)
          if (isTRUE(o_new$ll >= o$ll - 1e-12 * abs(o$ll))) break
        }
        step <- step / 2
        if (step < 1e-10) stop('Saturated model fit failed.')
      }
      tau <- tau + step * dt
      rho <- rho + step * dr
      o <- o_new
    }
    if (max(abs(o$gradient)) / n >= 1e-8) {
      warning('Saturated model fit did not converge.')
    }
    list(tau = tau, rho = rho, ll = o$ll, iter = iter, tau_ind = tau_ind)
  })
  tau_ind <- unlist(lapply(fits, function(x) x$tau_ind))
  ll_sat <- vapply(fits, function(x) x$ll, numeric(1))

  sat_terms <- function(TAU, RHO, DELTA = NULL, FULL = FALSE) {
    out <- list(ll = numeric(0), e = 0, v = 0, dad = 0)
    for (g in seq_len(dims$ngroups)) {
      it <- (g - 1) * nthr + seq_len(nthr)
      ir <- (g - 1) * P + seq_len(P)
      s <- sat(g, TAU[it], RHO[ir], 2, TRUE)
      htt <- s$htt / n
      htr <- s$htr / n
      d <- s$hrr / n
      w <- function(S) {
        if (length(eng[[g]]$WM)) {
          S %*% eng[[g]]$WM
        } else {
          S * rep(eng[[g]]$WV, each = nrow(S))
        }
      }
      tr <- cpp_saturated_traces(
        htt,
        htr,
        d,
        s$hschur / n,
        w(s$st),
        w(s$sr),
        n,
        FULL
      )
      out$ll <- c(out$ll, s$ll)
      out$e <- out$e + tr[1]
      out$v <- out$v + 2 * tr[2]
      if (!is.null(DELTA) && FULL) {
        Dt <- DELTA$tau[it, , drop = FALSE]
        Dr <- DELTA$rho[ir, , drop = FALSE]
        ht <- crossprod(Dt, htr %*% Dr)
        out$dad <- out$dad +
          crossprod(Dt, htt %*% Dt) +
          ht +
          t(ht) +
          crossprod(Dr, Dr * d)
      } else if (!is.null(DELTA)) {
        Dr <- DELTA[ir, , drop = FALSE]
        hd <- htr %*% Dr
        out$dad <- out$dad +
          crossprod(Dr, Dr * d) -
          crossprod(hd, solve(htt, hd))
      }
    }
    out
  }

  lc <- OBJ$lincon
  r <- cpp_rho(OBJ$theta, OBJ$engine$model, eng)
  jac <- if (ncol(lc$K) < nrow(lc$K)) r$jacobian %*% lc$K else r$jacobian
  invH <- solve(VARS$H)
  tpos <- which(theta_block(OBJ$model, nthr) == "tau")
  if (all(tpos %in% lc$free_idx)) {
    psi <- which(!(lc$free_idx %in% tpos))
    G <- (invH %*% VARS$J %*% invH)[psi, psi, drop = FALSE]
    P0 <- solve(invH[psi, psi, drop = FALSE], G)
    tz <- sat_terms(r$tau, r$rho, jac[, psi, drop = FALSE])
    cv <- sum(diag(tz$dad %*% G %*% P0))
  } else {
    P0 <- invH %*% VARS$J
    tz <- sat_terms(
      r$tau,
      r$rho,
      list(tau = r$tau_jacobian %*% lc$K, rho = jac),
      FULL = TRUE
    )
    cv <- sum(diag(tz$dad %*% P0 %*% P0 %*% invH))
  }
  model <- plrt_adjust(
    STAT = 2 * (ll_sat - tz$ll),
    DF = dims$nthr + dims$pairs - dims$npar,
    E = tz$e - sum(diag(P0)),
    V = tz$v + 2 * sum(P0 * t(P0)) - 4 * cv
  )

  # Independence model against saturated
  tb <- sat_terms(tau_ind, rep(0, dims$pairs))
  baseline <- plrt_adjust(
    STAT = 2 * (ll_sat - tb$ll),
    DF = dims$pairs,
    E = tb$e,
    V = tb$v
  )

  list(
    model = model,
    baseline = baseline,
    saturated = list(
      tau = unlist(lapply(fits, function(x) x$tau)),
      rho = unlist(lapply(fits, function(x) x$rho)),
      logl = sum(ll_sat),
      iterations = max(vapply(fits, function(x) x$iter, numeric(1)))
    ),
    independence = list(tau = tau_ind, logl = sum(tb$ll))
  )
}
