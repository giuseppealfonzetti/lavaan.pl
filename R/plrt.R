#' Pairwise likelihood ratio test
#'
#' @param OBJ Object of class `lavaanpl_internal` from [fit_lavaanpl()].
#' @param VARS Output of [compute_var()] for `OBJ`.
#' @param H1 Object of class `lavaanpl_internal` in which `OBJ` is nested, or `NULL` to test `OBJ` against the
#'   saturated model.
#' @return A list with the test `model` of `OBJ` against `H1` or the saturated model, as returned by [plrt_adjust()],
#'   and the timing `RTime`. Without `H1`, also the test `baseline` of the independence model against the saturated
#'   model and the fits `saturated` and `independence`.
#' @seealso [plrt_nested()], [plrt_fit()]
#' @keywords internal
compute_plrt <- function(OBJ, VARS = compute_var(OBJ), H1 = NULL) {
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
    pvalue.scaled = if (isTRUE(E > 0 && V > 0)) {
      stats::pchisq(c1 * sum(STAT), E^2 / (V / 2), lower.tail = FALSE)
    } else {
      NA
    }
  )
}


#' Test of nested models
#'
#' @param OBJ Object of class lavaanpl_internal, nested in `H1`.
#' @param H1 Object of class lavaanpl_internal.
#' @return The output of [plrt_adjust()] for `OBJ` against `H1`.
#' @seealso [compute_plrt()], [lavTestLRT()]
#' @keywords internal
plrt_nested <- function(OBJ, H1) {
  stopifnot(inherits(H1, "lavaanpl_internal"))
  if (
    !identical(OBJ$data$freq, H1$data$freq) ||
      !identical(OBJ$data$Y, H1$data$Y) ||
      !identical(OBJ$data$X, H1$data$X) ||
      !identical(OBJ$data$ufreq, H1$data$ufreq) ||
      !identical(OBJ$data$weights, H1$data$weights) ||
      !identical(OBJ$data$cluster, H1$data$cluster) ||
      H1$dims$npar <= OBJ$dims$npar
  ) {
    stop('OBJ is not nested in H1.')
  }
  pl <- function(FIT, PARAMETER_VECTOR) cpp_pl(PARAMETER_VECTOR, FIT$cpp_inputs$model, FIT$cpp_inputs$data)$pl
  lc0 <- OBJ$lincon
  lc1 <- H1$lincon
  emb <- function(x) {
    par <- vector_to_matrices(OBJ$model, as.vector(lc0$K %*% x + lc0$c))
    par <- Map(
      function(par, m) {
        q <- ncol(par$lambda)
        Lt <- par$lambda %*% solve(diag(q) - par$beta)
        nx <- ncol(OBJ$data$X)
        if (is.null(par$gamma) && nx > 0) {
          par$gamma <- matrix(0, q, nx)
          par$kappa <- matrix(0, nrow(par$lambda), nx)
        }
        thetapar <- m$PARAMETERIZATION == "theta"
        if (
          all(m$NU %in% 0) &&
            all(m$ALPHA %in% 0) &&
            all(item_scales(m) %in% 1)
        ) {
          R <- par$theta
          s <- if (thetapar) {
            1 / sqrt(diag(R))
          } else if (is.null(par$delta)) {
            1 / sqrt(rowSums((Lt %*% par$psi) * Lt) + diag(R))
          } else {
            par$delta
          }
          k <- OBJ$data$cat - 1
          par$tau <- rep(s, k) *
            (par$tau -
              rep(as.vector(par$nu + Lt %*% par$alpha), k))
          par$lambda <- s * par$lambda
          if (nx > 0) par$kappa <- s * par$kappa
          par$theta <- s * R * rep(s, each = length(s))
          par$nu[] <- 0
          par$alpha[] <- 0
          if (!thetapar) {
            diag(par$theta) <- 0
            par$delta <- rep(1, length(s))
          }
        } else if (!anyNA(m$ALPHA)) {
          par$nu <- as.vector(
            par$nu + Lt %*% (par$alpha - m$ALPHA)
          )
          par$alpha <- m$ALPHA
        }
        if (all(m$BETA %in% 0)) {
          A <- solve(diag(q) - par$beta)
          S <- A %*% par$psi %*% t(A)
          mk <- cbind(
            apply(!is.na(m$LAMBDA) & m$LAMBDA != 0, 2, match, x = TRUE),
            seq_len(q)
          )
          sc <- ifelse(
            !is.na(diag(m$PSI)),
            sqrt(diag(S) / diag(m$PSI)),
            ifelse(is.na(mk[, 1]), 1, m$LAMBDA[mk] / par$lambda[mk])
          )
          par$lambda <- par$lambda %*% diag(sc, q)
          par$psi <- S / tcrossprod(sc)
          par$alpha <- par$alpha / sc
          if (nx > 0) par$gamma <- (A %*% par$gamma) / sc
        }
        par
      },
      par,
      H1$model$MODELS
    )
    matrices_to_vector(H1$model, par)[lc1$free_idx]
  }
  free0 <- OBJ$parameter_vector[lc0$free_idx]
  vector1 <- as.vector(lc1$K %*% emb(free0) + lc1$c)
  pl0 <- pl(OBJ, OBJ$parameter_vector)
  if (abs(pl(H1, vector1) - pl0) > 1e-8 * abs(pl0)) {
    stop('OBJ is not nested in H1.')
  }
  N <- numDeriv::jacobian(emb, free0)
  HJ <- cpp_sandwich(vector1, H1$cpp_inputs$model, H1$cpp_inputs$data)
  Hf <- crossprod(lc1$K, HJ$H %*% lc1$K)
  M <- crossprod(lc1$K, HJ$J %*% lc1$K) %*%
    (solve(Hf) - N %*% solve(crossprod(N, Hf %*% N), t(N)))
  plrt_adjust(
    STAT = 2 * (pl(H1, H1$parameter_vector) - pl0),
    DF = H1$dims$npar - OBJ$dims$npar,
    E = sum(diag(M)),
    V = 2 * sum(M * t(M))
  )
}

#' Test against the saturated model
#'
#' @param OBJ Object of class lavaanpl_internal.
#' @param VARS Output of [compute_var()] for `OBJ`.
#' @return A list with the tests of the model `model` and of the independence model `baseline` against the
#'   saturated model, and the saturated and independence fits `saturated` and `independence`.
#' @seealso [compute_plrt()], [plrt_adjust()]
#' @keywords internal
plrt_fit <- function(OBJ, VARS) {
  dims <- OBJ$dims
  data <- OBJ$data
  n <- data$n
  nx <- ncol(data$X)
  nnum <- sum(data$cat == 1)
  nthr <- sum(data$cat - 1) + data$p * nx + 2 * nnum
  P <- data$pairs / data$ngroups
  inputs <- OBJ$cpp_inputs$data
  sat <- function(g, tau, rho, order, scores = FALSE) {
    cpp_saturated(tau, rho, inputs[[g]], order, scores)
  }

  Y <- OBJ$data$Y
  v <- (rowSums(!is.na(Y)) - 1 + OBJ$data$uweights) *
    (if (length(OBJ$data$weights)) OBJ$data$weights else 1)
  fits <- lapply(seq_len(data$ngroups), function(g) {
    rows <- OBJ$data$group == g - 1
    tau_ind <- unlist(lapply(seq_len(data$p), function(j) {
      pj <- vapply(
        seq_len(data$cat[j]) - 1,
        function(a) sum(v[rows & Y[, j] %in% a]),
        numeric(1)
      )
      stats::qnorm(cumsum(pj)[-data$cat[j]] / sum(pj))
    }))
    numeric_items <- which(data$cat == 1)
    centre <- vapply(numeric_items, function(j) {
      seen <- rows & !is.na(Y[, j])
      sum(v[seen] * Y[seen, j]) / sum(v[seen])
    }, numeric(1))
    spread <- vapply(seq_along(numeric_items), function(a) {
      seen <- rows & !is.na(Y[, numeric_items[a]])
      sqrt(sum(v[seen] * (Y[seen, numeric_items[a]] - centre[a])^2) / sum(v[seen]))
    }, numeric(1))
    tau_ind <- c(tau_ind, rep(0, data$p * nx), -centre / spread, 1 / spread)
    for (iter in seq_len(if (nx > 0 || nnum > 0) 100 else 0)) {
      o <- sat(g, tau_ind, rep(0, P), 2)
      if (max(abs(o$gradient[1:nthr])) / n < 1e-8) {
        break
      }
      dt <- as.vector(solve(o$htt, o$gradient[1:nthr]))
      step <- 1
      repeat {
        o_new <- sat(g, tau_ind + step * dt, rep(0, P), 0)
        if (isTRUE(o_new$ll >= o$ll - 1e-12 * abs(o$ll))) break
        step <- step / 2
        if (step < 1e-10) stop('Independence model fit failed.')
      }
      tau_ind <- tau_ind + step * dt
    }
    tau <- tau_ind
    R <- stats::cor(Y[rows, ], use = "pairwise")
    rho <- R[upper.tri(R)]
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
    for (g in seq_len(data$ngroups)) {
      it <- (g - 1) * nthr + seq_len(nthr)
      ir <- (g - 1) * P + seq_len(P)
      s <- sat(g, TAU[it], RHO[ir], 2, TRUE)
      htt <- s$htt / n
      htr <- s$htr / n
      d <- s$hrr / n
      w <- function(S) {
        if (length(inputs[[g]]$CLUSTER_WEIGHTS)) {
          S %*% inputs[[g]]$CLUSTER_WEIGHTS
        } else {
          S * rep(inputs[[g]]$JWEIGHTS, each = nrow(S))
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
  r <- cpp_reduced_form(OBJ$parameter_vector, OBJ$cpp_inputs$model, inputs)
  jac <- if (ncol(lc$K) < nrow(lc$K)) r$jacobian %*% lc$K else r$jacobian
  invH <- solve(VARS$H)
  tpos <- which(OBJ$blocks == "tau")
  if (nx == 0 && nnum == 0 && all(tpos %in% lc$free_idx)) {
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
    DF = dims$nthr + data$ngroups * (data$p * nx + 2 * nnum) + data$pairs - dims$npar,
    E = tz$e - sum(diag(P0)),
    V = tz$v + 2 * sum(P0 * t(P0)) - 4 * cv
  )

  # Independence model against saturated
  tb <- sat_terms(tau_ind, rep(0, data$pairs))
  baseline <- plrt_adjust(
    STAT = 2 * (ll_sat - tb$ll),
    DF = data$pairs,
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
