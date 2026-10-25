#' Pairwise likelihood ratio test
#' @param OBJ Object of class PlFaFit
#' @param VARS Output of [computeVar()] for `OBJ`
#' @param H1 Object of class PlFaFit in which `OBJ` is nested, or `NULL` to test `OBJ` against the
#'   saturated model
#' @return A list with the test of the model against `H1` or the saturated model `model` and, without
#'   `H1`, of the independence model against the saturated model `baseline`. Each holds the
#'   statistic `stat`, its degrees of freedom `df`, and the mean and variance adjusted `stat.scaled`
#'   (`stat` divided by `scaling.factor`), `df.scaled` and `pvalue.scaled`. Also, without `H1`, the
#'   saturated and independence fits `saturated` and `independence`, and the timing `RTime`.
#' @export
computePLRT <- function(OBJ, VARS = computeVar(OBJ), H1 = NULL) UseMethod('computePLRT')
#' Pairwise likelihood ratio test
#' @param OBJ Object of class PlFaFit
#' @param VARS Output of [computeVar()] for `OBJ`
#' @param H1 Object of class PlFaFit in which `OBJ` is nested, or `NULL` to test `OBJ` against the
#'   saturated model
#' @return A list with the test of the model against `H1` or the saturated model `model` and, without
#'   `H1`, of the independence model against the saturated model `baseline`. Each holds the
#'   statistic `stat`, its degrees of freedom `df`, and the mean and variance adjusted `stat.scaled`
#'   (`stat` divided by `scaling.factor`), `df.scaled` and `pvalue.scaled`. Also, without `H1`, the
#'   saturated and independence fits `saturated` and `independence`, and the timing `RTime`.
#' @export
computePLRT.PlFaFit <- function(OBJ, VARS = computeVar(OBJ), H1 = NULL){

  start_time <- Sys.time()
  dims <- OBJ$dims
  n <- dims$n
  nthr <- dims$nthr
  RcppParallel::setThreadOptions(numThreads = OBJ$cores)

  adjust <- function(STAT, DF, E, V){
    c1 <- E / (V / 2)
    list(stat = STAT, df = DF, stat.scaled = c1 * STAT, df.scaled = E^2 / (V / 2), scaling.factor = 1 / c1,
         pvalue.scaled = if(isTRUE(c1 > 0)) stats::pchisq(c1 * STAT, E^2 / (V / 2), lower.tail = FALSE) else NA)
  }

  # Nested models: H1 evaluated at the H0 estimates, with the exact restriction
  # given by the Jacobian of the embedding of OBJ into H1
  if(!is.null(H1)){
    stopifnot(inherits(H1, "PlFaFit"))
    if(!identical(OBJ$data$freq, H1$data$freq) || !identical(OBJ$data$ufreq, H1$data$ufreq) || H1$dims$npar <= dims$npar) {
      stop('OBJ is not nested in H1.')
    }
    nll <- function(FIT, THETA){
      cpp_multiThread_completePairwise(
        N = n, C_VEC = dims$cat, CONSTRMAT = FIT$model$CONSTRMAT, CONSTRLOGSD = FIT$model$CONSTRLOGSD,
        LLC = FIT$model$LLC, THETA = THETA, FREQ = FIT$data$freq, CORRFLAG = FIT$model$CORRFLAG,
        NTHR = FIT$dims$nthr, NLOAD = FIT$dims$nload, NCORR = FIT$dims$ncorr, NVAR = FIT$dims$nvar,
        UFREQ = FIT$data$ufreq
      )$iter_nll
    }
    emb <- function(theta){
      gp <- getPar(replace(OBJ, "theta", list(theta)), OPTION = "list")
      c(gp$thresholds,
        cpp_loadings_mat2vec(gp$loadings, H1$model$CONSTRMAT, H1$dims$nload),
        cpp_latvar_mat2vec(gp$latent_correlations, H1$model$CONSTRLOGSD, H1$dims$ncorr, H1$dims$nvar))
    }
    theta1 <- emb(OBJ$theta)
    nll0 <- nll(OBJ, OBJ$theta)
    if(abs(nll(H1, theta1) - nll0) > 1e-8 * abs(nll0)) stop('OBJ is not nested in H1.')
    K <- numDeriv::jacobian(emb, OBJ$theta)
    HJ <- cpp_sample_estimators_HJ(
      THETA = theta1, FREQ = H1$data$freq, DATA = H1$data$Y, C_VEC = H1$dims$cat,
      CONSTRMAT = H1$model$CONSTRMAT, CONSTRLOGSD = H1$model$CONSTRLOGSD, LLC = H1$model$LLC, N = n,
      CORRFLAG = H1$model$CORRFLAG, NTHR = H1$dims$nthr, NLOAD = H1$dims$nload, NCORR = H1$dims$ncorr,
      NVAR = H1$dims$nvar, UFREQ = H1$data$ufreq, UW = H1$data$uweights
    )
    M <- HJ$J %*% (solve(HJ$H) - K %*% solve(crossprod(K, HJ$H %*% K), t(K)))
    model <- adjust(STAT = 2 * (nll0 - nll(H1, H1$theta)),
                    DF   = H1$dims$npar - dims$npar,
                    E    = sum(diag(M)),
                    V    = 2 * sum(M * t(M)))
    end_time <- Sys.time()
    return(list(model = model, RTime = as.numeric(difftime(end_time, start_time, units = 'secs')[1])))
  }

  # Unique response patterns and their counts
  Y <- OBJ$data$Y
  Yu <- unique(Y)
  key <- match(do.call(paste, as.data.frame(Y)), do.call(paste, as.data.frame(Yu)))
  w <- as.numeric(tabulate(key))
  uw <- OBJ$data$uweights[match(seq_len(nrow(Yu)), key)]

  sat <- function(tau, rho, order, scores = FALSE){
    cpp_saturated(tau, rho, OBJ$data$freq, Yu, dims$cat, order, scores, OBJ$data$ufreq, uw)
  }

  # Independence model: thresholds from the univariate margins, with unit weights
  # given by the number of pairs and univariate terms each observed item enters
  v <- rowSums(!is.na(Y)) - 1 + OBJ$data$uweights
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
    B <- o$htt - o$htr %*% (t(o$htr) / o$hrr)
    dt <- as.vector(solve(B, gt - o$htr %*% (gr / o$hrr)))
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
  # and scores of the unique response patterns
  sat_terms <- function(TAU, RHO, DELTA = NULL){
    s <- sat(TAU, RHO, 2, TRUE)
    htt <- s$htt / n
    htr <- s$htr / n
    d <- s$hrr / n
    tr <- cpp_saturated_traces(htt, htr, d, s$st, s$sr, w, n)
    out <- list(ll = s$ll, e = tr[1], v = 2 * tr[2])
    if(!is.null(DELTA)){
      hd <- htr %*% DELTA
      out$dad <- crossprod(DELTA, DELTA * d) - crossprod(hd, solve(htt, hd))
    }
    out
  }

  # Model against saturated
  psi <- seq_len(dims$npar)[-(1:nthr)]
  rho_theta <- function(theta){
    gp <- getPar(replace(OBJ, "theta", list(theta)), OPTION = "list")
    S <- gp$loadings %*% gp$latent_correlations %*% t(gp$loadings)
    S[upper.tri(S)]
  }
  delta <- numDeriv::jacobian(rho_theta, OBJ$theta)[, psi, drop = FALSE]
  invH <- solve(VARS$H)
  G <- (invH %*% VARS$J %*% invH)[psi, psi, drop = FALSE]
  P0 <- solve(invH[psi, psi, drop = FALSE], G)
  tz <- sat_terms(OBJ$theta[1:nthr], rho_theta(OBJ$theta), delta)
  model <- adjust(STAT = 2 * (o$ll - tz$ll),
                  DF   = nthr + dims$pairs - dims$npar,
                  E    = tz$e - sum(diag(P0)),
                  V    = tz$v + 2 * sum(P0 * t(P0)) - 4 * sum(diag(tz$dad %*% G %*% P0)))

  # Independence model against saturated
  tb <- sat_terms(tau_ind, rep(0, dims$pairs))
  baseline <- adjust(STAT = 2 * (o$ll - tb$ll), DF = dims$pairs, E = tb$e, V = tb$v)

  end_time <- Sys.time()
  list(
    model = model,
    baseline = baseline,
    saturated = list(tau = tau, rho = rho, logl = o$ll, iterations = iter),
    independence = list(tau = tau_ind, logl = tb$ll),
    RTime = as.numeric(difftime(end_time, start_time, units = 'secs')[1])
  )
}
