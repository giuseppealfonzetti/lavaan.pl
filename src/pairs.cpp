#include "pairs.h"

pairs::Setup::Setup(
    const Eigen::Ref<const Eigen::VectorXd> THETA,
    const Eigen::Ref<const Eigen::MatrixXd> CONSTRMAT,
    const Eigen::Ref<const Eigen::VectorXd> CONSTRLOGSD,
    const int NTHR,
    const int NLOAD,
    const int NCORR,
    const int NVAR,
    const int CORRFLAG
): p(CONSTRMAT.rows()), q(CONSTRMAT.cols()), nthr(NTHR), nload(NLOAD), ncorr(NCORR), nvar(NVAR), m(NLOAD + NCORR + NVAR),
   direct(CORRFLAG == 2){

  Lam = params::loadings::theta2mat(THETA, CONSTRMAT, NTHR, NLOAD);
  sd  = params::latvar::theta2dmat(THETA, CONSTRLOGSD, NTHR, NLOAD, NCORR, NVAR, q).diagonal();
  tau = params::thresholds::theta2vec(THETA, NTHR);

  Z = Eigen::MatrixXd::Zero(q, q);
  ra.resize(ncorr);
  ca.resize(ncorr);
  for(int c = 1, a = 0; c < q && a < ncorr; c++){
    for(int r = 0; r < c; r++, a++){
      Z(r, c) = std::tanh(THETA(NTHR + NLOAD + a));
      ra[a] = r;
      ca[a] = c;
    }
  }
  if(direct){
    Sig = sd.array().square().matrix().asDiagonal();
    for(int a = 0; a < ncorr; a++){
      Sig(ra[a], ca[a]) = THETA(NTHR + NLOAD + a);
      Sig(ca[a], ra[a]) = THETA(NTHR + NLOAD + a);
    }
    Z.setZero();
  }
  C = (1 - Z.array().square()).sqrt().matrix();
  P = Eigen::MatrixXd::Ones(q, q);
  U = Eigen::MatrixXd::Zero(q, q);
  for(int c = 0; c < q; c++){
    for(int r = 1; r <= c; r++) P(r, c) = P(r - 1, c) * C(r - 1, c);
    for(int r = 0; r < c; r++) U(r, c) = Z(r, c) * P(r, c);
    U(c, c) = P(c, c);
  }

  alpha = Eigen::MatrixXd::Zero(q, ncorr);
  for(int a = 0; a < ncorr; a++){
    const int r = ra[a], c = ca[a];
    alpha(r, a) = P(r, c) * C(r, c) * C(r, c);
    for(int i = r + 1; i <= c; i++) alpha(i, a) = -U(i, c) * Z(r, c);
  }

  for(int h = 0; h < q; h++) if(!std::isfinite(CONSTRLOGSD(h))) sfree.push_back(h);

  JL.assign(p, Eigen::MatrixXd::Zero(q, nload));
  for(int h = 0, idx = 0; h < q; h++){
    for(int j = 0; j < p; j++){
      if(!std::isfinite(CONSTRMAT(j, h))){
        JL[j](h, idx) = 1;
        idx++;
      }
    }
  }
}

pairs::Pair::Pair(
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const int K,
    const int L,
    const double RHO
): k(K), l(L), ck(C_VEC(K)), cl(C_VEC(L)), rho(RHO){

  // Identifies quantities related to pair (k,l)
  tk = C_VEC.head(K).sum() - K;
  tl = C_VEC.head(L).sum() - L;
  mt = ck + cl - 2;

  // identify column index in freq table
  // i1: starting index item k
  offset = 0;
  for(int u = 1; u < K; u++) offset += C_VEC(u) * C_VEC.head(u).sum();
  // i2 starting index from i1 for item l
  offset += ck * C_VEC.head(L).sum();
}

pairs::Pair::Pair(
    const Setup &S,
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const int K,
    const int L,
    const bool GRADFLAG
): Pair(C_VEC, K, L, 0){

  if(S.direct){
    x = S.Lam.row(K).transpose();
    y = S.Lam.row(L).transpose();
    const Eigen::VectorXd Sx = S.Sig * x, Sy = S.Sig * y;
    rho = x.dot(Sy);
    if(GRADFLAG){
      grho.resize(S.m);
      grho.head(S.nload) = S.JL[K].transpose() * Sy + S.JL[L].transpose() * Sx;
      for(int a = 0; a < S.ncorr; a++) grho(S.nload + a) = x(S.ra[a]) * y(S.ca[a]) + x(S.ca[a]) * y(S.ra[a]);
      for(int s = 0; s < S.nvar; s++){
        const int h = S.sfree[s];
        grho(S.nload + S.ncorr + s) = 2 * S.sd(h) * S.sd(h) * x(h) * y(h);
      }
    }
    return;
  }

  x = S.sd.cwiseProduct(S.Lam.row(K).transpose());
  y = S.sd.cwiseProduct(S.Lam.row(L).transpose());
  X = S.U * x;
  Y = S.U * y;
  rho = X.dot(Y);

  if(GRADFLAG){
    grho.resize(S.m);
    const Eigen::VectorXd dUY = S.sd.cwiseProduct(S.U.transpose() * Y);
    const Eigen::VectorXd dUX = S.sd.cwiseProduct(S.U.transpose() * X);
    grho.head(S.nload) = S.JL[K].transpose() * dUY + S.JL[L].transpose() * dUX;
    for(int a = 0; a < S.ncorr; a++){
      grho(S.nload + a) = x(S.ca[a]) * S.alpha.col(a).dot(Y) + y(S.ca[a]) * S.alpha.col(a).dot(X);
    }
    for(int s = 0; s < S.nvar; s++){
      const int h = S.sfree[s];
      grho(S.nload + S.ncorr + s) = x(h) * S.U.col(h).dot(Y) + y(h) * S.U.col(h).dot(X);
    }
  }
}

double pairs::pattern_derivatives(
    const Eigen::Ref<const Eigen::VectorXd> TAU,
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const Pair &PR,
    const int SK,
    const int SL,
    const int ORDER,
    Eigen::VectorXd &G1,
    Eigen::MatrixXd &G2
){
  const int mt = PR.mt;
  const double rho = PR.rho;

  // identify thresholds
  const Eigen::VectorXd th = params::extract_thresholds(TAU, C_VEC, PR.k, PR.l, SK, SL);

  // compute pi
  const double pe = biprobs::compute_pi(C_VEC, th, rho, PR.k, PR.l, SK, SL) + 1e-8;
  if(ORDER == 0) return pe;

  const double r2 = 1 - rho * rho, sr = std::sqrt(r2);
  G1 = Eigen::VectorXd::Zero(mt + 1);
  if(ORDER == 2) G2 = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
  for(int e = 0; e < 2; e++){
    for(int f = 0; f < 2; f++){
      const double xv = th(2 * e), yv = th(2 * f + 1);
      const int ik = SK - e, il = SL - f;
      const int px = (ik >= 0 && ik < PR.ck - 1) ? ik : -1;
      const int py = (il >= 0 && il < PR.cl - 1) ? PR.ck - 1 + il : -1;
      const double sg = (e == f) ? 1 : -1;
      const double wx = (yv - rho * xv) / sr, wy = (xv - rho * yv) / sr;
      const double phx = R::dnorm(xv, 0, 1, 0), phy = R::dnorm(yv, 0, 1, 0);
      const double Fx = phx * R::pnorm(wx, 0, 1, 1, 0), Fy = phy * R::pnorm(wy, 0, 1, 1, 0);
      const double f2 = binorm::dbvnorm(xv, yv, rho, 0);
      if(px >= 0) G1(px) += sg * Fx;
      if(py >= 0) G1(py) += sg * Fy;
      G1(mt) += sg * f2;
      if(ORDER == 2){
        if(px >= 0){
          G2(px, px) += sg * (-xv * Fx - rho / sr * phx * R::dnorm(wx, 0, 1, 0));
          G2(px, mt) += -sg * f2 * (xv - rho * yv) / r2;
          G2(mt, px) += -sg * f2 * (xv - rho * yv) / r2;
        }
        if(py >= 0){
          G2(py, py) += sg * (-yv * Fy - rho / sr * phy * R::dnorm(wy, 0, 1, 0));
          G2(py, mt) += -sg * f2 * (yv - rho * xv) / r2;
          G2(mt, py) += -sg * f2 * (yv - rho * xv) / r2;
        }
        if(px >= 0 && py >= 0){
          G2(px, py) += sg * f2;
          G2(py, px) += sg * f2;
        }
        G2(mt, mt) += sg * f2 * (rho / r2 + xv * yv / r2 - rho * (xv * xv - 2 * rho * xv * yv + yv * yv) / (r2 * r2));
      }
    }
  }
  return pe;
}

Eigen::MatrixXd pairs::rho_hessian(const Setup &S, const Pair &PR){
  const int q = S.q, nload = S.nload, ncorr = S.ncorr, nvar = S.nvar;
  const Eigen::VectorXd &x = PR.x, &y = PR.y, &X = PR.X, &Y = PR.Y, &grho = PR.grho;

  if(S.direct){
    const Eigen::MatrixXd &Jk = S.JL[PR.k], &Jl = S.JL[PR.l];
    Eigen::MatrixXd hrho = Eigen::MatrixXd::Zero(S.m, S.m);
    const Eigen::MatrixXd ll = Jk.transpose() * S.Sig * Jl;
    hrho.topLeftCorner(nload, nload) = ll + ll.transpose();
    for(int a = 0; a < ncorr; a++){
      Eigen::VectorXd ek = Eigen::VectorXd::Zero(q), el = Eigen::VectorXd::Zero(q);
      ek(S.ra[a]) += y(S.ca[a]); ek(S.ca[a]) += y(S.ra[a]);
      el(S.ra[a]) += x(S.ca[a]); el(S.ca[a]) += x(S.ra[a]);
      const Eigen::VectorXd col = Jk.transpose() * ek + Jl.transpose() * el;
      hrho.block(0, nload + a, nload, 1) = col;
      hrho.block(nload + a, 0, 1, nload) = col.transpose();
    }
    for(int s = 0; s < nvar; s++){
      const int h = S.sfree[s], js = nload + ncorr + s;
      const double v = 2 * S.sd(h) * S.sd(h);
      const Eigen::VectorXd col = v * (Jk.row(h).transpose() * y(h) + Jl.row(h).transpose() * x(h));
      hrho.block(0, js, nload, 1) = col;
      hrho.block(js, 0, 1, nload) = col.transpose();
      hrho(js, js) = 2 * v * x(h) * y(h);
    }
    return hrho;
  }
  const Eigen::MatrixXd DJk = S.sd.asDiagonal() * S.JL[PR.k], DJl = S.sd.asDiagonal() * S.JL[PR.l];

  Eigen::MatrixXd dX(q, S.m), dY(q, S.m);
  dX.leftCols(nload) = S.U * DJk;
  dY.leftCols(nload) = S.U * DJl;
  for(int a = 0; a < ncorr; a++){
    dX.col(nload + a) = x(S.ca[a]) * S.alpha.col(a);
    dY.col(nload + a) = y(S.ca[a]) * S.alpha.col(a);
  }
  for(int s = 0; s < nvar; s++){
    dX.col(nload + ncorr + s) = x(S.sfree[s]) * S.U.col(S.sfree[s]);
    dY.col(nload + ncorr + s) = y(S.sfree[s]) * S.U.col(S.sfree[s]);
  }
  Eigen::MatrixXd hrho = dX.transpose() * dY;
  hrho += hrho.transpose().eval();

  const Eigen::VectorXd aY = S.alpha.transpose() * Y, aX = S.alpha.transpose() * X;
  const Eigen::VectorXd uY = S.U.transpose() * Y, uX = S.U.transpose() * X;
  for(int i = 0; i < nload; i++){
    for(int a = 0; a < ncorr; a++){
      const double v = DJk(S.ca[a], i) * aY(a) + DJl(S.ca[a], i) * aX(a);
      hrho(i, nload + a) += v;
      hrho(nload + a, i) += v;
    }
    for(int s = 0; s < nvar; s++){
      const int h = S.sfree[s];
      const double v = DJk(h, i) * uY(h) + DJl(h, i) * uX(h);
      hrho(i, nload + ncorr + s) += v;
      hrho(nload + ncorr + s, i) += v;
    }
  }
  for(int a = 0; a < ncorr; a++){
    for(int b = a; b < ncorr && S.ca[b] == S.ca[a]; b++){
      const int c = S.ca[a], i0 = S.ra[a], i1 = S.ra[b];
      Eigen::VectorXd dd = Eigen::VectorXd::Zero(q);
      if(i0 < i1){
        dd(i1) = -S.P(i1, c) * S.Z(i0, c) * S.C(i1, c) * S.C(i1, c);
        for(int i = i1 + 1; i <= c; i++) dd(i) = S.U(i, c) * S.Z(i0, c) * S.Z(i1, c);
      }else{
        dd(i0) = -2 * S.P(i0, c) * S.Z(i0, c) * S.C(i0, c) * S.C(i0, c);
        for(int i = i0 + 1; i <= c; i++) dd(i) = S.U(i, c) * (2 * S.Z(i0, c) * S.Z(i0, c) - 1);
      }
      const double v = x(c) * dd.dot(Y) + y(c) * dd.dot(X);
      hrho(nload + a, nload + b) += v;
      if(b != a) hrho(nload + b, nload + a) += v;
    }
  }
  for(int s = 0; s < nvar; s++){
    const int js = nload + ncorr + s;
    hrho(js, js) += grho(js);
    for(int a = 0; a < ncorr; a++){
      if(S.ca[a] != S.sfree[s]) continue;
      hrho(nload + a, js) += grho(nload + a);
      hrho(js, nload + a) += grho(nload + a);
    }
  }
  return hrho;
}

void pairs::univariate_contribution(
    const Eigen::Ref<const Eigen::VectorXd> TAU,
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const Eigen::Ref<const Eigen::VectorXd> UFREQ,
    const int ORDER,
    double &LL,
    Eigen::VectorXd &GRADIENT,
    Eigen::MatrixXd &HESSIAN,
    Eigen::MatrixXd &USCORE
){
  if(UFREQ.size() == 0) return;
  for(int k = 0, uk = 0; k < C_VEC.size(); uk += C_VEC(k), k++){
    const int ck = C_VEC(k), tk = uk - k;
    for(int a = 0; a < ck; a++){
      const double nc = UFREQ(uk + a);
      const bool up = a < ck - 1, dw = a > 0;
      const double hi = up ? TAU(tk + a) : 0, lo = dw ? TAU(tk + a - 1) : 0;
      const double pe = (up ? R::pnorm(hi, 0, 1, 1, 0) : 1) - (dw ? R::pnorm(lo, 0, 1, 1, 0) : 0) + 1e-8;
      LL += nc * log(pe);
      if(ORDER == 0) continue;

      const double ghi = up ? R::dnorm(hi, 0, 1, 0) : 0, glo = dw ? -R::dnorm(lo, 0, 1, 0) : 0;
      if(up) GRADIENT(tk + a) += nc * ghi / pe;
      if(dw) GRADIENT(tk + a - 1) += nc * glo / pe;
      if(USCORE.size() > 0){
        if(up) USCORE(tk + a, uk + a) = ghi / pe;
        if(dw) USCORE(tk + a - 1, uk + a) = glo / pe;
      }
      if(ORDER == 2){
        if(up) HESSIAN(tk + a, tk + a) += nc * (ghi * ghi / (pe * pe) + hi * ghi / pe);
        if(dw) HESSIAN(tk + a - 1, tk + a - 1) += nc * (glo * glo / (pe * pe) + lo * glo / pe);
        if(up && dw){
          HESSIAN(tk + a, tk + a - 1) += nc * ghi * glo / (pe * pe);
          HESSIAN(tk + a - 1, tk + a) += nc * ghi * glo / (pe * pe);
        }
      }
    }
  }
}

// Single pair contribution to
// 1. Log-likelihood
// 2. Gradient
void pairs::pair_contribution(
    const Setup &S,
    const Eigen::Ref<const Eigen::VectorXd> C_VEC,
    const Eigen::Ref<const Eigen::MatrixXd> FREQ,
    const int K,
    const int L,
    const int GRADFLAG,
    double &LL,
    Eigen::VectorXd &GRADIENT
){
  const Pair PR(S, C_VEC, K, L, GRADFLAG == 1);
  const int order = GRADFLAG == 1 ? 1 : 0;
  Eigen::VectorXd g, gpair = Eigen::VectorXd::Zero(PR.mt + 1);
  Eigen::MatrixXd G2;

  ////////////////////////////
  /* LIKELIHOOD COMPUTATION */
  ////////////////////////////
  for(int sk = 0; sk < PR.ck; sk++){
    for(int sl = 0; sl < PR.cl; sl++){

      // read frequency
      const double n_sksl = FREQ(4, PR.offset + sk * PR.cl + sl);

      const double pe = pattern_derivatives(S.tau, C_VEC, PR, sk, sl, order, g, G2);

      // update ll
      LL += n_sksl * log(pe);
      if(GRADFLAG == 1) gpair += n_sksl / pe * g;
    }
  }

  //////////////////////////
  /* GRADIENT COMPUTATION */
  /////////////////////////
  if(GRADFLAG == 1){
    GRADIENT.segment(PR.tk, PR.ck - 1) += gpair.head(PR.ck - 1);
    GRADIENT.segment(PR.tl, PR.cl - 1) += gpair.segment(PR.ck - 1, PR.cl - 1);
    GRADIENT.tail(S.m) += gpair(PR.mt) * PR.grho;
  }
}
