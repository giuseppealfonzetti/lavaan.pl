#include "pairs.h"

pairs::Data::Data(const Rcpp::List D):
  freq(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(D["FREQ"])), ufreq(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(D["UFREQ"])),
  stratum(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(D["S"])), w(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(D["W"])),
  x(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(D["X"])), y(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(D["Y"])),
  tables(Rcpp::as<bool>(D["TABLES"])){}

pairs::Grid::Grid(const int CK, const int CL, const double RHO, const int ORDER):
  ck(CK), cl(CL), order(ORDER), rho(RHO), r2(1 - RHO * RHO), sr(std::sqrt(1 - RHO * RHO)){
  const double inf = std::numeric_limits<double>::infinity();
  a.resize(ck + 1);
  b.resize(cl + 1);
  a(0) = b(0) = -inf;
  a(ck) = b(cl) = inf;
  F = Eigen::MatrixXd::Zero(ck + 1, cl + 1);
  if(order == 0) return;
  Fa = Fb = Fr = Eigen::MatrixXd::Zero(ck + 1, cl + 1);
  if(order == 2) Faa = Fbb = Fab = Far = Fbr = Frr = Eigen::MatrixXd::Zero(ck + 1, cl + 1);
}

pairs::Grid::Grid(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const int CK, const int CL,
                  const double RHO, const int ORDER, const double SHK, const double SHL): Grid(CK, CL, RHO, ORDER){
  fill(TAU, TK, TL, SHK, SHL, 1, ck, 1, cl);
}

void pairs::Grid::fill(const Eigen::Ref<const Eigen::VectorXd> TAU, const int TK, const int TL, const double SHK, const double SHL,
                       const int S0, const int S1, const int T0, const int T1){
  const double inf = std::numeric_limits<double>::infinity();
  a.segment(1, ck - 1) = TAU.segment(TK, ck - 1).array() - SHK;
  b.segment(1, cl - 1) = TAU.segment(TL, cl - 1).array() - SHL;
  for(int s = S0; s <= S1; s++){
    for(int t = T0; t <= T1; t++){
      if(s == ck && t == cl) F(s, t) = 1;
      else if(s == ck) F(s, t) = R::pnorm(b(t), 0, 1, 1, 0);
      else if(t == cl) F(s, t) = R::pnorm(a(s), 0, 1, 1, 0);
      else F(s, t) = binorm::pbvnorm(a(s), b(t), rho);
    }
  }
  if(order == 0) return;

  for(int s = S0; s <= S1; s++){
    for(int t = T0; t <= T1; t++){
      const bool fa = s < ck, fb = t < cl;
      if(fa){
        const double pa = R::dnorm(a(s), 0, 1, 0);
        const double w = fb ? (b(t) - rho * a(s)) / sr : inf;
        Fa(s, t) = pa * R::pnorm(w, 0, 1, 1, 0);
        if(order == 2) Faa(s, t) = -a(s) * Fa(s, t) - (fb ? rho / sr * pa * R::dnorm(w, 0, 1, 0) : 0);
      }
      if(fb){
        const double pb = R::dnorm(b(t), 0, 1, 0);
        const double w = fa ? (a(s) - rho * b(t)) / sr : inf;
        Fb(s, t) = pb * R::pnorm(w, 0, 1, 1, 0);
        if(order == 2) Fbb(s, t) = -b(t) * Fb(s, t) - (fa ? rho / sr * pb * R::dnorm(w, 0, 1, 0) : 0);
      }
      if(fa && fb){
        const double f2 = binorm::dbvnorm(a(s), b(t), rho, 0);
        Fr(s, t) = f2;
        if(order == 2){
          Fab(s, t) = f2;
          Far(s, t) = -f2 * (a(s) - rho * b(t)) / r2;
          Fbr(s, t) = -f2 * (b(t) - rho * a(s)) / r2;
          Frr(s, t) = f2 * (rho / r2 + a(s) * b(t) / r2 -
                            rho * (a(s) * a(s) - 2 * rho * a(s) * b(t) + b(t) * b(t)) / (r2 * r2));
        }
      }
    }
  }
}

double pairs::Grid::cell(const int SK, const int SL, Eigen::VectorXd &G1, Eigen::MatrixXd &G2) const{
  const double pe = F(SK + 1, SL + 1) - F(SK + 1, SL) - F(SK, SL + 1) + F(SK, SL) + 1e-8;
  if(order == 0) return pe;

  const int mt = ck + cl - 2;
  G1 = Eigen::VectorXd::Zero(mt + 1);
  if(order == 2) G2 = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
  for(int e = 0; e < 2; e++){
    for(int f = 0; f < 2; f++){
      const int s = SK + 1 - e, t = SL + 1 - f;
      const double sg = (e == f) ? 1 : -1;
      const int px = (s > 0 && s < ck) ? s - 1 : -1;
      const int py = (t > 0 && t < cl) ? ck - 2 + t : -1;
      if(px >= 0) G1(px) += sg * Fa(s, t);
      if(py >= 0) G1(py) += sg * Fb(s, t);
      G1(mt) += sg * Fr(s, t);
      if(order < 2) continue;
      if(px >= 0){
        G2(px, px) += sg * Faa(s, t);
        G2(px, mt) += sg * Far(s, t);
        G2(mt, px) += sg * Far(s, t);
      }
      if(py >= 0){
        G2(py, py) += sg * Fbb(s, t);
        G2(py, mt) += sg * Fbr(s, t);
        G2(mt, py) += sg * Fbr(s, t);
      }
      if(px >= 0 && py >= 0){
        G2(px, py) += sg * Fab(s, t);
        G2(py, px) += sg * Fab(s, t);
      }
      G2(mt, mt) += sg * Frr(s, t);
    }
  }
  return pe;
}

static Eigen::VectorXd cov_grad(const model::Setup &ST, const int LO, const Eigen::MatrixXd &DA, const Eigen::MatrixXd &DB,
                                const Eigen::VectorXd &XA, const Eigen::VectorXd &XB, const Eigen::VectorXd &UA,
                                const Eigen::VectorXd &UB){
  const model::Spec &S = ST.S;
  Eigen::VectorXd g;
  if(S.direct){
    g = DA.transpose() * (ST.Sig * XB) + DB.transpose() * (ST.Sig * XA);
    for(int a = 0; a < S.ncorr; a++) g(LO + a) = XA(S.ra[a]) * XB(S.ca[a]) + XA(S.ca[a]) * XB(S.ra[a]);
    for(int v = 0; v < S.nvar; v++){
      const int h = S.sfree[v];
      g(LO + S.ncorr + v) = 2 * ST.sd(h) * ST.sd(h) * XA(h) * XB(h);
    }
  }else{
    g = DA.transpose() * ST.sd.cwiseProduct(ST.U.transpose() * UB) + DB.transpose() * ST.sd.cwiseProduct(ST.U.transpose() * UA);
    for(int a = 0; a < S.ncorr; a++) g(LO + a) = XA(S.ca[a]) * ST.alpha.col(a).dot(UB) + XB(S.ca[a]) * ST.alpha.col(a).dot(UA);
    for(int v = 0; v < S.nvar; v++){
      const int h = S.sfree[v];
      g(LO + S.ncorr + v) = XA(h) * ST.U.col(h).dot(UB) + XB(h) * ST.U.col(h).dot(UA);
    }
  }
  return g;
}

void pairs::Pair::size_reduced(const int NX){
  nx = NX;
  mt = ck + cl - 2;
  numk = ck == 1;
  numl = cl == 1;
  lk = numk ? 1 : ck - 1;
  ll = numl ? 1 : cl - 1;
  ml = lk + ll;
  nsk = nx + numk;
  nsl = nx + numl;
  m = ml + 1 + nsk + nsl;
}

pairs::Pair::Pair(const model::Layout &L, const int R, const double RHO, const int NX):
  r(R), k(L.pk[R]), l(L.pl[R]), ck(L.cat(k)), cl(L.cat(l)), tk(L.tk(k)), tl(L.tk(l)), off(L.off[R]), s(0), plain(false), rho(RHO){
  size_reduced(NX);
}

pairs::Pair::Pair(const model::Setup &ST, const int R, const bool GRAD, const int K): r(R){
  const model::Spec &S = ST.S;
  plain = S.plain;
  k = K >= 0 ? K : S.L.pk[r];
  l = K >= 0 ? K : S.L.pl[r];
  ck = S.L.cat(k);
  cl = S.L.cat(l);
  tk = S.L.tk(k);
  tl = S.L.tk(l);
  size_reduced(S.nx);
  off = K >= 0 ? 0 : S.L.off[r];
  const int ir = K >= 0 ? -1 : S.resid(k, l);

  if(S.direct){
    x = ST.Lt.row(k).transpose();
    y = ST.Lt.row(l).transpose();
    rho = x.dot(ST.Sig * y) + ST.rc(k, l);
  }else{
    x = ST.sd.cwiseProduct(ST.Lt.row(k).transpose());
    y = ST.sd.cwiseProduct(ST.Lt.row(l).transpose());
    X = ST.U * x;
    Y = ST.U * y;
    rho = X.dot(Y) + ST.rc(k, l);
  }
  c = rho;
  if(!plain) rho *= ST.s(k) * ST.s(l);
  if(!GRAD) return;

  nk = S.lidx[k].size();
  const int nl = S.lidx[l].size();
  lo = nk + nl;
  mb = lo + S.ncorr + S.nvar;
  s = s0 = mb + S.nreg + (ir >= 0);
  idx = S.lidx[k];
  idx.insert(idx.end(), S.lidx[l].begin(), S.lidx[l].end());
  for(int a = S.oc; a < S.ob + S.nreg; a++) idx.push_back(a);
  if(ir >= 0) idx.push_back(ir);

  Dk = Eigen::MatrixXd::Zero(S.q, s);
  Dl = Eigen::MatrixXd::Zero(S.q, s);
  Dk.leftCols(nk) = ST.J[k].leftCols(nk);
  Dl.middleCols(nk, nl) = ST.J[l].leftCols(nl);
  Dk.middleCols(mb, S.nreg) = ST.J[k].rightCols(S.nreg);
  Dl.middleCols(mb, S.nreg) = ST.J[l].rightCols(S.nreg);
  grho = cov_grad(ST, lo, Dk, Dl, x, y, X, Y);
  if(ir >= 0) grho(s - 1) = 1;
  if(plain) return;

  auto add = [&](const int i){
    if(i < 0) return -1;
    idx.push_back(i);
    return s++;
  };
  const int pnk = add(S.nidx[k]), pnl = add(S.nidx[l]);
  pa = s;
  for(int a = 0; a < S.nmean; a++) add(S.oa + a);
  psk = add(S.scidx[k]);
  psl = add(S.scidx[l]);
  pg = s;
  for(int a = 0; a < S.ngam; a++) add(S.og + a);
  pkk = s;
  for(const int i : S.kidx[k]) add(i);
  pkl = s;
  for(const int i : S.kidx[l]) add(i);

  auto pad = [&](const Eigen::VectorXd &V){
    Eigen::VectorXd o = Eigen::VectorXd::Zero(s);
    o.head(V.size()) = V;
    return o;
  };
  gc = pad(grho);
  gmk = pad(Dk.transpose() * ST.kappa);
  gml = pad(Dl.transpose() * ST.kappa);
  if(pnk >= 0) gmk(pnk) = 1;
  if(pnl >= 0) gml(pnl) = 1;
  for(int a = 0; a < S.nmean; a++){
    gmk(pa + a) = ST.Lt(k, S.mfree[a]);
    gml(pa + a) = ST.Lt(l, S.mfree[a]);
  }
  gsk = gsl = Eigen::VectorXd::Zero(s);
  if(S.thitem[k]){
    gvk = pad(cov_grad(ST, lo, Dk, Dk, x, x, X, X));
    gsk = gvk;
    if(psk >= 0) gsk(psk) += 2 * ST.th(k);
    gsk *= -0.5 * std::pow(ST.s(k), 3);
  }else if(psk >= 0){
    gsk(psk) = ST.s(k);
  }
  if(S.thitem[l]){
    gvl = pad(cov_grad(ST, lo, Dl, Dl, y, y, Y, Y));
    gsl = gvl;
    if(psl >= 0) gsl(psl) += 2 * ST.th(l);
    gsl *= -0.5 * std::pow(ST.s(l), 3);
  }else if(psl >= 0){
    gsl(psl) = ST.s(l);
  }

  Gk = Gl = Eigen::MatrixXd::Zero(s, nx);
  Gk.topRows(s0) = Dk.transpose() * ST.Gam;
  Gl.topRows(s0) = Dl.transpose() * ST.Gam;
  for(int a = 0; a < S.ngam; a++){
    Gk(pg + a, S.gcov[a]) = ST.Lt(k, S.gfac[a]);
    Gl(pg + a, S.gcov[a]) = ST.Lt(l, S.gfac[a]);
  }
  for(std::size_t a = 0; a < S.kidx[k].size(); a++) Gk(pkk + a, S.kcov[k][a]) = 1;
  for(std::size_t a = 0; a < S.kidx[l].size(); a++) Gl(pkl + a, S.kcov[l][a]) = 1;

  const double sk = ST.s(k), sl = ST.s(l);
  Ju = Eigen::MatrixXd::Zero(m, mt + s);
  for(int a = 0; a < ck - 1; a++){
    Ju(a, a) = sk;
    Ju.row(a).tail(s) = (-sk * gmk + (ST.tau(tk + a) - ST.mu(k)) * gsk).transpose();
  }
  if(numk) Ju.row(0).tail(s) = (-sk * gmk - ST.mu(k) * gsk).transpose();
  for(int b = 0; b < cl - 1; b++){
    Ju(lk + b, ck - 1 + b) = sl;
    Ju.row(lk + b).tail(s) = (-sl * gml + (ST.tau(tl + b) - ST.mu(l)) * gsl).transpose();
  }
  if(numl) Ju.row(lk).tail(s) = (-sl * gml - ST.mu(l) * gsl).transpose();
  Ju.row(ml).tail(s) = (sl * c * gsk + sk * c * gsl + sk * sl * gc).transpose();
  const int ok = ml + 1, ol = ok + nsk;
  for(int a = 0; a < nx; a++){
    Ju.row(ok + a).tail(s) = (sk * Gk.col(a) + ST.Pi(k, a) * gsk).transpose();
    Ju.row(ol + a).tail(s) = (sl * Gl.col(a) + ST.Pi(l, a) * gsl).transpose();
  }
  if(numk) Ju.row(ok + nx).tail(s) = gsk.transpose();
  if(numl) Ju.row(ol + nx).tail(s) = gsl.transpose();
}

std::vector<int> pairs::Pair::local() const{
  std::vector<int> I;
  for(int a = 0; a < ck - 1; a++) I.push_back(tk + a);
  for(int a = 0; a < cl - 1; a++) I.push_back(tl + a);
  I.insert(I.end(), idx.begin(), idx.end());
  return I;
}

Eigen::VectorXd pairs::Pair::lift(const Eigen::Ref<const Eigen::VectorXd> V) const{
  if(!plain) return Ju.transpose() * V;
  Eigen::VectorXd o(mt + s);
  o.head(mt) = V.head(mt);
  o.tail(s) = V(mt) * grho;
  return o;
}

void pairs::Pair::slots(const Eigen::Ref<const Eigen::VectorXd> X, const double YK, const double YL, Eigen::VectorXd &WK,
                        Eigen::VectorXd &WL) const{
  WK.resize(nsk);
  WL.resize(nsl);
  WK.head(nx) = -X;
  WL.head(nx) = -X;
  if(numk) WK(nx) = YK;
  if(numl) WL(nx) = YL;
}

void pairs::Pair::add_score(const double N, const Eigen::Ref<const Eigen::VectorXd> G, const Eigen::Ref<const Eigen::VectorXd> WK,
                            const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::VectorXd &OUT) const{
  const int L = ml + 1;
  OUT.head(L) += N * G;
  if(nsk > 0) OUT.segment(L, nsk) += N * G.head(lk).sum() * WK;
  if(nsl > 0) OUT.segment(L + nsk, nsl) += N * G.segment(lk, ll).sum() * WL;
}

void pairs::Pair::add_info(const double N, const Eigen::Ref<const Eigen::MatrixXd> M, const Eigen::Ref<const Eigen::VectorXd> WK,
                           const Eigen::Ref<const Eigen::VectorXd> WL, Eigen::MatrixXd &OUT) const{
  const int L = ml + 1;
  OUT.topLeftCorner(L, L) += N * M;
  if(nsk + nsl == 0) return;
  const Eigen::VectorXd u = M.leftCols(lk).rowwise().sum(), v = M.middleCols(lk, ll).rowwise().sum();
  const double uu = u.head(lk).sum(), uv = u.segment(lk, ll).sum(), vv = v.segment(lk, ll).sum();
  OUT.block(0, L, L, nsk) += N * u * WK.transpose();
  OUT.block(L, 0, nsk, L) += N * WK * u.transpose();
  OUT.block(0, L + nsk, L, nsl) += N * v * WL.transpose();
  OUT.block(L + nsk, 0, nsl, L) += N * WL * v.transpose();
  OUT.block(L, L, nsk, nsk) += N * uu * WK * WK.transpose();
  OUT.block(L, L + nsk, nsk, nsl) += N * uv * WK * WL.transpose();
  OUT.block(L + nsk, L, nsl, nsk) += N * uv * WL * WK.transpose();
  OUT.block(L + nsk, L + nsk, nsl, nsl) += N * vv * WL * WL.transpose();
}

void pairs::Pair::add_scale(const double N, const Eigen::Ref<const Eigen::VectorXd> SCALE, Eigen::VectorXd &GU, Eigen::MatrixXd *M) const{
  const int ok = ml + 1 + nx, ol = ml + 1 + nsk + nx;
  if(numk){
    if(GU.size() > 0) GU(ok) += N / SCALE(k);
    if(M) (*M)(ok, ok) += N / (SCALE(k) * SCALE(k));
  }
  if(numl){
    if(GU.size() > 0) GU(ol) += N / SCALE(l);
    if(M) (*M)(ol, ol) += N / (SCALE(l) * SCALE(l));
  }
}

static void regression_curvature(const model::Setup &ST, const pairs::Pair &PR, const int K,
                                 const Eigen::MatrixXd &D, const Eigen::VectorXd &W, Eigen::MatrixXd &HRHO){
  const model::Spec &S = ST.S;
  const Eigen::VectorXd Aw = ST.A * W;
  for(int b = 0; b < S.nreg; b++){
    const int i = S.bi[b], j = S.bj[b];
    const Eigen::VectorXd col = Aw(j) * D.row(i).head(PR.lo).transpose();
    HRHO.block(0, PR.mb + b, PR.lo, 1) += col;
    HRHO.block(PR.mb + b, 0, 1, PR.lo) += col.transpose();
    for(int c = 0; c < S.nreg; c++){
      const int r = S.bi[c], t = S.bj[c];
      HRHO(PR.mb + b, PR.mb + c) += ST.Lt(K, r) * ST.A(t, i) * Aw(j) + ST.Lt(K, i) * ST.A(j, r) * Aw(t);
    }
  }
}

static Eigen::MatrixXd cov_hessian(const model::Setup &ST, const pairs::Pair &PR, const Eigen::MatrixXd &Dk,
                                   const Eigen::MatrixXd &Dl, const Eigen::VectorXd &x, const Eigen::VectorXd &y,
                                   const Eigen::VectorXd &X, const Eigen::VectorXd &Y,
                                   const Eigen::Ref<const Eigen::VectorXd> grho, const int K, const int L){
  const model::Spec &S = ST.S;
  const int q = S.q, lo = PR.lo, ncorr = S.ncorr, nvar = S.nvar;

  if(S.direct){
    const Eigen::MatrixXd ll = Dk.transpose() * ST.Sig * Dl;
    Eigen::MatrixXd hrho = ll + ll.transpose();
    for(int a = 0; a < ncorr; a++){
      Eigen::VectorXd ek = Eigen::VectorXd::Zero(q), el = Eigen::VectorXd::Zero(q);
      ek(S.ra[a]) += y(S.ca[a]); ek(S.ca[a]) += y(S.ra[a]);
      el(S.ra[a]) += x(S.ca[a]); el(S.ca[a]) += x(S.ra[a]);
      const Eigen::VectorXd col = Dk.transpose() * ek + Dl.transpose() * el;
      hrho.col(lo + a) += col;
      hrho.row(lo + a) += col.transpose();
    }
    for(int v = 0; v < nvar; v++){
      const int h = S.sfree[v], js = lo + ncorr + v;
      const double w = 2 * ST.sd(h) * ST.sd(h);
      const Eigen::VectorXd col = w * (Dk.row(h).transpose() * y(h) + Dl.row(h).transpose() * x(h));
      hrho.col(js) += col;
      hrho.row(js) += col.transpose();
      hrho(js, js) += 2 * w * x(h) * y(h);
    }
    if(S.nreg > 0){
      regression_curvature(ST, PR, K, Dk, ST.Sig * y, hrho);
      regression_curvature(ST, PR, L, Dl, ST.Sig * x, hrho);
    }
    return hrho;
  }

  const Eigen::MatrixXd DJk = ST.sd.asDiagonal() * Dk, DJl = ST.sd.asDiagonal() * Dl;
  Eigen::MatrixXd dX = ST.U * DJk, dY = ST.U * DJl;
  for(int a = 0; a < ncorr; a++){
    dX.col(lo + a) = x(S.ca[a]) * ST.alpha.col(a);
    dY.col(lo + a) = y(S.ca[a]) * ST.alpha.col(a);
  }
  for(int v = 0; v < nvar; v++){
    dX.col(lo + ncorr + v) = x(S.sfree[v]) * ST.U.col(S.sfree[v]);
    dY.col(lo + ncorr + v) = y(S.sfree[v]) * ST.U.col(S.sfree[v]);
  }
  Eigen::MatrixXd hrho = dX.transpose() * dY;
  hrho += hrho.transpose().eval();

  const Eigen::VectorXd aY = ST.alpha.transpose() * Y, aX = ST.alpha.transpose() * X;
  const Eigen::VectorXd uY = ST.U.transpose() * Y, uX = ST.U.transpose() * X;
  for(int a = 0; a < ncorr; a++){
    const Eigen::VectorXd v = DJk.row(S.ca[a]).transpose() * aY(a) + DJl.row(S.ca[a]).transpose() * aX(a);
    hrho.col(lo + a) += v;
    hrho.row(lo + a) += v.transpose();
  }
  for(int s = 0; s < nvar; s++){
    const int h = S.sfree[s];
    const Eigen::VectorXd v = DJk.row(h).transpose() * uY(h) + DJl.row(h).transpose() * uX(h);
    hrho.col(lo + ncorr + s) += v;
    hrho.row(lo + ncorr + s) += v.transpose();
  }
  for(int a = 0; a < ncorr; a++){
    for(int b = a; b < ncorr && S.ca[b] == S.ca[a]; b++){
      const int c = S.ca[a], i0 = S.ra[a], i1 = S.ra[b];
      Eigen::VectorXd dd = Eigen::VectorXd::Zero(q);
      if(i0 < i1){
        dd(i1) = -ST.P(i1, c) * ST.Z(i0, c) * ST.C(i1, c) * ST.C(i1, c);
        for(int i = i1 + 1; i <= c; i++) dd(i) = ST.U(i, c) * ST.Z(i0, c) * ST.Z(i1, c);
      }else{
        dd(i0) = -2 * ST.P(i0, c) * ST.Z(i0, c) * ST.C(i0, c) * ST.C(i0, c);
        for(int i = i0 + 1; i <= c; i++) dd(i) = ST.U(i, c) * (2 * ST.Z(i0, c) * ST.Z(i0, c) - 1);
      }
      const double v = x(c) * dd.dot(Y) + y(c) * dd.dot(X);
      hrho(lo + a, lo + b) += v;
      if(b != a) hrho(lo + b, lo + a) += v;
    }
  }
  for(int s = 0; s < nvar; s++){
    const int js = lo + ncorr + s;
    hrho(js, js) += grho(js);
    for(int a = 0; a < ncorr; a++){
      if(S.ca[a] != S.sfree[s]) continue;
      hrho(lo + a, js) += grho(lo + a);
      hrho(js, lo + a) += grho(lo + a);
    }
  }
  if(S.nreg > 0){
    regression_curvature(ST, PR, K, Dk, ST.sd.cwiseProduct(uY), hrho);
    regression_curvature(ST, PR, L, Dl, ST.sd.cwiseProduct(uX), hrho);
  }
  return hrho;
}

Eigen::MatrixXd pairs::rho_hessian(const model::Setup &ST, const Pair &PR){
  return cov_hessian(ST, PR, PR.Dk, PR.Dl, PR.x, PR.y, PR.X, PR.Y, PR.grho, PR.k, PR.l);
}

Eigen::MatrixXd pairs::curvature(const model::Setup &ST, const Pair &PR, const Eigen::Ref<const Eigen::VectorXd> GU){
  const model::Spec &S = ST.S;
  const int mt = PR.mt, s = PR.s, s0 = PR.s0, ck = PR.ck, cl = PR.cl;
  const double gr = GU(PR.ml);
  auto hmu = [&](const int K, const Eigen::MatrixXd &D, const double A, const Eigen::VectorXd &WX){
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(s, s);
    Eigen::MatrixXd h0 = Eigen::MatrixXd::Zero(s0, s0);
    const Eigen::VectorXd w = A * ST.kappa + ST.Gam * WX;
    regression_curvature(ST, PR, K, D, w, h0);
    H.topLeftCorner(s0, s0) = h0;
    for(int a = 0; a < S.nmean; a++){
      const Eigen::VectorXd col = A * D.row(S.mfree[a]).transpose();
      H.block(0, PR.pa + a, s0, 1) += col;
      H.block(PR.pa + a, 0, 1, s0) += col.transpose();
    }
    for(int a = 0; a < S.ngam; a++){
      const Eigen::VectorXd col = WX(S.gcov[a]) * D.row(S.gfac[a]).transpose();
      H.block(0, PR.pg + a, s0, 1) += col;
      H.block(PR.pg + a, 0, 1, s0) += col.transpose();
    }
    return H;
  };
  auto hsc = [&](const int K, const Eigen::MatrixXd &D, const Eigen::VectorXd &x, const Eigen::VectorXd &X,
                 const Eigen::VectorXd &gv, const int ps){
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(s, s);
    const double sc = ST.s(K);
    if(!S.thitem[K]){
      if(ps >= 0) H(ps, ps) = sc;
      return H;
    }
    Eigen::VectorXd gt = gv;
    H.topLeftCorner(s0, s0) = cov_hessian(ST, PR, D, D, x, x, X, X, gv.head(s0), K, K);
    if(ps >= 0){
      gt(ps) += 2 * ST.th(K);
      H(ps, ps) += 4 * ST.th(K);
    }
    return Eigen::MatrixXd(0.75 * std::pow(sc, 5) * gt * gt.transpose() - 0.5 * std::pow(sc, 3) * H);
  };

  const Eigen::MatrixXd Hsk = hsc(PR.k, PR.Dk, PR.x, PR.X, PR.gvk, PR.psk);
  const Eigen::MatrixXd Hsl = hsc(PR.l, PR.Dl, PR.y, PR.Y, PR.gvl, PR.psl);
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(mt + s, mt + s);
  auto side = [&](const int K, const int TK, const int CK, const int LOC, const int LN, const int HEAD, const int OX, const bool NUM,
                  const Eigen::MatrixXd &D, const Eigen::VectorXd &gm, const Eigen::MatrixXd &Gx, const Eigen::VectorXd &gs,
                  const Eigen::MatrixXd &Hs){
    const Eigen::VectorXd gt = GU.segment(LOC, LN), wx = GU.segment(OX, S.nx);
    const double a = gt.sum();
    double b = ST.Pi.row(K).dot(wx);
    if(NUM){
      b += -ST.mu(K) * gt(0) + GU(OX + S.nx);
    }else{
      b += gt.dot((ST.tau.segment(TK, CK - 1).array() - ST.mu(K)).matrix());
    }
    Eigen::VectorXd e = Eigen::VectorXd::Zero(mt + s), g = Eigen::VectorXd::Zero(mt + s);
    if(!NUM) e.segment(HEAD, CK - 1) = gt;
    e.tail(s) = -a * gm + Gx * wx;
    g.tail(s) = gs;
    W += g * e.transpose() + e * g.transpose();
    W.bottomRightCorner(s, s) += b * Hs + ST.s(K) * hmu(K, D, -a, wx);
  };
  side(PR.k, PR.tk, ck, 0, PR.lk, 0, PR.ml + 1, PR.numk, PR.Dk, PR.gmk, PR.Gk, PR.gsk, Hsk);
  side(PR.l, PR.tl, cl, PR.lk, PR.ll, ck - 1, PR.ml + 1 + PR.nsk, PR.numl, PR.Dl, PR.gml, PR.Gl, PR.gsl, Hsl);
  if(gr != 0){
    const double sk = ST.s(PR.k), sl = ST.s(PR.l), c = PR.c;
    const Eigen::VectorXd &gk = PR.gsk, &gl = PR.gsl, &gc = PR.gc;
    Eigen::MatrixXd R = c * (gk * gl.transpose()) + sl * (gk * gc.transpose()) + sk * (gl * gc.transpose());
    R += R.transpose().eval();
    R += c * sl * Hsk + c * sk * Hsl;
    R.topLeftCorner(s0, s0) += sk * sl * rho_hessian(ST, PR);
    W.bottomRightCorner(s, s) += gr * R;
  }
  return W;
}

double pairs::unit(const Reduced &R, const Data &D, const Pair &PR, Grid &G, const int U, const int ORDER, Work &W){
  const int k = PR.k, l = PR.l, t = D.stratum(U);
  const double yk = D.y(U, k), yl = D.y(U, l);
  const double shk = D.x.col(t).dot(R.pis.row(k).transpose()), shl = D.x.col(t).dot(R.pis.row(l).transpose());
  if(!PR.numk && !PR.numl){
    const int sk = yk, sl = yl;
    G.fill(R.tau, PR.tk, PR.tl, shk, shl, std::max(sk, 1), sk + 1, std::max(sl, 1), sl + 1);
    const double pe = G.cell(sk, sl, W.gc, W.G2);
    if(ORDER > 0) W.g = W.gc / pe;
    if(ORDER == 2){
      W.M.noalias() = W.gc * W.gc.transpose() / (pe * pe);
      W.M -= W.G2 / pe;
    }
    return std::log(pe);
  }
  const double rho = PR.rho, r2 = 1 - rho * rho, r = std::sqrt(r2), log2pi = std::log(2 * M_PI);
  if(PR.numk && PR.numl){
    const double zk = R.scale(k) * yk + R.loc(k) - shk, zl = R.scale(l) * yl + R.loc(l) - shl;
    const double q = zk * zk - 2 * rho * zk * zl + zl * zl, r4 = r2 * r2;
    if(ORDER > 0){
      W.g.resize(3);
      W.g << -(zk - rho * zl) / r2, -(zl - rho * zk) / r2, (rho + zk * zl) / r2 - rho * q / r4;
    }
    if(ORDER == 2){
      W.M.resize(3, 3);
      W.M(0, 0) = W.M(1, 1) = 1 / r2;
      W.M(0, 1) = W.M(1, 0) = -rho / r2;
      W.M(0, 2) = W.M(2, 0) = -zl / r2 + 2 * rho * (zk - rho * zl) / r4;
      W.M(1, 2) = W.M(2, 1) = -zk / r2 + 2 * rho * (zl - rho * zk) / r4;
      W.M(2, 2) = -((1 + rho * rho) / r4 + 4 * rho * zk * zl / r4 - q / r4 - 4 * rho * rho * q / (r4 * r2));
    }
    return -log2pi - 0.5 * std::log(r2) - 0.5 * q / r2 + std::log(R.scale(k)) + std::log(R.scale(l));
  }
  const int num = PR.numk ? k : l;
  const int cat = PR.numk ? PR.cl : PR.ck, to = PR.numk ? PR.tl : PR.tk, s = PR.numk ? yl : yk;
  const double z = R.scale(num) * (PR.numk ? yk : yl) + R.loc(num) - (PR.numk ? shk : shl), sho = PR.numk ? shl : shk;
  const bool top = s == cat - 1, bot = s == 0;
  const double hi = top ? 0 : R.tau(to + s) - sho, lo = bot ? 0 : R.tau(to + s - 1) - sho;
  const double uh = (hi - rho * z) / r, ul = (lo - rho * z) / r;
  const double ph = top ? 0 : R::dnorm(uh, 0, 1, 0), pl = bot ? 0 : R::dnorm(ul, 0, 1, 0);
  const double pr = (top ? 1 : R::pnorm(uh, 0, 1, 1, 0)) - (bot ? 0 : R::pnorm(ul, 0, 1, 1, 0)) + 1e-8;
  const double r3 = r2 * r;
  Eigen::Vector4d du, dl;
  du << -rho / r, 0, 1 / r, (rho * hi - z) / r3;
  dl << -rho / r, 1 / r, 0, (rho * lo - z) / r3;
  const Eigen::Vector4d dd = ph * du - pl * dl;
  const int ob = PR.numk ? 1 : 0, iz = PR.numk ? 0 : PR.lk;
  const int idx[4] = {iz, bot ? -1 : ob + s - 1, top ? -1 : ob + s, PR.ml};
  if(ORDER > 0){
    W.g = Eigen::VectorXd::Zero(PR.ml + 1);
    for(int a = 0; a < 4; a++) if(idx[a] >= 0) W.g(idx[a]) += dd(a) / pr;
    W.g(iz) -= z;
  }
  if(ORDER == 2){
    Eigen::Matrix4d H = Eigen::Matrix4d::Zero();
    if(!top){
      H += -ph * uh * du * du.transpose();
      H(0, 3) += -ph / r3;
      H(3, 0) += -ph / r3;
      H(2, 3) += ph * rho / r3;
      H(3, 2) += ph * rho / r3;
      H(3, 3) += ph * (hi / r3 + 3 * rho * (rho * hi - z) / (r3 * r2));
    }
    if(!bot){
      H += pl * ul * dl * dl.transpose();
      H(0, 3) += pl / r3;
      H(3, 0) += pl / r3;
      H(1, 3) += -pl * rho / r3;
      H(3, 1) += -pl * rho / r3;
      H(3, 3) += -pl * (lo / r3 + 3 * rho * (rho * lo - z) / (r3 * r2));
    }
    const Eigen::Matrix4d M4 = dd * dd.transpose() / (pr * pr) - H / pr;
    W.M = Eigen::MatrixXd::Zero(PR.ml + 1, PR.ml + 1);
    for(int a = 0; a < 4; a++){
      for(int b = 0; b < 4; b++) if(idx[a] >= 0 && idx[b] >= 0) W.M(idx[a], idx[b]) += M4(a, b);
    }
    W.M(iz, iz) += 1;
  }
  return -0.5 * log2pi - 0.5 * z * z + std::log(R.scale(num)) + std::log(pr);
}

void pairs::moments(const model::Layout &LAY, const Reduced &R, const Data &D, const Pair &PR, const int ORDER, double &LL,
                    Eigen::VectorXd &GU, Eigen::MatrixXd &M, Eigen::MatrixXd *SC){
  const int L = PR.ml + 1, nc = PR.ck * PR.cl;
  Work W;
  Grid G(PR.ck, PR.cl, PR.rho, ORDER);
  if(D.tables && !PR.numk && !PR.numl){
    Eigen::VectorXd g, gt(ORDER > 0 ? L : 0), gc(SC ? PR.m : 0);
    Eigen::MatrixXd G2, Mt(ORDER == 2 ? L : 0, ORDER == 2 ? L : 0);
    for(int t = 0; t < D.x.cols(); t++){
      G.fill(R.tau, PR.tk, PR.tl, D.x.col(t).dot(R.pis.row(PR.k).transpose()), D.x.col(t).dot(R.pis.row(PR.l).transpose()),
             1, PR.ck, 1, PR.cl);
      PR.slots(D.x.col(t), 0, 0, W.wk, W.wl);
      const int o = t * LAY.ncell + PR.off;
      if(ORDER > 0) gt.setZero();
      if(ORDER == 2) Mt.setZero();
      for(int sk = 0; sk < PR.ck; sk++){
        for(int sl = 0; sl < PR.cl; sl++){
          const int c = sk * PR.cl + sl;
          const double n = D.freq(o + c);
          if(n == 0 && !SC) continue;
          const double pe = G.cell(sk, sl, g, G2);
          if(SC){
            const int col = t * nc + c;
            if(PR.plain){
              SC->col(col).head(PR.mt) = g.head(PR.mt) / pe;
              SC->col(col).tail(PR.s) = g(PR.mt) / pe * PR.grho;
            }else{
              gc.setZero();
              PR.add_score(1 / pe, g, W.wk, W.wl, gc);
              SC->col(col) = PR.Ju.transpose() * gc;
            }
          }
          if(n == 0) continue;
          LL += n * std::log(pe);
          if(ORDER > 0) gt += n / pe * g;
          if(ORDER == 2) Mt += n * (g * g.transpose() / (pe * pe) - G2 / pe);
        }
      }
      if(ORDER > 0) PR.add_score(1, gt, W.wk, W.wl, GU);
      if(ORDER == 2) PR.add_info(1, Mt, W.wk, W.wl, M);
    }
    return;
  }
  for(int u = 0; u < D.y.rows(); u++){
    const double yk = D.y(u, PR.k), yl = D.y(u, PR.l);
    if(std::isnan(yk) || std::isnan(yl)) continue;
    const double n = D.w(u);
    LL += n * unit(R, D, PR, G, u, ORDER, W);
    if(ORDER == 0) continue;
    PR.slots(D.x.col(D.stratum(u)), yk, yl, W.wk, W.wl);
    PR.add_score(n, W.g, W.wk, W.wl, GU);
    if(ORDER == 2) PR.add_info(n, W.M, W.wk, W.wl, M);
    PR.add_scale(n, R.scale, GU, ORDER == 2 ? &M : nullptr);
  }
}

void pairs::pair_contribution(const model::Setup &ST, const Data &D, const int R, const int ORDER, double &LL,
                              Eigen::VectorXd &GRADIENT){
  const Pair PR(ST, R, ORDER == 1);
  Eigen::VectorXd gpair = Eigen::VectorXd::Zero(ORDER == 1 ? PR.m : 0);
  Eigen::MatrixXd M;
  moments(ST.S.L, Reduced(ST), D, PR, ORDER, LL, gpair, M);
  if(ORDER != 1) return;
  if(PR.plain){
    GRADIENT.segment(PR.tk, PR.ck - 1) += gpair.head(PR.ck - 1);
    GRADIENT.segment(PR.tl, PR.cl - 1) += gpair.segment(PR.ck - 1, PR.cl - 1);
    for(int a = 0; a < PR.s; a++) GRADIENT(PR.idx[a]) += gpair(PR.mt) * PR.grho(a);
    return;
  }
  const std::vector<int> I = PR.local();
  const Eigen::VectorXd gl = PR.lift(gpair);
  for(std::size_t a = 0; a < I.size(); a++) GRADIENT(I[a]) += gl(a);
}

void pairs::univariate_contribution(const Eigen::Ref<const Eigen::VectorXd> TAU, const model::Layout &L,
                                    const Eigen::Ref<const Eigen::VectorXd> UFREQ, const int ORDER, double &LL,
                                    Eigen::VectorXd &GRADIENT, Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE,
                                    const Eigen::Ref<const Eigen::VectorXd> ITEMW){
  if(UFREQ.size() == 0) return;
  for(int k = 0, uk = 0; k < L.p; uk += L.cat(k), k++){
    const int ck = L.cat(k), tk = L.tk(k);
    const double wk = ITEMW.size() > 0 ? ITEMW(k) : 1;
    if(wk == 0) continue;
    for(int a = 0; a < ck; a++){
      const double nc = wk * UFREQ(uk + a);
      if(nc == 0) continue;
      const bool up = a < ck - 1, dw = a > 0;
      const double hi = up ? TAU(tk + a) : 0, lo = dw ? TAU(tk + a - 1) : 0;
      const double pe = (up ? R::pnorm(hi, 0, 1, 1, 0) : 1) - (dw ? R::pnorm(lo, 0, 1, 1, 0) : 0) + 1e-8;
      LL += nc * std::log(pe);
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

void pairs::univariate(const model::Setup &ST, const Data &D, const int ORDER, double &LL, Eigen::VectorXd &GRADIENT,
                       Eigen::MatrixXd &HESSIAN, Eigen::MatrixXd &USCORE, const Eigen::Ref<const Eigen::VectorXd> ITEMW){
  const model::Spec &S = ST.S;
  if(S.plain || D.ufreq.size() == 0){
    univariate_contribution(ST.taus, S.L, D.ufreq, ORDER, LL, GRADIENT, HESSIAN, USCORE, ITEMW);
    return;
  }
  const int nt = S.L.nthr, uc = S.L.cat.sum();
  std::vector<Pair> items;
  items.reserve(ORDER > 0 ? S.p : 0);
  for(int k = 0; k < S.p && ORDER > 0; k++) items.emplace_back(ST, -1, true, k);
  for(int t = 0; t < D.x.cols(); t++){
    const auto xt = D.x.col(t);
    const Eigen::VectorXd shift = ST.pis * xt;
    Eigen::VectorXd taut = ST.taus;
    for(int k = 0; k < S.p; k++) taut.segment(S.L.tk(k), S.L.cat(k) - 1).array() -= shift(k);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nt);
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(ORDER == 2 ? nt : 0, ORDER == 2 ? nt : 0);
    Eigen::MatrixXd us = Eigen::MatrixXd::Zero(USCORE.size() > 0 ? nt : 0, uc);
    univariate_contribution(taut, S.L, D.ufreq.segment(t * uc, uc), ORDER, LL, g, h, us, ITEMW);
    if(ORDER == 0) continue;
    for(int k = 0, uk = 0; k < S.p; uk += S.L.cat(k), k++){
      const Pair &PR = items[k];
      const int m = PR.ck - 1, tk = PR.tk;
      const std::vector<int> I = PR.local();
      Eigen::MatrixXd Jk = PR.Ju.topRows(m);
      if(S.nx > 0) Jk.rowwise() -= Eigen::RowVectorXd(xt.transpose() * PR.Ju.middleRows(PR.ml + 1, S.nx));
      const Eigen::VectorXd gl = Jk.transpose() * g.segment(tk, m);
      for(std::size_t a = 0; a < I.size(); a++) GRADIENT(I[a]) += gl(a);
      for(int c = 0; c < S.L.cat(k) && us.size() > 0; c++){
        const Eigen::VectorXd sl = Jk.transpose() * us.block(tk, uk + c, m, 1);
        for(std::size_t a = 0; a < I.size(); a++) USCORE(I[a], t * uc + uk + c) += sl(a);
      }
      if(ORDER < 2) continue;
      Eigen::VectorXd gu = Eigen::VectorXd::Zero(PR.m), gt = Eigen::VectorXd::Zero(PR.ml + 1), wk, wl;
      gt.head(m) = g.segment(tk, m);
      PR.slots(xt, 0, 0, wk, wl);
      PR.add_score(1, gt, wk, wl, gu);
      const Eigen::MatrixXd HL = Jk.transpose() * h.block(tk, tk, m, m) * Jk - curvature(ST, PR, gu);
      for(std::size_t a = 0; a < I.size(); a++){
        for(std::size_t b = 0; b < I.size(); b++) HESSIAN(I[a], I[b]) += HL(a, b);
      }
    }
  }
}
