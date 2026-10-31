#include "model.h"

model::Layout::Layout(const Eigen::Ref<const Eigen::VectorXd> C_VEC): p(C_VEC.size()), P(0), nthr(0), ncell(0){
  cat = C_VEC.cast<int>();
  tk.resize(p);
  for(int k = 0; k < p; k++){
    tk(k) = nthr;
    nthr += cat(k) - 1;
  }
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++){
      pk.push_back(k);
      pl.push_back(l);
      off.push_back(ncell);
      ncell += cat(k) * cat(l);
      P++;
    }
  }
}

model::Spec::Spec(const Rcpp::List MODEL, const Eigen::Ref<const Eigen::VectorXd> C_VEC): L(C_VEC){
  constrmat = Rcpp::as<Eigen::MatrixXd>(MODEL["CONSTRMAT"]);
  constrlogsd = Rcpp::as<Eigen::VectorXd>(MODEL["CONSTRLOGSD"]);
  regmat = Rcpp::as<Eigen::MatrixXd>(MODEL["REGMAT"]);
  rescov = Rcpp::as<Eigen::MatrixXd>(MODEL["RESCOV"]);
  constrint = Rcpp::as<Eigen::VectorXd>(MODEL["CONSTRINT"]);
  constrmean = Rcpp::as<Eigen::VectorXd>(MODEL["CONSTRMEAN"]);
  constrlogscale = Rcpp::as<Eigen::VectorXd>(MODEL["CONSTRLOGSCALE"]);
  thetapar = Rcpp::as<bool>(MODEL["THETAPAR"]);
  corrflag = Rcpp::as<int>(MODEL["CORRFLAG"]);
  direct = corrflag == 2;
  p = constrmat.rows();
  q = constrmat.cols();
  ol = L.nthr;

  lidx.assign(p, std::vector<int>());
  lfac.assign(p, std::vector<int>());
  nload = 0;
  for(int h = 0; h < q; h++){
    for(int j = 0; j < p; j++){
      if(!std::isfinite(constrmat(j, h))){
        lidx[j].push_back(ol + nload);
        lfac[j].push_back(h);
        nload++;
      }
    }
  }

  oc = ol + nload;
  for(int c = 1; c < q && corrflag > 0; c++){
    for(int r = 0; r < c; r++){
      ra.push_back(r);
      ca.push_back(c);
    }
  }
  ncorr = ra.size();

  ov = oc + ncorr;
  for(int h = 0; h < q; h++) if(!std::isfinite(constrlogsd(h))) sfree.push_back(h);
  nvar = sfree.size();

  ob = ov + nvar;
  for(int j = 0; j < regmat.cols(); j++){
    for(int i = 0; i < regmat.rows(); i++){
      if(!std::isfinite(regmat(i, j))){
        bi.push_back(i);
        bj.push_back(j);
      }
    }
  }
  nreg = bi.size();

  orc = ob + nreg;
  nres = 0;
  resid = Eigen::MatrixXi::Constant(p, p, -1);
  for(int c = 1; c < rescov.cols(); c++){
    for(int r = 0; r < c; r++){
      if(!std::isfinite(rescov(r, c))){
        resid(r, c) = resid(c, r) = orc + nres;
        nres++;
      }
    }
  }
  on = orc + nres;
  nidx.assign(p, -1);
  nint = 0;
  for(int j = 0; j < p; j++) if(!std::isfinite(constrint(j))) nidx[j] = on + nint++;
  oa = on + nint;
  for(int h = 0; h < q; h++) if(!std::isfinite(constrmean(h))) mfree.push_back(h);
  nmean = mfree.size();
  os = oa + nmean;
  scidx.assign(p, -1);
  nscale = 0;
  for(int j = 0; j < p; j++) if(!std::isfinite(constrlogscale(j))) scidx[j] = os + nscale++;
  d = os + nscale;
  plain = !thetapar && nint + nmean + nscale == 0 && (constrint.array() == 0).all() && (constrmean.array() == 0).all() &&
    (constrlogscale.array() == 0).all();
}

model::Setup::Setup(const Eigen::Ref<const Eigen::VectorXd> THETA, const Spec &SPEC): S(SPEC){
  const int p = S.p, q = S.q;
  if(THETA.size() != S.d) Rcpp::stop("THETA has length %d, the model has %d parameters.", (int)THETA.size(), S.d);

  tau = THETA.head(S.L.nthr);
  Lam = S.constrmat;
  for(int j = 0; j < p; j++){
    for(std::size_t a = 0; a < S.lidx[j].size(); a++) Lam(j, S.lfac[j][a]) = THETA(S.lidx[j][a]);
  }
  sd = S.constrlogsd;
  for(int s = 0; s < S.nvar; s++) sd(S.sfree[s]) = THETA(S.ov + s);
  sd = sd.array().exp();

  Z = Eigen::MatrixXd::Zero(q, q);
  if(!S.direct) for(int a = 0; a < S.ncorr; a++) Z(S.ra[a], S.ca[a]) = std::tanh(THETA(S.oc + a));
  C = (1 - Z.array().square()).sqrt().matrix();
  P = Eigen::MatrixXd::Ones(q, q);
  U = Eigen::MatrixXd::Zero(q, q);
  for(int c = 0; c < q; c++){
    for(int r = 1; r <= c; r++) P(r, c) = P(r - 1, c) * C(r - 1, c);
    for(int r = 0; r < c; r++) U(r, c) = Z(r, c) * P(r, c);
    U(c, c) = P(c, c);
  }
  alpha = Eigen::MatrixXd::Zero(q, S.ncorr);
  for(int a = 0; a < S.ncorr && !S.direct; a++){
    const int r = S.ra[a], c = S.ca[a];
    alpha(r, a) = P(r, c) * C(r, c) * C(r, c);
    for(int i = r + 1; i <= c; i++) alpha(i, a) = -U(i, c) * Z(r, c);
  }
  if(S.direct){
    Sig = sd.array().square().matrix().asDiagonal();
    for(int a = 0; a < S.ncorr; a++) Sig(S.ra[a], S.ca[a]) = Sig(S.ca[a], S.ra[a]) = THETA(S.oc + a);
  }

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(q, q);
  if(S.regmat.size() > 0) B = S.regmat;
  for(int b = 0; b < S.nreg; b++) B(S.bi[b], S.bj[b]) = THETA(S.ob + b);
  A = (Eigen::MatrixXd::Identity(q, q) - B).inverse();
  Lt = Lam * A;

  if(S.rescov.size() > 0){
    res = S.rescov;
    res.diagonal().setZero();
    for(int c = 1; c < p; c++){
      for(int r = 0; r < c; r++) if(S.resid(r, c) >= 0) res(r, c) = res(c, r) = THETA(S.resid(r, c));
    }
  }

  taus = tau;
  if(!S.plain){
    Eigen::VectorXd nu = S.constrint, ls = S.constrlogscale;
    kappa = S.constrmean;
    for(int j = 0; j < p; j++){
      if(S.nidx[j] >= 0) nu(j) = THETA(S.nidx[j]);
      if(S.scidx[j] >= 0) ls(j) = THETA(S.scidx[j]);
    }
    for(int a = 0; a < S.nmean; a++) kappa(S.mfree[a]) = THETA(S.oa + a);
    mu = nu + Lt * kappa;
    const Eigen::MatrixXd Psi = S.direct ? Sig : Eigen::MatrixXd(sd.asDiagonal() * U.transpose() * U * sd.asDiagonal());
    v = (Lt * Psi).cwiseProduct(Lt).rowwise().sum();
    if(S.thetapar){
      th = (2 * ls).array().exp();
      s = (v + th).array().rsqrt();
    }else{
      s = ls.array().exp();
    }
    for(int j = 0; j < p; j++){
      for(int c = 0; c < S.L.cat(j) - 1; c++) taus(S.L.tk(j) + c) = s(j) * (tau(S.L.tk(j) + c) - mu(j));
    }
  }

  J.resize(p);
  for(int j = 0; j < p; j++){
    const int nl = S.lidx[j].size();
    J[j].resize(q, nl + S.nreg);
    for(int a = 0; a < nl; a++) J[j].col(a) = A.row(S.lfac[j][a]).transpose();
    for(int b = 0; b < S.nreg; b++) J[j].col(nl + b) = Lt(j, S.bi[b]) * A.row(S.bj[b]).transpose();
  }
}
