#include "workers.h"

void workers::PL::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t i = BEGIN; i < END; i++) pairs::pair_contribution(st, data, rs[i], order, ll, gradient);
}

void workers::PL::join(const PL &O){
  ll += O.ll;
  gradient += O.gradient;
}

void workers::Hess::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t r = BEGIN; r < END; r++){
    const pairs::Pair PR(st, r, true);
    const int mt = PR.mt, s = PR.s;
    const std::vector<int> I = PR.local();

    const bool keep = !score.empty() && data.tables && !PR.numk && !PR.numl && data.x.cols() * st.S.L.ncell * (mt + s) <= 5e7;
    Eigen::MatrixXd sc(keep ? mt + s : 0, keep ? data.x.cols() * PR.ck * PR.cl : 0);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(PR.m, PR.m);
    Eigen::VectorXd gu = Eigen::VectorXd::Zero(PR.m);
    double ll = 0;
    pairs::moments(st.S.L, pairs::Reduced(st), data, PR, 2, ll, gu, M, keep ? &sc : nullptr);
    const Eigen::VectorXd gsum = PR.lift(gu);

    Eigen::MatrixXd HL(mt + s, mt + s);
    if(PR.plain){
      HL.topLeftCorner(mt, mt) = M.topLeftCorner(mt, mt);
      HL.topRightCorner(mt, s) = M.col(mt).head(mt) * PR.grho.transpose();
      HL.bottomLeftCorner(s, mt) = HL.topRightCorner(mt, s).transpose();
      HL.bottomRightCorner(s, s) = M(mt, mt) * PR.grho * PR.grho.transpose() - gu(mt) * pairs::rho_hessian(st, PR);
    }else{
      HL = PR.Ju.transpose() * M * PR.Ju - pairs::curvature(st, PR, gu);
    }
    for(int i = 0; i < mt + s; i++){
      gradient(I[i]) += gsum(i);
      for(int j = 0; j < mt + s; j++) hessian(I[i], I[j]) += HL(i, j);
    }
    if(!score.empty()){
      score[r] = std::move(sc);
      sidx[r] = I;
      pgrad[r] = gsum;
    }
  }
}

void workers::Hess::join(const Hess &O){
  hessian += O.hessian;
  gradient += O.gradient;
}

void workers::Cases::operator()(std::size_t BEGIN, std::size_t END){
  const model::Layout &L = st.S.L;
  const pairs::Reduced R(st);
  pairs::Work W;
  for(int r = 0; r < L.P; r++){
    const int k = L.pk[r], l = L.pl[r], ck = L.cat(k), cl = L.cat(l);
    const std::vector<int> &I = sidx[r];
    if(score[r].size() > 0){
      for(std::size_t u = BEGIN; u < END; u++){
        if(std::isnan(data.y(u, k)) || std::isnan(data.y(u, l))) continue;
        const int c = data.stratum(u) * ck * cl + data.y(u, k) * cl + data.y(u, l);
        for(std::size_t a = 0; a < I.size(); a++) out(I[a], u) += score[r](a, c);
      }
      continue;
    }
    const pairs::Pair PR(st, r, true);
    pairs::Grid G(ck, cl, PR.rho, 1);
    Eigen::VectorXd gc(PR.m);
    for(std::size_t u = BEGIN; u < END; u++){
      if(std::isnan(data.y(u, k)) || std::isnan(data.y(u, l))) continue;
      pairs::unit(R, data, PR, G, u, 1, W);
      gc.setZero();
      PR.slots(data.x.col(data.stratum(u)), data.y(u, k), data.y(u, l), W.wk, W.wl);
      PR.add_score(1, W.g, W.wk, W.wl, gc);
      PR.add_scale(1, R.scale, gc, nullptr);
      const Eigen::VectorXd v = PR.lift(gc);
      for(std::size_t a = 0; a < I.size(); a++) out(I[a], u) += v(a);
    }
  }
  if(uw.size() == 0) return;
  for(std::size_t u = BEGIN; u < END; u++){
    if(uw(u) == 0) continue;
    const int ou = data.stratum(u) * L.cat.sum();
    for(int k = 0, uk = 0; k < L.p; uk += L.cat(k), k++){
      if(!std::isnan(data.y(u, k))) out.col(u) += uw(u) * uscore.col(ou + uk + data.y(u, k));
    }
  }
}

static std::vector<int> saturated_index(const model::Layout &L, const pairs::Pair &PR, const int NTHR, const int NX){
  std::vector<int> idx(PR.m);
  const int ot = NTHR + NX * L.p, os = ot + L.nnum;
  for(int a = 0; a < PR.ck - 1; a++) idx[a] = PR.tk + a;
  for(int a = 0; a < PR.cl - 1; a++) idx[PR.lk + a] = PR.tl + a;
  if(PR.numk) idx[0] = ot + L.inum(PR.k);
  if(PR.numl) idx[PR.lk] = ot + L.inum(PR.l);
  idx[PR.ml] = os + L.nnum + PR.r;
  const int ok = PR.ml + 1, ol = ok + PR.nsk;
  for(int c = 0; c < NX; c++){
    idx[ok + c] = NTHR + c * L.p + PR.k;
    idx[ol + c] = NTHR + c * L.p + PR.l;
  }
  if(PR.numk) idx[ok + NX] = os + L.inum(PR.k);
  if(PR.numl) idx[ol + NX] = os + L.inum(PR.l);
  return idx;
}

void workers::Saturated::operator()(std::size_t BEGIN, std::size_t END){
  const int nthr = reduced.tau.size(), nx = reduced.pis.cols();
  Eigen::VectorXd g;
  Eigen::MatrixXd G2;
  for(std::size_t r = BEGIN; r < END; r++){
    const pairs::Pair PR(L, r, rho(r), nx);
    const int rr = PR.ml, m = PR.m;
    const std::vector<int> idx = saturated_index(L, PR, nthr, nx);
    Eigen::VectorXd gu = Eigen::VectorXd::Zero(m);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m, m);
    pairs::moments(L, reduced, data, PR, order, ll, gu, M);
    if(order > 0){
      for(int i = 0; i < m; i++) gradient(idx[i]) += gu(i);
    }
    if(order == 2){
      const double hr = M(rr, rr);
      for(int i = 0; i < m; i++){
        if(i == rr) continue;
        for(int j = 0; j < m; j++){
          if(j == rr) continue;
          htt(idx[i], idx[j]) += M(i, j);
          hschur(idx[i], idx[j]) += M(i, j) - M(i, rr) * M(j, rr) / hr;
        }
        htr(idx[i], r) += M(i, rr);
      }
      hrr(r) += hr;
    }
    if(scores && data.tables && !PR.numk && !PR.numl){
      for(int t = 0; t < data.x.cols(); t++){
        const pairs::Grid Gs(reduced.tau, PR.tk, PR.tl, PR.ck, PR.cl, PR.rho, 1, data.x.col(t).dot(reduced.pis.row(PR.k).transpose()),
                             data.x.col(t).dot(reduced.pis.row(PR.l).transpose()));
        for(int sk = 0; sk < PR.ck; sk++){
          for(int sl = 0; sl < PR.cl; sl++){
            const double pe = Gs.cell(sk, sl, g, G2);
            patscore.col(t * L.ncell + PR.off + sk * PR.cl + sl).head(PR.mt + 1) = g / pe;
          }
        }
      }
    }
  }
}

void workers::Saturated::join(const Saturated &O){
  ll += O.ll;
  gradient += O.gradient;
  htt += O.htt;
  htr += O.htr;
  hrr += O.hrr;
  hschur += O.hschur;
}

void workers::SaturatedCases::operator()(std::size_t BEGIN, std::size_t END){
  const int nthr = reduced.tau.size(), nx = reduced.pis.cols(), ni = st.rows();
  pairs::Work W;
  for(int r = 0; r < L.P; r++){
    const pairs::Pair PR(L, r, rho(r), nx);
    const int k = PR.k, l = PR.l;
    const bool tab = data.tables && !PR.numk && !PR.numl;
    pairs::Grid G(PR.ck, PR.cl, PR.rho, 1);
    const std::vector<int> idx = saturated_index(L, PR, nthr, nx);
    Eigen::VectorXd gu(PR.m);
    for(std::size_t u = BEGIN; u < END; u++){
      if(std::isnan(data.y(u, k)) || std::isnan(data.y(u, l))) continue;
      const int t = data.stratum(u);
      if(tab){
        const int sk = data.y(u, k), sl = data.y(u, l);
        W.g = patscore.col(t * L.ncell + PR.off + sk * PR.cl + sl).head(PR.ml + 1);
      }else{
        pairs::unit(reduced, data, PR, G, u, 1, W);
      }
      gu.setZero();
      PR.slots(data.x.col(t), data.y(u, k), data.y(u, l), W.wk, W.wl);
      PR.add_score(1, W.g, W.wk, W.wl, gu);
      PR.add_scale(1, reduced.scale, gu, nullptr);
      for(int i = 0; i < PR.m; i++){
        if(idx[i] < ni) st(idx[i], u) += gu(i);
      }
      sr(r, u) = gu(PR.ml);
    }
  }
  if(uw.size() == 0) return;
  for(std::size_t u = BEGIN; u < END; u++){
    if(uw(u) == 0) continue;
    const int ou = data.stratum(u) * L.cat.sum();
    for(int k = 0, uk = 0; k < L.p; uk += L.cat(k), k++){
      if(!std::isnan(data.y(u, k))) st.col(u) += uw(u) * uscore.col(ou + uk + data.y(u, k));
    }
  }
}
