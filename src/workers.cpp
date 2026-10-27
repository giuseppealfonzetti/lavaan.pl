#include "workers.h"

void workers::PL::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t i = BEGIN; i < END; i++) pairs::pair_contribution(st, freq, rs[i], order, ll, gradient);
}

void workers::PL::join(const PL &O){
  ll += O.ll;
  gradient += O.gradient;
}

void workers::Hess::operator()(std::size_t BEGIN, std::size_t END){
  Eigen::VectorXd g;
  Eigen::MatrixXd G2;
  for(std::size_t r = BEGIN; r < END; r++){
    const pairs::Pair PR(st, r, true);
    const pairs::Grid G(st.tau, PR.tk, PR.tl, PR.ck, PR.cl, PR.rho, 2);
    const Eigen::MatrixXd hrho = pairs::rho_hessian(st, PR);
    const double cap = pairs::cell_cap(st, freq, PR);
    const int mt = PR.mt, s = PR.s;

    std::vector<int> I;
    for(int a = 0; a < PR.ck - 1; a++) I.push_back(PR.tk + a);
    for(int a = 0; a < PR.cl - 1; a++) I.push_back(PR.tl + a);
    I.insert(I.end(), PR.idx.begin(), PR.idx.end());

    Eigen::MatrixXd sc(mt + s, PR.ck * PR.cl);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
    Eigen::VectorXd gl(mt + s), gsum = Eigen::VectorXd::Zero(mt + s);
    double w = 0;
    for(int sk = 0; sk < PR.ck; sk++){
      for(int sl = 0; sl < PR.cl; sl++){
        const int c = sk * PR.cl + sl;
        const double n = freq(PR.off + c);
        const double pe = G.cell(sk, sl, g, G2);
        const double ns = std::min(n, cap * pe);
        gl.head(mt) = g.head(mt) / pe;
        gl.tail(s) = g(mt) / pe * PR.grho;
        sc.col(c) = ns < n ? Eigen::VectorXd::Zero(mt + s) : gl;
        if(n == 0) continue;
        gsum += ns * gl;
        M += (ns < n ? 0 : n) * g * g.transpose() / (pe * pe) - ns * G2 / pe;
        w += ns * g(mt) / pe;
      }
    }

    Eigen::MatrixXd HL(mt + s, mt + s);
    HL.topLeftCorner(mt, mt) = M.topLeftCorner(mt, mt);
    HL.topRightCorner(mt, s) = M.col(mt).head(mt) * PR.grho.transpose();
    HL.bottomLeftCorner(s, mt) = HL.topRightCorner(mt, s).transpose();
    HL.bottomRightCorner(s, s) = M(mt, mt) * PR.grho * PR.grho.transpose() - w * hrho;
    for(int i = 0; i < mt + s; i++){
      gradient(I[i]) += gsum(i);
      for(int j = 0; j < mt + s; j++) hessian(I[i], I[j]) += HL(i, j);
    }
    if(!score.empty()){
      score[r] = sc;
      sidx[r] = I;
    }
  }
}

void workers::Hess::join(const Hess &O){
  hessian += O.hessian;
  gradient += O.gradient;
}

void workers::Cases::operator()(std::size_t BEGIN, std::size_t END){
  for(std::size_t u = BEGIN; u < END; u++){
    Eigen::VectorXd col = Eigen::VectorXd::Zero(out.rows());
    for(int r = 0; r < L.P; r++){
      const int k = L.pk[r], l = L.pl[r];
      if(std::isnan(y(u, k)) || std::isnan(y(u, l))) continue;
      const int c = y(u, k) * L.cat(l) + y(u, l);
      const std::vector<int> &I = sidx[r];
      for(std::size_t a = 0; a < I.size(); a++) col(I[a]) += score[r](a, c);
    }
    if(uw.size() > 0 && uw(u) > 0){
      for(int k = 0, uk = 0; k < L.p; uk += L.cat(k), k++){
        if(!std::isnan(y(u, k))) col.head(L.nthr) += uw(u) * uscore.col(uk + y(u, k));
      }
    }
    out.col(u) = col;
  }
}

void workers::Saturated::operator()(std::size_t BEGIN, std::size_t END){
  Eigen::VectorXd g;
  Eigen::MatrixXd G2;
  for(std::size_t r = BEGIN; r < END; r++){
    const int k = L.pk[r], l = L.pl[r], ck = L.cat(k), cl = L.cat(l), mt = ck + cl - 2;
    const pairs::Grid G(tau, L.tk(k), L.tk(l), ck, cl, rho(r), order);
    std::vector<int> idx(mt);
    for(int a = 0; a < ck - 1; a++) idx[a] = L.tk(k) + a;
    for(int a = 0; a < cl - 1; a++) idx[ck - 1 + a] = L.tk(l) + a;

    Eigen::VectorXd gpair = Eigen::VectorXd::Zero(mt + 1);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
    for(int sk = 0; sk < ck; sk++){
      for(int sl = 0; sl < cl; sl++){
        const int col = L.off[r] + sk * cl + sl;
        const double n = freq(col);
        const double pe = G.cell(sk, sl, g, G2);
        ll += n * std::log(pe);
        if(order > 0) gpair += n / pe * g;
        if(order == 2) M += n * (g * g.transpose() / (pe * pe) - G2 / pe);
        if(scores) patscore.col(col).head(mt + 1) = g / pe;
      }
    }
    if(order > 0){
      for(int i = 0; i < mt; i++) gradient(idx[i]) += gpair(i);
      gradient(tau.size() + r) += gpair(mt);
    }
    if(order == 2){
      for(int i = 0; i < mt; i++){
        for(int j = 0; j < mt; j++){
          htt(idx[i], idx[j]) += M(i, j);
          hschur(idx[i], idx[j]) += M(i, j) - M(i, mt) * M(j, mt) / M(mt, mt);
        }
        htr(idx[i], r) += M(i, mt);
      }
      hrr(r) += M(mt, mt);
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
  for(std::size_t i = BEGIN; i < END; i++){
    for(int r = 0; r < L.P; r++){
      const int k = L.pk[r], l = L.pl[r], ck = L.cat(k), cl = L.cat(l);
      if(std::isnan(y(i, k)) || std::isnan(y(i, l))) continue;
      const int col = L.off[r] + y(i, k) * cl + y(i, l);
      for(int a = 0; a < ck - 1; a++) st(L.tk(k) + a, i) += patscore(a, col);
      for(int b = 0; b < cl - 1; b++) st(L.tk(l) + b, i) += patscore(ck - 1 + b, col);
      sr(r, i) = patscore(ck + cl - 2, col);
    }
    if(uw.size() > 0 && uw(i) > 0){
      for(int k = 0, uk = 0; k < L.p; uk += L.cat(k), k++){
        if(!std::isnan(y(i, k))) st.col(i) += uw(i) * uscore.col(uk + y(i, k));
      }
    }
  }
}
