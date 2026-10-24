#include "parallelWorkers.h"

void parallelWorkers::pairsGradient::operator()(std::size_t BEGIN, std::size_t END){
  // Run along the pairs identified by indexes in index_vector
  // computing their contribution to nll, gradient and eventually Hessian.
  int d = theta.size();
  //subset_gradient.resize(d); subset_gradient.setZero();
  for (int h = BEGIN; h < END; h++){

    // identify corresponding column in pairs_table
    const int col = index_vector[h];

    // identify the pair
    const int k = items_pairs(0, col);
    const int l = items_pairs(1, col);

    // initialize empty pair-output
    double pair_ll = 0;
    Eigen::VectorXd pair_gradient = Eigen::VectorXd::Zero(d);

    // computation of log-likelihood, gradient
    // rearrange parameters
    const int q = constrmat.cols();
    const Eigen::MatrixXd Lam             = params::loadings::theta2mat(theta, constrmat, llc, nthr, nload);
    const Eigen::MatrixXd Ru              = params::latvar::theta2cmat(theta, nthr, nload, ncorr, nvar, q);
    const Eigen::MatrixXd Du              = params::latvar::theta2dmat(theta, constrsd, nthr, nload, ncorr, nvar, q);
    const Eigen::MatrixXd Sigma_u         = Du * Ru * Du;
    const Eigen::VectorXd tau             = params::thresholds::theta2vec(theta, nthr);
    pairs::pair_contribution_extended(constrmat, constrsd, llc, c_vec, theta,
                                      Lam, Ru, Du, Sigma_u, tau,
                                      corrFLAG, nthr, nload, ncorr, nvar, k, l, pairs_table, silentFLAG, gradFLAG, pair_ll, pair_gradient);

    // update
                                      {
                                        subset_ll += pair_ll;
                                        subset_gradient += pair_gradient;
                                      }

  }
}

void parallelWorkers::pairsGradient::join(const pairsGradient &RHS){
  subset_ll += RHS.subset_ll;
  for(int i = 0; i < subset_gradient.size(); i++){
    subset_gradient(i) += RHS.subset_gradient(i);

  }
}




void parallelWorkers::patternwiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){
  // Run along the pairs identified by indexes in index_vector
  // computing their contribution to nll, gradient and eventually Hessian.
  const int d = theta.size();
  const int q = constrmat.cols();
  const int p = constrmat.rows();
  const int m = nload + ncorr + nvar;

  // rearrange parameters
  const Eigen::MatrixXd Lam = params::loadings::theta2mat(theta, constrmat, llc, nthr, nload);
  const Eigen::VectorXd sd  = params::latvar::theta2dmat(theta, constrsd, nthr, nload, ncorr, nvar, q).diagonal();
  const Eigen::VectorXd tau = params::thresholds::theta2vec(theta, nthr);

  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(q, q);
  std::vector<int> ra(ncorr), ca(ncorr);
  for(int c = 1, a = 0; c < q && a < ncorr; c++){
    for(int r = 0; r < c; r++, a++){
      Z(r, c) = std::tanh(theta(nthr + nload + a));
      ra[a] = r;
      ca[a] = c;
    }
  }
  const Eigen::MatrixXd C = (1 - Z.array().square()).sqrt().matrix();
  Eigen::MatrixXd P = Eigen::MatrixXd::Ones(q, q);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(q, q);
  for(int c = 0; c < q; c++){
    for(int r = 1; r <= c; r++) P(r, c) = P(r - 1, c) * C(r - 1, c);
    for(int r = 0; r < c; r++) U(r, c) = Z(r, c) * P(r, c);
    U(c, c) = P(c, c);
  }

  Eigen::MatrixXd alpha = Eigen::MatrixXd::Zero(q, ncorr);
  for(int a = 0; a < ncorr; a++){
    const int r = ra[a], c = ca[a];
    alpha(r, a) = P(r, c) * C(r, c) * C(r, c);
    for(int i = r + 1; i <= c; i++) alpha(i, a) = -U(i, c) * Z(r, c);
  }

  std::vector<int> sfree;
  for(int h = 0; h < q; h++) if(!std::isfinite(constrsd(h))) sfree.push_back(h);

  std::vector<Eigen::MatrixXd> JL(p, Eigen::MatrixXd::Zero(q, nload));
  for(int h = 0, idx = 0; h < q; h++){
    for(int j = 0; j < p; j++){
      if(!std::isfinite(constrmat(j, h))){
        JL[j](h, idx) = 1;
        idx++;
      }
    }
  }
  for(const auto &lc : llc){
    Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(nload);
    for(std::size_t v = 1; v < lc.size(); v++) row += lc.at(v).at(0) * JL[lc.at(v).at(1) - 1].row(lc.at(v).at(2) - 1);
    JL[lc.at(0).at(0) - 1].row(lc.at(0).at(1) - 1) = row;
  }

  std::vector<int> pk, pl, off;
  for(int k = 1, o = 0; k < p; k++){
    for(int l = 0; l < k; l++){
      pk.push_back(k);
      pl.push_back(l);
      off.push_back(o);
      o += c_vec(k) * c_vec(l);
    }
  }

  for(std::size_t r = BEGIN; r < END; r++){
    // identify the pair
    const int k = pk[r], l = pl[r];
    const int ck = c_vec(k), cl = c_vec(l);
    const int tk = c_vec.head(k).sum() - k, tl = c_vec.head(l).sum() - l;
    const int mt = ck + cl - 2;

    const Eigen::VectorXd x = sd.cwiseProduct(Lam.row(k).transpose());
    const Eigen::VectorXd y = sd.cwiseProduct(Lam.row(l).transpose());
    const Eigen::VectorXd X = U * x, Y = U * y;
    // implied model correlation
    const double rho = X.dot(Y);
    const Eigen::MatrixXd DJk = sd.asDiagonal() * JL[k], DJl = sd.asDiagonal() * JL[l];

    Eigen::MatrixXd dX(q, m), dY(q, m);
    dX.leftCols(nload) = U * DJk;
    dY.leftCols(nload) = U * DJl;
    for(int a = 0; a < ncorr; a++){
      dX.col(nload + a) = x(ca[a]) * alpha.col(a);
      dY.col(nload + a) = y(ca[a]) * alpha.col(a);
    }
    for(int s = 0; s < nvar; s++){
      dX.col(nload + ncorr + s) = x(sfree[s]) * U.col(sfree[s]);
      dY.col(nload + ncorr + s) = y(sfree[s]) * U.col(sfree[s]);
    }
    const Eigen::VectorXd grho = dX.transpose() * Y + dY.transpose() * X;
    Eigen::MatrixXd hrho = dX.transpose() * dY;
    hrho += hrho.transpose().eval();

    const Eigen::VectorXd aY = alpha.transpose() * Y, aX = alpha.transpose() * X;
    const Eigen::VectorXd uY = U.transpose() * Y, uX = U.transpose() * X;
    for(int i = 0; i < nload; i++){
      for(int a = 0; a < ncorr; a++){
        const double v = DJk(ca[a], i) * aY(a) + DJl(ca[a], i) * aX(a);
        hrho(i, nload + a) += v;
        hrho(nload + a, i) += v;
      }
      for(int s = 0; s < nvar; s++){
        const int h = sfree[s];
        const double v = DJk(h, i) * uY(h) + DJl(h, i) * uX(h);
        hrho(i, nload + ncorr + s) += v;
        hrho(nload + ncorr + s, i) += v;
      }
    }
    for(int a = 0; a < ncorr; a++){
      for(int b = a; b < ncorr && ca[b] == ca[a]; b++){
        const int c = ca[a], i0 = ra[a], i1 = ra[b];
        Eigen::VectorXd dd = Eigen::VectorXd::Zero(q);
        if(i0 < i1){
          dd(i1) = -P(i1, c) * Z(i0, c) * C(i1, c) * C(i1, c);
          for(int i = i1 + 1; i <= c; i++) dd(i) = U(i, c) * Z(i0, c) * Z(i1, c);
        }else{
          dd(i0) = -2 * P(i0, c) * Z(i0, c) * C(i0, c) * C(i0, c);
          for(int i = i0 + 1; i <= c; i++) dd(i) = U(i, c) * (2 * Z(i0, c) * Z(i0, c) - 1);
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
        if(ca[a] != sfree[s]) continue;
        hrho(nload + a, js) += grho(nload + a);
        hrho(js, nload + a) += grho(nload + a);
      }
    }

    const double r2 = 1 - rho * rho, sr = std::sqrt(r2);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
    double w = 0;
    for(int sk = 0; sk < ck; sk++){
      for(int sl = 0; sl < cl; sl++){
        const int col = off[r] + sk * cl + sl;
        const double nc = freq(4, col);
        // identify thresholds
        const Eigen::VectorXd th = params::extract_thresholds(tau, c_vec, k, l, sk, sl);
        // compute pi
        const double pe = biprobs::compute_pi(c_vec, th, rho, k, l, sk, sl) + 1e-8;

        Eigen::VectorXd g = Eigen::VectorXd::Zero(mt + 1);
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
        for(int e = 0; e < 2; e++){
          for(int f = 0; f < 2; f++){
            const double xv = th(2 * e), yv = th(2 * f + 1);
            const int ik = sk - e, il = sl - f;
            const int px = (ik >= 0 && ik < ck - 1) ? ik : -1;
            const int py = (il >= 0 && il < cl - 1) ? ck - 1 + il : -1;
            const double sg = (e == f) ? 1 : -1;
            const double wx = (yv - rho * xv) / sr, wy = (xv - rho * yv) / sr;
            const double phx = R::dnorm(xv, 0, 1, 0), phy = R::dnorm(yv, 0, 1, 0);
            const double Fx = phx * R::pnorm(wx, 0, 1, 1, 0), Fy = phy * R::pnorm(wy, 0, 1, 1, 0);
            const double f2 = binorm::dbvnorm(xv, yv, rho, 0);
            if(px >= 0){
              g(px) += sg * Fx;
              G(px, px) += sg * (-xv * Fx - rho / sr * phx * R::dnorm(wx, 0, 1, 0));
              G(px, mt) += -sg * f2 * (xv - rho * yv) / r2;
              G(mt, px) += -sg * f2 * (xv - rho * yv) / r2;
            }
            if(py >= 0){
              g(py) += sg * Fy;
              G(py, py) += sg * (-yv * Fy - rho / sr * phy * R::dnorm(wy, 0, 1, 0));
              G(py, mt) += -sg * f2 * (yv - rho * xv) / r2;
              G(mt, py) += -sg * f2 * (yv - rho * xv) / r2;
            }
            if(px >= 0 && py >= 0){
              G(px, py) += sg * f2;
              G(py, px) += sg * f2;
            }
            g(mt) += sg * f2;
            G(mt, mt) += sg * f2 * (rho / r2 + xv * yv / r2 - rho * (xv * xv - 2 * rho * xv * yv + yv * yv) / (r2 * r2));
          }
        }

        // update
        M += nc * (g * g.transpose() / (pe * pe) - G / pe);
        w += nc * g(mt) / pe;

        Eigen::VectorXd gcol = Eigen::VectorXd::Zero(d);
        gcol.segment(tk, ck - 1) = g.head(ck - 1);
        gcol.segment(tl, cl - 1) = g.segment(ck - 1, cl - 1);
        gcol.tail(m) = g(mt) * grho;
        gcol /= pe;
        gradmat.col(col) = gcol;
        gradient += nc * gcol;
      }
    }

    std::vector<int> idx(mt);
    for(int i = 0; i < ck - 1; i++) idx[i] = tk + i;
    for(int i = 0; i < cl - 1; i++) idx[ck - 1 + i] = tl + i;
    for(int i = 0; i < mt; i++){
      for(int j = 0; j < mt; j++) hessian(idx[i], idx[j]) += M(i, j);
      hessian.row(idx[i]).tail(m) += M(i, mt) * grho.transpose();
      hessian.col(idx[i]).tail(m) += M(i, mt) * grho;
    }
    hessian.bottomRightCorner(m, m) += M(mt, mt) * grho * grho.transpose() - w * hrho;
  }
}

void parallelWorkers::patternwiseOuterProd::join(const patternwiseOuterProd &RHS){
  gradient += RHS.gradient;
  hessian += RHS.hessian;
}


void parallelWorkers::casewiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){

  const int d = gradmat.rows();
  const int p = data.cols();
  //subset_gradient.resize(d); subset_gradient.setZero();
  for (int i = BEGIN; i < END; i++){

    Eigen::VectorXd gradienti = Eigen::VectorXd::Zero(d);

    for(int k = 1; k < p; k++){
      const int ck = c_vec(k);

      // identify column index in freq table
      // i1: starting index item k
      int i1 = 0;
      if(k > 1){
        for(unsigned int u = 1; u < k; u++){
          const int cu = c_vec(u);
          i1 += cu * c_vec.segment(0,u).sum();
        }
      }



      const int sk = data(i, k);

      for(int l = 0; l < k; l++){
        const int cl = c_vec(l);

        // i2 starting index from i1 dor item l
        int i2 = 0;
        if(l > 0){
          i2 = c_vec.segment(0,l).sum() * c_vec(k);
        }

        // i3: starting index from i2 for cat sk
        unsigned int i3 = sk * cl;


        const int sl = data(i, l);

        // final column index for pairs_tab. Print to check
        const int r = i1 + i2 + i3 + sl;

        gradienti += gradmat.col(r);






      }
    }

    gradient += gradienti;
    variability += gradienti * gradienti.transpose();

  }
}

void parallelWorkers::casewiseOuterProd::join(const casewiseOuterProd &RHS){
  gradient += RHS.gradient;
  variability += RHS.variability;
}
