#include "parallelWorkers.h"

void parallelWorkers::pairsGradient::operator()(std::size_t BEGIN, std::size_t END){
  // Run along the pairs identified by indexes in index_vector
  // computing their contribution to nll, gradient and eventually Hessian.
  int d = theta.size();

  // rearrange parameters
  const pairs::Setup S(theta, constrmat, constrsd, nthr, nload, ncorr, nvar, corrFLAG);

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
    pairs::pair_contribution(S, c_vec, pairs_table, k, l, gradFLAG, pair_ll, pair_gradient);

    // update
    subset_ll += pair_ll;
    subset_gradient += pair_gradient;
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
  const int p = constrmat.rows();

  // rearrange parameters
  const pairs::Setup S(theta, constrmat, constrsd, nthr, nload, ncorr, nvar, corrFLAG);
  const int m = S.m;

  std::vector<int> pk, pl;
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++){
      pk.push_back(k);
      pl.push_back(l);
    }
  }

  Eigen::VectorXd g;
  Eigen::MatrixXd G;
  for(std::size_t r = BEGIN; r < END; r++){
    // identify the pair
    const pairs::Pair PR(S, c_vec, pk[r], pl[r], true);
    const int ck = PR.ck, cl = PR.cl, tk = PR.tk, tl = PR.tl, mt = PR.mt;
    const Eigen::MatrixXd hrho = pairs::rho_hessian(S, PR);

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
    double w = 0;
    for(int sk = 0; sk < ck; sk++){
      for(int sl = 0; sl < cl; sl++){
        const int col = PR.offset + sk * cl + sl;
        const double nc = freq(4, col);
        // compute pi
        const double pe = pairs::pattern_derivatives(S.tau, c_vec, PR, sk, sl, 2, g, G);

        // update
        M += nc * (g * g.transpose() / (pe * pe) - G / pe);
        w += nc * g(mt) / pe;

        Eigen::VectorXd gcol = Eigen::VectorXd::Zero(d);
        gcol.segment(tk, ck - 1) = g.head(ck - 1);
        gcol.segment(tl, cl - 1) = g.segment(ck - 1, cl - 1);
        gcol.tail(m) = g(mt) * PR.grho;
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
      hessian.row(idx[i]).tail(m) += M(i, mt) * PR.grho.transpose();
      hessian.col(idx[i]).tail(m) += M(i, mt) * PR.grho;
    }
    hessian.bottomRightCorner(m, m) += M(mt, mt) * PR.grho * PR.grho.transpose() - w * hrho;
  }
}

void parallelWorkers::patternwiseOuterProd::join(const patternwiseOuterProd &RHS){
  gradient += RHS.gradient;
  hessian += RHS.hessian;
}


void parallelWorkers::casewiseOuterProd::operator()(std::size_t BEGIN, std::size_t END){

  const int d = gradmat.rows();
  const int p = data.cols();
  for (int i = BEGIN; i < END; i++){

    Eigen::VectorXd gradienti = Eigen::VectorXd::Zero(d);

    for(int k = 1; k < p; k++){

      // identify column index in freq table
      // i1: starting index item k
      int i1 = 0;
      if(k > 1){
        for(unsigned int u = 1; u < k; u++){
          const int cu = c_vec(u);
          i1 += cu * c_vec.segment(0,u).sum();
        }
      }



      if(std::isnan(data(i, k))) continue;
      const int sk = data(i, k);

      for(int l = 0; l < k; l++){
        if(std::isnan(data(i, l))) continue;
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

    if(uw.size() > 0 && uw(i) > 0){
      for(int k = 0, uk = 0; k < p; uk += c_vec(k), k++){
        if(!std::isnan(data(i, k))) gradienti.head(uscore.rows()) += uw(i) * uscore.col(uk + data(i, k));
      }
    }

    const double gwi = gw.size() > 0 ? gw(i) : 1;
    gradient += gwi * gradienti;
    if(cl.size() > 0){
      ztot.col(cl(i)) += gwi * gradienti;
    }else{
      const double jwi = jw.size() > 0 ? jw(i) : 1;
      variability += jwi * jwi * gradienti * gradienti.transpose();
    }

  }
}

void parallelWorkers::casewiseOuterProd::join(const casewiseOuterProd &RHS){
  gradient += RHS.gradient;
  variability += RHS.variability;
  ztot += RHS.ztot;
}


void parallelWorkers::saturatedPairs::operator()(std::size_t BEGIN, std::size_t END){
  const int p = c_vec.size();
  const int nthr = tau.size();

  std::vector<int> pk, pl;
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++){
      pk.push_back(k);
      pl.push_back(l);
    }
  }

  Eigen::VectorXd g;
  Eigen::MatrixXd G;
  for(std::size_t r = BEGIN; r < END; r++){
    // identify the pair
    const pairs::Pair PR(c_vec, pk[r], pl[r], rho(r));
    const int ck = PR.ck, cl = PR.cl, mt = PR.mt;

    std::vector<int> idx(mt);
    for(int i = 0; i < ck - 1; i++) idx[i] = PR.tk + i;
    for(int i = 0; i < cl - 1; i++) idx[ck - 1 + i] = PR.tl + i;

    Eigen::VectorXd gpair = Eigen::VectorXd::Zero(mt + 1);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(mt + 1, mt + 1);
    for(int sk = 0; sk < ck; sk++){
      for(int sl = 0; sl < cl; sl++){
        const int col = PR.offset + sk * cl + sl;
        const double nc = freq(4, col);
        const double pe = pairs::pattern_derivatives(tau, c_vec, PR, sk, sl, order, g, G);

        ll += nc * log(pe);
        if(order > 0) gpair += nc / pe * g;
        if(order == 2) M += nc * (g * g.transpose() / (pe * pe) - G / pe);
        if(scores) patscore.col(col).head(mt + 1) = g / pe;
      }
    }

    if(order > 0){
      for(int i = 0; i < mt; i++) gradient(idx[i]) += gpair(i);
      gradient(nthr + r) += gpair(mt);
    }
    if(order == 2){
      for(int i = 0; i < mt; i++){
        for(int j = 0; j < mt; j++) htt(idx[i], idx[j]) += M(i, j);
        htr(idx[i], r) += M(i, mt);
      }
      hrr(r) += M(mt, mt);
    }
  }
}

void parallelWorkers::saturatedPairs::join(const saturatedPairs &RHS){
  ll += RHS.ll;
  gradient += RHS.gradient;
  htt += RHS.htt;
  htr += RHS.htr;
  hrr += RHS.hrr;
}


void parallelWorkers::saturatedScores::operator()(std::size_t BEGIN, std::size_t END){
  const int p = data.cols();

  std::vector<pairs::Pair> prs;
  for(int k = 1; k < p; k++){
    for(int l = 0; l < k; l++) prs.emplace_back(c_vec, k, l, 0);
  }

  for(std::size_t i = BEGIN; i < END; i++){
    for(std::size_t r = 0; r < prs.size(); r++){
      const pairs::Pair &PR = prs[r];
      if(std::isnan(data(i, PR.k)) || std::isnan(data(i, PR.l))) continue;
      const int col = PR.offset + data(i, PR.k) * PR.cl + data(i, PR.l);
      for(int a = 0; a < PR.ck - 1; a++) st(PR.tk + a, i) += patscore(a, col);
      for(int b = 0; b < PR.cl - 1; b++) st(PR.tl + b, i) += patscore(PR.ck - 1 + b, col);
      sr(r, i) = patscore(PR.mt, col);
    }
    if(uw.size() > 0 && uw(i) > 0){
      for(int k = 0, uk = 0; k < p; uk += c_vec(k), k++){
        if(!std::isnan(data(i, k))) st.col(i) += uw(i) * uscore.col(uk + data(i, k));
      }
    }
  }
}
