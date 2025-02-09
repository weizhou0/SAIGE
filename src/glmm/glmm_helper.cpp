#include <armadillo>

struct SA_NULL {
    arma::mat XV;
    arma::mat XVX;
    arma::mat XXVX_inv;
    arma::mat XVX_inv;
    arma::vec S_a;
    arma::mat XVX_inv_XV;
    arma::vec V;
};

SA_NULL ScoreTest_NULL_Model(const arma::vec& mu, const arma::vec& mu2, const arma::vec& y, const arma::mat& X) {
    SA_NULL result;

    result.V = mu2;
    arma::vec res = y - mu;
    result.XV = X.each_col() % result.V;
    result.XVX = X.t() * result.XV;
    result.XVX_inv = arma::inv(result.XVX);
    result.XXVX_inv = X * result.XVX_inv;
    result.XVX_inv_XV = result.XXVX_inv.each_col() % result.V;
    result.S_a = arma::sum(X.each_col() % res, 0).t();

    return result;
}


void getAIScore(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, const arma::ivec& fixtauVec,
                const arma::vec& Sigma_iY, const arma::mat& Sigma_iX, const arma::mat& cov, arma::mat& AI, arma::vec& YPAPY, arma::vec& Trace,
                arma::vec& PY, unsigned int k1, int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, bool LOCO,
                const arma::ivec& Ivec_start_indices, const arma::mat& eMat) {

    fixtauVec.print("fixtauVec");
    int q2 = arma::sum(fixtauVec == 0);
    arma::uvec idxtau = arma::find(fixtauVec == 0);
    arma::vec tau0;
    AI.set_size(k1, k1);
    YPAPY.zeros(k1);
    Trace.zeros(k1);
    arma::mat Sigma_iXt = Sigma_iX.t();
    arma::mat Xmatt = Xmat.t();

    PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));

    unsigned int n = PY.n_elem;
    arma::vec PAPY_1, PAPY, APY;
    arma::mat APYmat(n, k1);

    arma::vec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;

    if(eMat.n_nrows == 0){
        if(Ivec_start_indices.n_elem == 0){
            for (int i = 0; i < k1; i++) {
                if (fixtauVec(i) == 0) {
                    if (i == 0) {
                        APY = PY;
                    } else if (i == 1) {
                        APY = getCrossprodMatAndKin(PY, LOCO);
                    } 
                    APYmat.col(i) = APY;
                    PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                    PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                    for (int j = 0; j <= i; j++) {
                        AI(i, j) = arma::dot(APYmat.col(j), PAPY);
                        if (j != i) {
                            AI(j, i) = AI(i, j);
                        }
                    }
                    YPAPY(i) = arma::dot(PY, APYmat.col(i));
                }
            }

        }else{
            for (int i = 0; i < k1; i++) {
                if (fixtauVec(i) == 0) {
                    if (i == 0) {
                        APY = PY;
                    } else if (i == 1) {
                        if(g_isGRM){
                            GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                            APY = g_I_longl_mat * GRM_I_bvec;
                          }else{
                            APY = g_I_longl_mat * Ibvec;
                          }
                    } 
                    APYmat.col(i) = APY;
                    PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                    PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                    for (int j = 0; j <= i; j++) {
                        AI(i, j) = arma::dot(APYmat.col(j), PAPY);
                        if (j != i) {
                            AI(j, i) = AI(i, j);
                        }
                    }
                    YPAPY(i) = arma::dot(PY, APYmat.col(i));
                }
            }

        }        
    }else{  //with eMat
        if(Ivec_start_indices.n_elem == 0){
            for(int i=0; i<k1; i++){
                if(fixtauVec(i) == 0){
                    if(i==0){
                            APY = PY1;
                    }else if(i==1){
                            if(g_isGRM){
                                    APY = getCrossprodMatAndKin(PY1, LOCO);
                            }else{
                                    APY = getprod_eMat(PY1);
                            }
                    }else if(i==2){
                            if(g_isGRM){
                                    APY = getCrossprodMatAndKin_eMat(PY1, LOCO);
                            }
                    }
                    APYmat.col(i) = APY;
                    PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                    PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                    for (int j = 0; j <= i; j++) {
                        AI(i, j) = arma::dot(APYmat.col(j), PAPY);
                        if (j != i) {
                            AI(j, i) = AI(i, j);
                        }
                    }
                    YPAPY(i) = arma::dot(PY, APYmat.col(i));
                }
            }

        }else{ //with eMat and Ivec_start_indices
            for (int i = 0; i < k1; i++) {
                if (fixtauVec(i) == 0) {
                    if (i == 0) {
                        APY = PY;
                    } else if (i == 1) {
                        if(g_isGRM){
                            APY = g_I_longl_mat * GRM_I_bvec;
                          }else{
                            APY = getprodImatImattbVec(PY1);
                          }
                    } else if (i == 2) {
                        if(g_isGRM){
                            APY = getCrossprodMatAndKin_eMat_Imat(PY1, LOCO);
                        }else{
                           APY = getCrossprodMatAndI_eMat_Imat(PY1, LOCO);
                        }
                    }
                    APYmat.col(i) = APY;
                    PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                    PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                    for (int j = 0; j <= i; j++) {
                        AI(i, j) = arma::dot(APYmat.col(j), PAPY);
                        if (j != i) {
                            AI(j, i) = AI(i, j);
                        }
                    }
                    YPAPY(i) = arma::dot(PY, APYmat.col(i));
                }
            }     
    }
    GetTrace(Sigma_iX, Xmat, wVec, tauVec, fixtauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO, Trace, Ivec_start_indices);
}


void GetTrace(arma::mat Sigma_iX, arma::mat& Xmat, arma::vec& wVec, arma::vec& tauVec, arma::ivec & fixtauVec, arma::mat& cov1, int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, bool LOCO, arma::vec& trace, const arma::uvec & Ivec_start_indices) {

    set_seed(200);

    int q2 = arma::sum(fixtauVec == 0);
    arma::uvec idxtau = arma::find(fixtauVec == 0);

    idxtau.print("idxtau");

    arma::mat Sigma_iXt = Sigma_iX.t();
    int Nnomissing = wVec.n_elem;
    unsigned int k1 = g_num_Kmat;
    arma::mat temp_mat(nrun, k1);
    arma::mat temp_mat_update(nrun, q2);
    arma::vec temp_vec(nrun);
    arma::vec temp_vec_double(Nnomissing);
    temp_mat.zeros();
    temp_mat_update.zeros();

    arma::vec Sigma_iu;
    arma::vec Pu;
    arma::mat Au_mat(Nnomissing, k1);
    arma::vec uVec, uVec0;

    int nrun_trace_start = 0;
    int nrun_trace_end = nrun;
    arma::vec traceCV(q2);
    traceCV.fill(traceCVcutoff + 0.1);

    arma::vec traceCVsub;
    bool isConverge = false;

    arma::vec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
    while (!isConverge) {
        for (int i = nrun_trace_start; i < nrun_trace_end; i++) {

            uVec0 = arma::randn(Nnomissing);
            uVec = uVec0 * 2 - 1;

            Sigma_iu = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, uVec, maxiterPCG, tolPCG, LOCO);
            Pu = Sigma_iu - Sigma_iX * (cov1 * (Sigma_iXt * uVec));

            if (fixtauVec(0) == 0) {
                Au_mat.col(0) = uVec;
                temp_mat(i, 0) = dot(Au_mat.col(0), Pu);
            }

            if (Ivec_start_indices.n_elem == 0) {
                if(eMat.n_rows != 0){
                    if (fixtauVec(1) == 0) {
                        if (g_isGRM) {
                            temp_vec_double = getCrossprodMatAndKin(uVec, LOCO);
                        } else {
                            temp_vec_double = getCrossprodMatAndI_eMat(uVec, LOCO);
                        }
                        Au_mat.col(1) = temp_vec_double;
                        temp_mat(i, 1) = dot(temp_vec_double, Pu);
                    }
       
                    if (fixtauVec(2) == 0) {
                         if (g_isGRM) {
                            temp_vec_double = getCrossprodMatAndKin_eMat(uVec, LOCO);
                        } 
                        Au_mat.col(2) = temp_vec_double;
                        temp_mat(i, 2) = dot(temp_vec_double, Pu);
                    }
                }else{ //no eMat
                    if (g_isGRM) {
                        if (fixtauVec(1) == 0) {
                                temp_vec_double = getCrossprodMatAndKin(uVec, LOCO);
                        } 
                        Au_mat.col(1) = temp_vec_double;
                        temp_mat(i, 1) = dot(temp_vec_double, Pu);
                    }
                    
                }

            } else { //with Ivec_start_indices
                if (g_isGRM) {
                    Ibvec = g_I_longl_mat_t * uVec;
                    GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                }
                    
                if (g_isGRM) {
                    if (fixtauVec(1) == 0) {
                        temp_vec_double = g_I_longl_mat * GRM_I_bvec;
                        Au_mat.col(1) = temp_vec_double;
                        temp_mat(i, 1) = dot(temp_vec_double, Pu);
                    }
                    if(eMat.n_rows != 0){
                        if (fixtauVec(2) == 0) {
                            temp_vec_double = getCrossprodMatAndKin_eMat_Imat(uVec, LOCO);
                            Au_mat.col(2) = temp_vec_double;
                            temp_mat(i, 2) = dot(temp_vec_double, Pu);
                        }
                    }

                } else {
                    if (fixtauVec(1) == 0) {
                        temp_vec_double = getprodImatImattbVec(uVec);
                        Au_mat.col(1) = temp_vec_double;
                        temp_mat(i, 1) = dot(temp_vec_double, Pu);
                    }
                    if(eMat.n_rows != 0){
                        if (fixtauVec(2) == 0) {
                            temp_vec_double = getCrossprodMatAndI_eMat_Imat(uVec, LOCO);
                            Au_mat.col(2) = temp_vec_double;
                            temp_mat(i, 2) = dot(temp_vec_double, Pu);
                        }
                    }
                }
            }

            uVec.clear();
            Pu.clear();
            Sigma_iu.clear();
        }

        temp_mat_update = temp_mat.cols(idxtau);

        std::cout << "dim temp_mat_update" << temp_mat_update.n_rows << " " << temp_mat_update.n_cols << std::endl;

        for (int k = 0; k < q2; k++) {
            temp_vec = temp_mat_update.col(k);
            traceCV(k) = calCV(temp_vec);
        }


        if (arma::any(traceCV > traceCVcutoff)) {
            isConverge = false;
        } else {
            isConverge = true;
        }


        if (!isConverge) {
            nrun_trace_start = nrun_trace_end;
            nrun_trace_end = nrun_trace_end + 10;
            temp_mat.resize(nrun_trace_end, k1);
            temp_mat_update.resize(nrun_trace_end, q2);
            Rcout << "CV for trace random estimator using " << nrun_trace_start << " runs is " << traceCV << " > " << traceCVcutoff << std::endl;
            Rcout << "try " << nrun_trace_end << " runs" << std::endl;
        }
    }

    Au_mat.clear();
    Pu.clear();
    Sigma_iu.clear();
    uVec.clear();
    temp_vec.clear();

    trace = arma::mean(temp_mat_update, 0).t();
    temp_mat.clear();
    temp_mat_update.clear();
}

// [[Rcpp::export]]
void set_g_omp_num_threads(unsigned int t_omp_num_threads){
g_omp_num_threads = t_omp_num_threads;
}

// [[Rcpp::export]]
arma::vec getprodImatImattbVec(arma::vec & bVec){
  omp_set_num_threads(g_omp_num_threads);
  auto n = g_I_start_indices.n_elem - 1;
  arma::vec resultVec(bVec.n_elem, arma::fill::zeros);
  #pragma omp parallel
  {
          auto thread_idx = omp_get_thread_num();
          for(int j = thread_idx; j < n; j += g_omp_num_threads) {
                  double sum = 0;
                  size_t start = g_I_start_indices[j];
                  size_t end = g_I_start_indices[j + 1];
                  for(size_t k = start; k < end; k++) {
                          sum += bVec[k];
                  }
                  for(size_t k = start; k < end; k++) {
                          resultVec[k] += sum;
                  }
          }
  }
  return(resultVec);
}




void getCoefficients_multiV(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, int maxiterPCG, double tolPCG, bool LOCO, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov, arma::vec& alpha, arma::vec& eta) {
    int Nnomissing = wVec.n_elem;
    std::cout << "before Sigma_iY" << std::endl;
    std::cout << "Yvec.n_elem " << Yvec.n_elem << std::endl;
    Sigma_iY = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, Yvec, maxiterPCG, tolPCG, LOCO);
    std::cout << "after Sigma_iY" << std::endl;
    int colNumX = Xmat.n_cols;
    Sigma_iX.set_size(Nnomissing, colNumX);
    arma::vec XmatVecTemp;
    for (int i = 0; i < colNumX; i++) {
        XmatVecTemp = Xmat.col(i);
        Sigma_iX.col(i) = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG, LOCO);
    }
    arma::mat Xmatt = Xmat.t();
    try {
        cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
    } catch (const std::exception& e) {
        cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
        std::cout << "inv_sympd failed, inverted with pinv" << std::endl;
    }
    cov.print("cov");
    arma::mat Sigma_iXt = Sigma_iX.t();
    arma::vec SigmaiXtY = Sigma_iXt * Yvec;
    alpha = cov * SigmaiXtY;
    eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
}




void Get_Coef(const arma::vec& y, const arma::mat& X, const arma::vec& tau, const arma::vec& offset, const arma::vec& var_weights, 
                    GLMResults* fit0ptr, arma::vec& alpha, arma::vec& eta, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov, 
                     arma::vec& Y, int maxiterPCG, double tolPCG, int maxiter, bool verbose, bool LOCO) {
    double tol_coef = 0.1;
    arma::vec mu = fit0ptr->linkinv(eta); 
    arma::vec mu_eta = fit0ptr->mu_eta(eta);
    Y = eta - offset + (y - mu) / mu_eta;

    if (var_weights.is_empty()) {
        var_weights = arma::ones(y.n_elem);
    } 
    arma::vec sqrtW = mu_eta / arma::sqrt(1 / var_weights % fit0ptr->variance(mu));
    arma::vec W = arma::square(sqrtW);

    for (int i = 0; i < maxiter; ++i) {
        // Assuming getCoefficients_multiV is a function that returns a struct with alpha, eta, cov, Sigma_iY, Sigma_iX
        getCoefficients_multiV(Y, X, W, tau, maxiterPCG, tolPCG, LOCO, Sigma_iY, Sigma_iX, cov, alpha, eta);
        eta = eta + offset;

     
        std::cout << "Tau:" << std::endl;
        tau.print();
        std::cout << "Fixed-effect coefficients:" << std::endl;
        alpha.print();


        mu = fit0ptr->linkinv(eta); // Assuming fit0ptr is the family$linkinv(eta)
        mu_eta = fit0ptr->mu_eta(eta); // Assuming fit0ptr is the family$mu.eta(eta)

        Y = eta - offset + (y - mu) / mu_eta;

        sqrtW = mu_eta / arma::sqrt(1 / var_weights % fit0ptr->variance(mu));
        W = arma::square(sqrtW);

        if (arma::max(arma::abs(alpha - alpha0) / (arma::abs(alpha) + arma::abs(alpha0) + tol_coef)) < tol_coef) {
            break;
        }
        alpha0 = alpha;
    }

    std::cout << "alpha" << std::endl;
    alpha.print();
}
   

void fitglmmaiRPCG_multiV_updateTau(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, const arma::ivec& fixtauVec, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov,
    arma::vec& alpha, int nrun, int maxiterPCG, double tolPCG, double tol, double traceCVcutoff, bool LOCO) {
    
    unsigned int k1 = g_num_Kmat;
    int q2 = arma::sum(fixtauVec == 0);
    arma::uvec idxtau = arma::find(fixtauVec == 0);
    arma::vec tau0 = tauVec;

    getAIScore(Yvec, Xmat, wVec, tauVec, fixtauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO, Ivec_start_indices, eMat);

    arma::vec YPAPY = re["YPAPY"];
    arma::vec Trace = re["Trace"];
    arma::vec score1 = YPAPY - Trace;
    arma::mat AI1 = re["AI"];
    arma::vec Dtau;

    try {
        Dtau = arma::solve(AI1, score1, arma::solve_opts::allow_ugly);
    } catch (std::runtime_error&) {
        std::cout << "arma::solve(AI, score): AI seems singular, using less variant components matrix is suggested." << std::endl;
        Dtau.zeros();
    }

    arma::vec Dtau_k1(k1, arma::fill::zeros);
    int i2 = 0;
    for (int i = 0; i < k1; i++) {
        if (fixtauVec(i) == 0) {
            Dtau_k1(i) = Dtau(i2);
            i2++;
        }
    }

    tau0 = tauVec;
    tauVec = tauVec + Dtau_k1;
    
    tauVec.elem(arma::find(tauVec < tol && tau0 < tol)).zeros();
    double step = 1.0;
    while (arma::any(tauVec < 0.0)) {
        step *= 0.5;
        tauVec = tau0 + step * Dtau_k1;
        tauVec.elem(arma::find(tauVec < tol && tau0 < tol)).zeros();
    }
    tauVec.elem(arma::find(tauVec < tol)).zeros();
    
}
