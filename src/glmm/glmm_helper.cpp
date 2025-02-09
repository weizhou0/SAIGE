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

void readPhenoFile(
    const std::string& phenoFile,
    const std::string& phenoCol,
    const std::vector<std::string>& covarColList,
    const std::vector<std::string>& categoricCovarCol,
    const std::vector<std::string>& envCovarCol,
    const bool& MaleOnly,
    const bool& FemaleOnly,
    const std::string& sexCol,
    const std::string& FemaleCode,
    const std::string& MaleCode,
    std::vector<std::string>& sampleIDsinPheno,
    arma::mat& X,
    std::vector<std::string>& Xname,
    arma::mat& eMat,
    arma::vec& pheno
) {
    std::ifstream file(phenoFile);
    if (!file.is_open()) {
        throw std::runtime_error("ERROR! Unable to open pheno file");
    }

    std::string line;
    std::getline(file, line); // Read the header line
    std::istringstream headerStream(line);
    std::vector<std::string> columnNames;
    std::string columnName;
    char delimiter = (line.find('\t') != std::string::npos) ? '\t' : ' ';
    while (std::getline(headerStream, columnName, delimiter)) {
        columnNames.push_back(columnName);
    }

    // Find the index of the phenotype column
    auto phenoColIt = std::find(columnNames.begin(), columnNames.end(), phenoCol);
    if (phenoColIt == columnNames.end()) {
        throw std::runtime_error("ERROR! Phenotype column not found in pheno file");
    }
    size_t phenoColIndex = std::distance(columnNames.begin(), phenoColIt);

    // Find elements that are in covarColList but not in categoricCovarCol
    std::vector<std::string> numericCovarCol;
    std::set_difference(
        covarColList.begin(), covarColList.end(),
        categoricCovarCol.begin(), categoricCovarCol.end(),
        std::back_inserter(numericCovarCol)
    );

    // Find the indices of the covariate columns
    std::vector<size_t> categoricCovarIndices;
    if (!categoricCovarCol.empty()) {
        for (const auto& col : categoricCovarCol) {
            auto it = std::find(columnNames.begin(), columnNames.end(), col);
            if (it == columnNames.end()) {
                throw std::runtime_error("ERROR! Covariate column not found in pheno file: " + col);
            }
            categoricCovarIndices.push_back(std::distance(columnNames.begin(), it));
        }
    }

    // Find the indices of the covariate columns
    std::vector<size_t> numericCovarIndices;
    if (!numericCovarCol.empty()) {
        for (const auto& col : numericCovarCol) {
            auto it = std::find(columnNames.begin(), columnNames.end(), col);
            if (it == columnNames.end()) {
                throw std::runtime_error("ERROR! Covariate column not found in pheno file: " + col);
            }
            numericCovarIndices.push_back(std::distance(columnNames.begin(), it));
        }
    }

    // Find the indices of the env covariate columns
    std::vector<size_t> envCovarIndices;
    if(!envCovarCol.empty()) {
        for (const auto& col : envCovarCol) {
            auto it = std::find(columnNames.begin(), columnNames.end(), col);
            if (it == columnNames.end()) {
                throw std::runtime_error("ERROR! Covariate column not found in pheno file: " + col);
            }
            envCovarIndices.push_back(std::distance(columnNames.begin(), it));
        }
    }

    std::vector<std::vector<double>> covariateData_numeric;
    std::vector<std::vector<double>> covariateData_env;
    std::vector<std::vector<std::string>> covariateData_categorical; 

    std::vector<double> phenotypeValues;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        // Find the index of the sample ID column
        auto sampleIDColIt = std::find(columnNames.begin(), columnNames.end(), "sampleID");
        if (sampleIDColIt == columnNames.end()) {
            throw std::runtime_error("ERROR! Sample ID column not found in pheno file");
        }
        size_t sampleIDColIndex = std::distance(columnNames.begin(), sampleIDColIt);

        std::string id;
        for (size_t i = 0; i <= sampleIDColIndex; ++i) {
            std::getline(iss, id, delimiter);
        }

        // Check if the line contains missing data (NA)
        bool hasMissingData = false;

        std::vector<std::string> categoric_covariates(categoricCovarCol.size());
        std::vector<double> env_covariates(envCovarCol.size());
        std::vector<double> numeric_covariates(numericCovarCol.size());

        // Read categorical covariates
        for (size_t i = 0; i < categoricCovarIndices.size(); ++i) {
            std::string value;
            for (size_t j = 0; j <= categoricCovarIndices[i]; ++j) {
                std::getline(iss, value, delimiter);
            }
            if (value == "NA") {
                hasMissingData = true;
                break;
            }
            categoric_covariates[i] = value;
        }

        // Read numeric covariates
        for (size_t i = 0; i < numericCovarIndices.size(); ++i) {
            std::string value;
            for (size_t j = 0; j <= numericCovarIndices[i]; ++j) {
                std::getline(iss, value, delimiter);
            }
            if (value == "NA") {
                hasMissingData = true;
                break;
            }
            numeric_covariates[i] = std::stod(value);
        }

        // Read environmental covariates
        for (size_t i = 0; i < envCovarIndices.size(); ++i) {
            std::string value;
            for (size_t j = 0; j <= envCovarIndices[i]; ++j) {
                std::getline(iss, value, delimiter);
            }
            if (value == "NA") {
                hasMissingData = true;
                break;
            }
            env_covariates[i] = std::stod(value);
        }

        std::string value;
        for (size_t i = 0; i <= phenoColIndex; ++i) {
            std::getline(iss, value, delimiter);
        }
        if (value == "NA") {
            hasMissingData = true;
        }

        if (FemaleOnly && MaleOnly) {
            throw std::runtime_error("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job");
        } else if (FemaleOnly || MaleOnly) {
            auto sexColIt = std::find(columnNames.begin(), columnNames.end(), sexCol);
            if (sexColIt == std::string::npos) {
                throw std::runtime_error("ERROR! Sex column not found in pheno file");
            }
            size_t sexColIndex = std::distance(columnNames.begin(), sexColIt);

            std::string sexval;
            for (size_t i = 0; i <= sexColIndex; ++i) {
                std::getline(iss, sexval, delimiter);
            }
            if (FemaleOnly) {
                std::cout << "Female-specific model will be fitted. Samples coded as " << FemaleCode << " in the column " << sexCol << " in the phenotype file will be included" << std::endl;
                if (sexval != FemaleCode) {
                    hasMissingData = true;
                }
            } else if (MaleOnly) {
                std::cout << "Male-specific model will be fitted. Samples coded as " << MaleCode << " in the column " << sexCol << " in the phenotype file will be included" << std::endl;
                if (sexval != MaleCode) {
                    hasMissingData = true;
                }
            }
        }

        if (!hasMissingData) {
            sampleIDsinPheno.push_back(id);
            covariateData_categorical.push_back(categoric_covariates);
            covariateData_numeric.push_back(numeric_covariates);
            covariateData_env.push_back(env_covariates);
            phenotypeValues.push_back(std::stod(value));
        }
    }

    file.close();

    // Convert covariateData to arma::mat X
    // First convert the categorical covariates to factors
    std::vector<std::vector<int>> factorCovariates;
    if (!covariateData_categorical.empty()) {
        for (const auto& covariate : covariateData_categorical) {
            std::vector<int> factor = asFactor(covariate);
            size_t numLevels = *std::max_element(factor.begin(), factor.end()) + 1;
            for (size_t level = 1; level < numLevels; ++level) {
                std::vector<int> binaryVector(factor.size(), 0);
                for (size_t i = 0; i < factor.size(); ++i) {
                    if (factor[i] == level) {
                        binaryVector[i] = 1;
                    }
                }
                factorCovariates.push_back(binaryVector);
            }
        }
    }


    // Combine factor covariates and numeric covariates
    std::vector<std::vector<double>> combinedCovariates;
    Xname.push_back("Intercept");

    // Initialize combinedCovariates with the size of covariateData_numeric
    if (!covariateData_numeric.empty()) {
        combinedCovariates.resize(covariateData_numeric.size());
        for (size_t i = 0; i < covariateData_numeric.size(); ++i) {
            combinedCovariates[i] = covariateData_numeric[i];
            Xname.push_back(numeric_covariates[i]);
        }
    }

    // Add factor covariates to combinedCovariates if not empty
    std::string factorCovariateName;
    if (!factorCovariates.empty()) {
        for (size_t i = 0; i < factorCovariates.size(); ++i) {
            for (size_t j = 0; j < factorCovariates[i].size(); ++j) {
                combinedCovariates[j].push_back(static_cast<double>(factorCovariates[i][j]));
                factorCovariateName = categoricCovarCol[i] + std::to_string(j);
                Xname.push_back(factorCovariateName);

            }
        }
    }

    // Convert combinedCovariates to arma::mat X
    if (!combinedCovariates.empty()) {
        X.set_size(combinedCovariates.size(), combinedCovariates[0].size() + 1);
        for (size_t i = 0; i < combinedCovariates.size(); ++i) {
            X(i, 0) = 1.0; // Intercept column with 1s
            for (size_t j = 0; j < combinedCovariates[i].size(); ++j) {
                X(i, j + 1) = combinedCovariates[i][j];
            }
        }
    } else {
        X.set_size(sampleIDsinPheno.size(), 1);
        X.fill(1.0); // Only intercept column with 1s
    }

    // Convert env_covariates to arma::mat eMat if not empty
    if (!covariateData_env.empty()) {
        std::vector<std::vector<double>> combinedEnvCovariates;
        for (size_t i = 0; i < covariateData_env.size(); ++i) {
            combinedEnvCovariates.push_back(covariateData_env[i]);
        }

        // Convert combinedEnvCovariates to arma::mat eMat
        eMat.set_size(combinedEnvCovariates.size(), combinedEnvCovariates[0].size());
        for (size_t i = 0; i < combinedEnvCovariates.size(); ++i) {
            for (size_t j = 0; j < combinedEnvCovariates[i].size(); ++j) {
                eMat(i, j) = combinedEnvCovariates[i][j];
            }
        }
    }

    // Convert phenotypeValues to arma::vec pheno
    pheno.set_size(phenotypeValues.size());
    for (size_t i = 0; i < phenotypeValues.size(); ++i) {
        pheno(i) = phenotypeValues[i];
    }
    std::cout << sampleIDsinPheno.size() << " samples found in pheno file" << std::endl;
}

void readSampleIDsFromFamFile(const std::string& famFile, std::vector<std::string>& sampleIDsinfamFile) {
    if (!famFile.empty()) {
        std::ifstream famfile(famFile);
        if (!famfile.is_open()) {
            throw std::runtime_error("ERROR! Unable to open fam file");
        }
        std::string line;
        while (std::getline(famfile, line)) {
            std::istringstream iss(line);
            std::string id;
            iss >> id >> id; // Skip the first column and read the second column
            sampleIDsinfamFile.push_back(id);
        }
        famfile.close();
        std::cout << sampleIDsinfamFile.size() << " samples found in fam file" << std::endl;
    }
}


void readSampleIDsFromSingleColFile(const std::string& sampleIDFile, std::set<std::string>& sampleIDs, const std::string sampleIDFileName) {
    if (!sampleIDFile.empty()) {
        std::ifstream file(sampleIDFile);
        if (!file.is_open()) {
            throw std::runtime_error("ERROR! Unable to open sample ID file", sampleIDFile);
        }
        std::string line;
        while (std::getline(file, line)) {
            sampleIDs.insert(line);
        }
        file.close();
        std::cout << sampleIDs.size() << " samples found in sample ID file " << sampleIDFileName << std::endl;
    }
}



  
    arma::vec getPCG1ofSigmaAndVector(
        arma::vec& wVec,  
        arma::vec& tauVec, 
        arma::vec& bVec, 
        int maxiterPCG, 
        double tolPCG, 
        bool LOCO, 
        const arma::sp_mat & spSigma, 
        const bool isGRM,
        const bool isspGRM,
        const arma::sp_mat& spGRM,
        const arma::ivec& Ivec_start_indices,
        const arma::mat& eMat,
        const arma::vec& REEt_diagVec,
        const NullGENO::NullGenoClass* ptr_gNULLGENOobj,
        unsigned int omp_num_threads, 
        const arma::mat& EEt_sqrtEigenMat
 ) {
        int Nnomissing = wVec.n_elem;
        arma::vec xVec(Nnomissing);
        xVec.zeros();
        std::cout << "g_isStoreSigma " << g_isStoreSigma << std::endl;
        if(spSigma.n_rows > 0){
            xVec = arma::spsolve(spSigma, bVec);
        }else{
            arma::vec rVec = bVec;
            arma::vec r1Vec;
            arma::vec crossProdVec(Nnomissing);
            arma::vec zVec(Nnomissing);
            arma::vec minvVec(Nnomissing);
            minvVec = 1/getDiagOfSigma(wVec, tauVec, LOCO, isGRM, isspGRM, spSigma, Ivec_start_indices, eMat, REEt_diagVec, EEt_sqrtEigenMat, ptr_gNULLGENOobj, omp_num_threads);

            zVec = minvVec % rVec;

            double sumr2 = sum(rVec % rVec);
            arma::vec z1Vec(Nnomissing);
            arma::vec pVec = zVec;

            int iter = 0;
            arma::colvec ApVec;
            while (sumr2 > tolPCG && iter < maxiterPCG) {
                    iter = iter + 1;
                    ApVec = getCrossprod(pVec, wVec, tauVec, LOCO, omp_num_threads, Ivec_start_indices, isGRM, isspGRM, spGRM, eMat, EEt_sqrtEigenMat, ptr_gNULLGENOobj);

                    arma::vec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
                    double a = preA(0);
                    xVec = xVec + a * pVec;
                    r1Vec = rVec - a * ApVec;
                    z1Vec = minvVec % r1Vec;

                    arma::vec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                    double bet = Prebet(0);
                    pVec = z1Vec + bet * pVec;
                    zVec = z1Vec;
                    rVec = r1Vec;
                    sumr2 = sum(rVec % rVec);
            }

            if (iter >= maxiterPCG){
                    cout << "pcg did not converge. You may increase maxiter number." << endl;
            }
            cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
        }
        return xVec;
    }


arma::vec getDiagOfSigma(
    arma::vec& wVec,
    arma::vec& tauVec,
    bool LOCO,
    const bool isGRM,
    const bool isspGRM,
    const arma::sp_mat& spGRM,
    const arma::ivec& Ivec_start_indices,
    const arma::mat& eMat,
    const arma::vec& REEt_diagVec,
    const arma::mat& EEt_sqrtEigenMat,
    const NullGENO::NullGenoClass* ptr_gNULLGENOobj,
    unsigned int omp_num_threads
) {
    int Nnomissing = wVec.n_elem;
    arma::vec diagVec(Nnomissing);
    arma::vec diagVec0(Nnomissing);
    diagVec0.zeros();
    arma::sp_vec diagVecV0;
    arma::vec diagVecG, diagVecV, diagVecG_I, diagVecV_I;
    unsigned int tauind = 0;

    if (Ivec_start_indices.n_elem == 0) {
        if (isGRM) {
            if (!(ptr_gNULLGENOobj->setKinDiagtoOne)) {
                if (!isspGRM) {
                    if (!LOCO) {
                        int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
                        diagVec = tauVec(1) * (*ptr_gNULLGENOobj->Get_Diagof_StdGeno()) / MminMAF + tauVec(0) / wVec;
                    } else {
                        diagVec = tauVec(1) * (*ptr_gNULLGENOobj->Get_Diagof_StdGeno_LOCO());
                        int Msub_MAFge_minMAFtoConstructGRM_in_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_in();
                        int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
                        diagVec = diagVec / (Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b) + tauVec(0) / wVec;
                    }
                    tauind = tauind + 2;
                } else {
                    diagVec = tauVec(0) / wVec;
                    tauind = tauind + 1;
                    diagVecG = spGRM.diag();
                    diagVec = diagVec + tauVec(tauind) * diagVecG;
                    tauind = tauind + 1;
                }
            } else {
                diagVec = tauVec(1) + tauVec(0) / wVec;
                tauind = tauind + 2;
            }
        } else {
            diagVec = tauVec(0) / wVec;
            tauind = tauind + 1;
        }
    } else {
        diagVec = tauVec(0) / wVec;
        tauind = tauind + 1;
        if (isGRM && isspGRM) {
            diagVecG = spGRM.diag();
            //diagVecG_I = diagVecG.elem(g_I_longl_vec);
            diagVecG_I = getprodImatbVec(diagVecG, omp_num_threads, Ivec_start_indices);

            diagVec = diagVec + tauVec(tauind) * diagVecG_I;
            tauind = tauind + 1;
        } else {
            diagVec = diagVec + tauVec(tauind);
            tauind = tauind + 1;
        }
    }

    if (Ivec_start_indices.n_elem == 0) {
        if (isGRM) {
            if (!(ptr_gNULLGENOobj->setKinDiagtoOne)) {
                if (!isspGRM) {
                    if (!LOCO) {
                        int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
                        diagVec0.zeros();
                        for (unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
                            diagVec0 = diagVec0 + ((*ptr_gNULLGENOobj->Get_Diagof_StdGeno()) / MminMAF) % EEt_sqrtEigenMat.col(i);
                        }
                    } else {
                        int Msub_MAFge_minMAFtoConstructGRM_in_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_in();
                        int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
                        diagVec0.zeros();
                        for (unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
                            diagVec0 = diagVec0 + ((*ptr_gNULLGENOobj->Get_Diagof_StdGeno_LOCO()) / (Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b)) % EEt_sqrtEigenMat.col(i);
                        }
                    }
                } else {
                    diagVec0.zeros();
                    for (unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
                        diagVec0 = diagVec0 + (g_spGRM.diag()) % (EEt_sqrtEigenMat.col(i));
                    }
                }
                diagVec = diagVec + tauVec(tauind) * diagVec0;
            } else {
                diagVec = tauVec(1) + tauVec(0) / wVec;
                tauind = tauind + 2;
                diagVec0.zeros();
                for (unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
                    diagVec0 = diagVec0 + (EEt_sqrtEigenMat.col(i));
                }
                diagVec = diagVec + tauVec(tauind) * diagVec0;
            }
        }
        tauind = tauind + 1;
    } else {
        if (isGRM && isspGRM) {
            diagVecG = spGRM.diag();
            //diagVecG_I = diagVecG.elem(g_I_longl_vec);
            diagVecG_I = getprodImatbVec(diagVecG, omp_num_threads, Ivec_start_indices);

            diagVec0.zeros();
            for (unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
                diagVec0 = diagVec0 + diagVecG_I % (EEt_sqrtEigenMat.col(i));
            }
            diagVec = diagVec + tauVec(tauind) * diagVec0;
        } else {
            diagVec = diagVec + tauVec(tauind) * REEt_diagVec;
        }
        tauind = tauind + 1;
    }

    for (unsigned int i = 0; i < Nnomissing; i++) {
        if (diagVec(i) < 1e-4) {
            diagVec(i) = 1e-4;
        }
    }

    return diagVec;
}

arma::colvec getCrossprod(
    arma::colvec& bVec,
    arma::vec& wVec,
    arma::vec& tauVec,
    bool LOCO,
    unsigned int omp_num_threads, 
    const arma::ivec& Ivec_start_indices,
    const bool isGRM,
    const bool isspGRM,
    const arma::sp_mat& spGRM,
    const arma::mat& eMat,
    const arma::mat& EEt_sqrtEigenMat, 
    const NullGENO::NullGenoClass* ptr_gNULLGENOobj
) {
    arma::colvec crossProdVec;
    arma::vec crossProd1, GRM_I_bvec, Ibvec, V_I_bvec, V_T_bvec;

    unsigned int tau_ind = 0;
    crossProdVec = tauVec(0) * (bVec % (1 / wVec));
    tau_ind = tau_ind + 1;


    if (Ivec_start_indices.n_rows == 0) { // it must have specified GRM
        if (isGRM) {
            crossProd1 = getCrossprodMatAndKin(bVec, LOCO, isspGRM, spGRM, ptr_gNULLGENOobj);
            crossProdVec = crossProdVec + tauVec(1) * crossProd1;
        }
    } else {
        if (isGRM) {
            Ibvec = g_I_longl_mat_t * bVec; //need to update
            GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO, isspGRM, spGRM, ptr_gNULLGENOobj);
            crossProd1 = g_I_longl_mat * GRM_I_bvec;
            crossProdVec = crossProdVec + tauVec(1) * crossProd1;
        } else {
            crossProd1 = getprodImatImattbVec(bVec, omp_num_threads, Ivec_start_indices);
            crossProdVec = crossProdVec + tauVec(1) * crossProd1;
        }
    }


    tau_ind = tau_ind + 1;
if(eMat.n_rows > 0){
    if (Ivec_start_indices.n_rows == 0) {
        if (isGRM) {
            crossProd1 = getCrossprodMatAndKin_eMat(bVec, LOCO, EEt_sqrtEigenMat, isspGRM, spGRM, ptr_gNULLGENOobj);
        } else {
            crossProd1 = getCrossprodMatAndI_eMat(bVec, LOCO, EEt_sqrtEigenMat);
        }
    } else {
        if (isGRM) {
            crossProd1 = getCrossprodMatAndKin_eMat_Imat(bVec,  LOCO, EEt_sqrtEigenMat, omp_num_threads, Ivec_start_indices, isspGRM, spGRM, ptr_gNULLGENOobj);
           
        } else {
            crossProd1 = getCrossprodMatAndI_eMat_Imat(bVec, LOCO, EEt_sqrtEigenMat, omp_num_threads, Ivec_start_indices);
        }
    }
    crossProdVec = crossProdVec + tauVec(tau_ind) * crossProd1;
    tau_ind = tau_ind + 1;
}
    return crossProdVec;
}

arma::vec getCrossprodMatAndKin(arma::colvec& bVec, bool LOCO,  const bool isspGRM, const arma::sp_mat& spGRM,  const NullGENO::NullGenoClass* ptr_gNULLGENOobj){
    arma::vec crossProdVec;
    if(isspGRM){
        crossProdVec = spGRM * bVec;
    }else{
        if(!LOCO){
          crossProdVec = parallelCrossProd(bVec, ptr_gNULLGENOobj);
        }else{
          crossProdVec = parallelCrossProd_LOCO(bVec, ptr_gNULLGENOobj);
        }
    }
    return(crossProdVec);
}

arma::vec parallelCrossProd(arma::colvec & bVec, const NullGENO::NullGenoClass* ptr_gNULLGENOobj) {
          int Msub_mafge1perc = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
          CorssProd CorssProd(bVec);
  
          parallelReduce(0, Msub_mafge1perc, CorssProd);

          return CorssProd.m_bout/(CorssProd.Msub_mafge1perc);
  }

arma::vec parallelCrossProd_LOCO(arma::colvec & bVec, const NullGENO::NullGenoClass* ptr_gNULLGENOobj) {

    int numberMarker_full = 0;
    arma::vec outvec = parallelCrossProd_full(bVec, numberMarker_full);

    CorssProd CorssProd(bVec);
    int startIndex = ptr_gNULLGENOobj->getStartIndex();
    int endIndex = ptr_gNULLGENOobj->getEndIndex();

    parallelReduce(startIndex, endIndex + 1, CorssProd);

    outvec = outvec - CorssProd.m_bout;

    int markerNum = numberMarker_full - CorssProd.Msub_mafge1perc;
    return outvec / markerNum;
}

struct CorssProd_usingSubMarker : public Worker
{
    // source vectors
    arma::colvec & m_bVec;
    unsigned int m_N;
    unsigned int m_M_Submarker;
    unsigned int m_M;
    arma::ivec subMarkerIndex ;

    // product that I have accumulated
    arma::vec m_bout;


    // constructors
    CorssProd_usingSubMarker(arma::colvec & y)
        : m_bVec(y) {

        //m_Msub = ptr_gNULLGENOobj->getMsub();
        subMarkerIndex = getSubMarkerIndex();
        m_M_Submarker = subMarkerIndex.n_elem;
        m_N = ptr_gNULLGENOobj->getNnomissing();
        m_bout.zeros(m_N);
    }
    CorssProd_usingSubMarker(const CorssProd_usingSubMarker& CorssProd_usingSubMarker, Split)
        : m_bVec(CorssProd_usingSubMarker.m_bVec)
    {

        m_N = CorssProd_usingSubMarker.m_N;
        //m_M = CorssProd_usingSubMarker.m_M;
        m_M_Submarker = CorssProd_usingSubMarker.m_M_Submarker;
        subMarkerIndex = CorssProd_usingSubMarker.subMarkerIndex;
        m_bout.zeros(m_N);

    }

       // process just the elements of the range I've been asked to
    void operator()(std::size_t begin, std::size_t end) {
        arma::vec vec;
        double val1;
        int j;
        for(unsigned int i = begin; i < end; i++){
            j = subMarkerIndex[i];
//                      std::cout << "j: " << j << std::endl;
            ptr_gNULLGENOobj->Get_OneSNP_StdGeno(j, &vec);
            val1 = dot(vec,  m_bVec);
            m_bout += val1 * (vec);
        }
    }

    // join my value with that of another InnerProduct
    void join(const  CorssProd_usingSubMarker & rhs) {
    m_bout += rhs.m_bout;
    }
};

struct CorssProd : public Worker
{
    // source vectors
    arma::colvec & m_bVec;
    unsigned int m_N;
    unsigned int m_M;

    // product that I have accumulated
    arma::vec m_bout;
    int Msub_mafge1perc;

    // constructors
    CorssProd(arma::colvec & y)
        : m_bVec(y) {

        m_M = ptr_gNULLGENOobj->getM();
        m_N = ptr_gNULLGENOobj->getNnomissing();
        m_bout.zeros(m_N);
        Msub_mafge1perc=0;
        //ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
    }
    CorssProd(const CorssProd& CorssProd, Split)
        : m_bVec(CorssProd.m_bVec)
    {

        m_N = CorssProd.m_N;
        m_M = CorssProd.m_M;
        m_bout.zeros(m_N);
        Msub_mafge1perc=0;
        //CorssProd.Msub_mafge1perc;

    }
    // process just the elements of the range I've been asked to
    void operator()(std::size_t begin, std::size_t end) {
        arma::vec vec;
        for(unsigned int i = begin; i < end; i++){
            //if(ptr_gNULLGENOobj->alleleFreqVec[i] >= minMAFtoConstructGRM && ptr_gNULLGENOobj->alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                ptr_gNULLGENOobj->Get_OneSNP_StdGeno(i, &vec);
                double val1 = dot(vec,  m_bVec);
                m_bout += val1 * (vec) ;
                Msub_mafge1perc += 1;

        }
    }

    // join my value with that of another InnerProduct
    void join(const CorssProd & rhs) {
        m_bout += rhs.m_bout;
        Msub_mafge1perc += rhs.Msub_mafge1perc;
    }
};

arma::vec getprodImatImattbVec(arma::vec & bVec, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices){
    omp_set_num_threads(omp_num_threads);
    auto n = Ivec_start_indices.n_elem - 1;
    arma::vec resultVec(bVec.n_elem, arma::fill::zeros);
    #pragma omp parallel
    {
            auto thread_idx = omp_get_thread_num();
            for(int j = thread_idx; j < n; j += g_omp_num_threads) {
                    double sum = 0;
                    size_t start = Ivec_start_indices[j];
                    size_t end = Ivec_start_indices[j + 1];
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

arma::vec getprodImatbVec(arma::vec & bVec, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices){
    omp_set_num_threads(omp_num_threads);
    int N = Ivec_start_indices[Ivec_start_indices.n_elem]; // -1?
    auto n = Ivec_start_indices.n_elem - 1;

    arma::vec resultVec(N, arma::fill::zeros);
    #pragma omp parallel
    {
            auto thread_idx = omp_get_thread_num();
            for(int j = thread_idx; j < n; j += g_omp_num_threads) {
                    double sum = 0;
                    size_t start = Ivec_start_indices[j];
                    size_t end = Ivec_start_indices[j + 1];
                    
                    for(size_t k = start; k < end; k++) {
                        resultVec[k] = bVec[j];
                    }

            }
    }
    return(resultVec);
}

arma::vec getprodImat_t_bVec(arma::vec & bVec, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices){
    omp_set_num_threads(omp_num_threads);
    auto n = Ivec_start_indices.n_elem - 1;
    arma::vec resultVec(n, arma::fill::zeros);
    #pragma omp parallel
    {
            auto thread_idx = omp_get_thread_num();
            for(int j = thread_idx; j < n; j += g_omp_num_threads) {
                    double sum = 0;
                    size_t start = Ivec_start_indices[j];
                    size_t end = Ivec_start_indices[j + 1];
                    for(size_t k = start; k < end; k++) {
                        resultVec[j] += bVec[k];
                    }
                    //for(size_t k = start; k < end; k++) {
                    //        resultVec[k] = sum;
                    //}
            }
    }
    return(resultVec);
}




arma::vec getCrossprodMatAndKin_eMat(arma::colvec& bVec, bool LOCO,  const arma::mat& EEt_sqrtEigenMat, const bool isspGRM, const arma::sp_mat& spGRM,  const NullGENO::NullGenoClass* ptr_gNULLGENOobj){
    arma::colvec bVec_1, crossProdVecGRM_1, crossProdVecGRM, crossProdVec;
    crossProdVec.zeros(bVec.n_elem);
    for(unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++){
            bVec_1 = EEt_sqrtEigenMat.col(i) % bVec;
            crossProdVecGRM = getCrossprodMatAndKin(bVec_1, LOCO, isspGRM, spGRM, ptr_gNULLGENOobj);
            crossProdVecGRM_1 = EEt_sqrtEigenMat.col(i) % crossProdVecGRM;
            crossProdVec = crossProdVec + crossProdVecGRM_1;
            crossProdVecGRM.zeros();
            crossProdVecGRM_1.zeros();
    }
    return(crossProdVec);
}

arma::vec getCrossprodMatAndI_eMat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat){
    arma::colvec bVec_1, crossProdVecGRM_1, crossProdVecGRM, crossProdVec;
    crossProdVec.zeros(bVec.n_elem);
    for(unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++){
            arma::colvec sqrtEigenVec = EEt_sqrtEigenMat.col(i);
            bVec_1 = sqrtEigenVec % bVec;
            crossProdVecGRM = bVec_1;
            crossProdVecGRM_1 = sqrtEigenVec % crossProdVecGRM;
            crossProdVec = crossProdVec + crossProdVecGRM_1;
    }
    return(crossProdVec);
}


arma::vec getCrossprodMatAndI_eMat_Imat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices) {
    arma::colvec crossProdVec(bVec.n_elem, arma::fill::zeros);
    omp_set_num_threads(omp_num_threads);

    auto n = Ivec_start_indices.n_elem - 1;
#pragma omp parallel
{
    auto thread_idx = omp_get_thread_num();
    auto g_omp_num_threads = omp_get_num_threads();
    for(unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++) {
        for(int j = thread_idx; j < n; j += g_omp_num_threads) {
            double sum = 0;
            size_t start = g_I_start_indices[j];
            size_t end = g_I_start_indices[j + 1];
            for(size_t k = start; k < end; k++) {
                sum += EEt_sqrtEigenMat.at(k, i) * bVec[k];
            }
            for(size_t k = start; k < end; k++) {
                crossProdVec[k] += EEt_sqrtEigenMat.at(k, i) * sum;
            }
        }
    }
}
    return crossProdVec;
}

arma::vec getCrossprodMatAndKin_eMat_Imat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices, const bool isspGRM, const arma::sp_mat& spGRM,  const NullGENO::NullGenoClass* ptr_gNULLGENOobj){
    arma::colvec bVec_1, crossProdVecGRM_1, IbVec, GRM_I_bvec, crossProdVecGRM, crossProdVec;
    crossProdVec.zeros(bVec.n_elem);
    for(unsigned int i = 0; i < EEt_sqrtEigenMat.n_cols; i++){
            bVec_1 = EEt_sqrtEigenMat.col(i) % bVec;
            IbVec = getprodImat_t_bVec(bVec_1,omp_num_threads, Ivec_start_indices);
            GRM_I_bvec = getCrossprodMatAndKin(IbVec, LOCO, , isspGRM, spGRM, ptr_gNULLGENOobj);
            crossProdVecGRM = getprodImatbVec(GRM_I_bvec, omp_num_threads, Ivec_start_indices);

            crossProdVecGRM_1 = EEt_sqrtEigenMat.col(i) % crossProdVecGRM;
            crossProdVec = crossProdVec + crossProdVecGRM_1;
            crossProdVecGRM.zeros();
            crossProdVecGRM_1.zeros();
    }
    return(crossProdVec);
}
