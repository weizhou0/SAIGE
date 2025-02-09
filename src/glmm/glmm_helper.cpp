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


void getAIScore(const arma::fvec& Yvec, const arma::fmat& Xmat, const arma::fvec& wVec, const arma::fvec& tauVec, const arma::ivec& fixtauVec,
                const arma::fvec& Sigma_iY, const arma::fmat& Sigma_iX, const arma::fmat& cov, arma::fmat& AI, arma::fvec& YPAPY, arma::fvec& Trace,
                arma::fvec& PY, unsigned int k1, int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO,
                const arma::ivec& Ivec_start_indices, const arma::fmat& eMat) {

    fixtauVec.print("fixtauVec");
    int q2 = arma::sum(fixtauVec == 0);
    arma::uvec idxtau = arma::find(fixtauVec == 0);
    arma::fvec tau0;
    AI.set_size(k1, k1);
    YPAPY.zeros(k1);
    Trace.zeros(k1);
    arma::fmat Sigma_iXt = Sigma_iX.t();
    arma::fmat Xmatt = Xmat.t();

    PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));

    unsigned int n = PY.n_elem;
    arma::fvec PAPY_1, PAPY, APY;
    arma::fmat APYmat(n, k1);

    arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;

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

    Trace = GetTrace_multiV_eMat(Sigma_iX, Xmat, wVec, tauVec, fixtauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);
}
