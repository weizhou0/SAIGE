#ifndef GLMM_HELPER_H
#define GLMM_HELPER_H

#include <armadillo>

// Structure to hold the result of the ScoreTest_NULL_Model function
struct SA_NULL {
    arma::mat XV;
    arma::mat XVX;
    arma::mat XXVX_inv;
    arma::mat XVX_inv;
    arma::vec S_a;
    arma::mat XVX_inv_XV;
    arma::vec V;
};

// Function to perform the Score Test for the NULL Model
SA_NULL ScoreTest_NULL_Model(const arma::vec& mu, const arma::vec& mu2, const arma::vec& y, const arma::mat& X);

// Function to calculate the AI score
void getAIScore(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, const arma::ivec& fixtauVec,
                const arma::vec& Sigma_iY, const arma::mat& Sigma_iX, const arma::mat& cov, arma::mat& AI, arma::vec& YPAPY, arma::vec& Trace,
                arma::vec& PY, unsigned int k1, int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, bool LOCO,
                const arma::ivec& Ivec_start_indices, const arma::mat& eMat);

// Function to get the trace
void GetTrace(arma::mat Sigma_iX, arma::mat& Xmat, arma::vec& wVec, arma::vec& tauVec, arma::ivec & fixtauVec, arma::mat& cov1, int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, bool LOCO, arma::vec& trace, const arma::uvec & Ivec_start_indices);

// Function to set the number of OpenMP threads
void set_g_omp_num_threads(unsigned int t_omp_num_threads);

// Function to get the product of Imat and Imatt with bVec
arma::vec getprodImatImattbVec(arma::vec & bVec);

// Function to get the coefficients for multiV
void getCoefficients_multiV(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, int maxiterPCG, double tolPCG, bool LOCO, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov, arma::vec& alpha, arma::vec& eta);

// Function to get the coefficients
void Get_Coef(const arma::vec& y, const arma::mat& X, const arma::vec& tau, const arma::vec& offset, const arma::vec& var_weights, 
              GLMResults* fit0ptr, arma::vec& alpha, arma::vec& eta, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov, 
              arma::vec& Y, int maxiterPCG, double tolPCG, int maxiter, bool verbose, bool LOCO);

// Function to update tau for fitglmmaiRPCG_multiV
void fitglmmaiRPCG_multiV_updateTau(const arma::vec& Yvec, const arma::mat& Xmat, const arma::vec& wVec, const arma::vec& tauVec, const arma::ivec& fixtauVec, arma::vec& Sigma_iY, arma::mat& Sigma_iX, arma::mat& cov,
                                    arma::vec& alpha, int nrun, int maxiterPCG, double tolPCG, double tol, double traceCVcutoff, bool LOCO);

// Function to read phenotype file
void readPhenoFile(const std::string& phenoFile, const std::string& phenoCol, const std::vector<std::string>& covarColList, const std::vector<std::string>& categoricCovarCol, const std::vector<std::string>& envCovarCol, const bool& MaleOnly, const bool& FemaleOnly, const std::string& sexCol, const std::string& FemaleCode, const std::string& MaleCode, std::vector<std::string>& sampleIDsinPheno, arma::mat& X, std::vector<std::string>& Xname, arma::mat& eMat, arma::vec& pheno);

// Function to read sample IDs from fam file
void readSampleIDsFromFamFile(const std::string& famFile, std::vector<std::string>& sampleIDsinfamFile);

// Function to read sample IDs from a single column file
void readSampleIDsFromSingleColFile(const std::string& sampleIDFile, std::set<std::string>& sampleIDs, const std::string sampleIDFileName);

// Function to get the PCG of Sigma and Vector
arma::vec getPCG1ofSigmaAndVector(arma::vec& wVec, arma::vec& tauVec, arma::vec& bVec, int maxiterPCG, double tolPCG, bool LOCO, const arma::sp_mat & spSigma, const bool isGRM, const bool isspGRM, const arma::sp_mat& spGRM, const arma::ivec& Ivec_start_indices, const arma::mat& eMat, const arma::vec& REEt_diagVec, const NullGENO::NullGenoClass* ptr_gNULLGENOobj, unsigned int omp_num_threads, const arma::mat& EEt_sqrtEigenMat);

// Function to get the diagonal of Sigma
arma::vec getDiagOfSigma(arma::vec& wVec, arma::vec& tauVec, bool LOCO, const bool isGRM, const bool isspGRM, const arma::sp_mat& spGRM, const arma::ivec& Ivec_start_indices, const arma::mat& eMat, const arma::vec& REEt_diagVec, const arma::mat& EEt_sqrtEigenMat, const NullGENO::NullGenoClass* ptr_gNULLGENOobj, unsigned int omp_num_threads);

// Function to get the cross product
arma::colvec getCrossprod(arma::colvec& bVec, arma::vec& wVec, arma::vec& tauVec, bool LOCO, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices, const bool isGRM, const bool isspGRM, const arma::sp_mat& spGRM, const arma::mat& eMat, const arma::mat& EEt_sqrtEigenMat, const NullGENO::NullGenoClass* ptr_gNULLGENOobj);

// Function to get the cross product of Mat and Kin
arma::vec getCrossprodMatAndKin(arma::colvec& bVec, bool LOCO, const bool isspGRM, const arma::sp_mat& spGRM, const NullGENO::NullGenoClass* ptr_gNULLGENOobj);

// Function to get the product of Imat and bVec
arma::vec getprodImatbVec(arma::vec & bVec, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices);

// Function to get the product of Imat_t and bVec
arma::vec getprodImat_t_bVec(arma::vec & bVec, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices);

// Function to get the cross product of Mat and Kin with eMat
arma::vec getCrossprodMatAndKin_eMat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat, const bool isspGRM, const arma::sp_mat& spGRM, const NullGENO::NullGenoClass* ptr_gNULLGENOobj);

// Function to get the cross product of Mat and I with eMat
arma::vec getCrossprodMatAndI_eMat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat);

// Function to get the cross product of Mat and I with eMat and Imat
arma::vec getCrossprodMatAndI_eMat_Imat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices);

// Function to get the cross product of Mat and Kin with eMat and Imat
arma::vec getCrossprodMatAndKin_eMat_Imat(arma::colvec& bVec, bool LOCO, const arma::mat& EEt_sqrtEigenMat, unsigned int omp_num_threads, const arma::ivec& Ivec_start_indices, const bool isspGRM, const arma::sp_mat& spGRM, const NullGENO::NullGenoClass* ptr_gNULLGENOobj);

#endif // GLMM_HELPER_H