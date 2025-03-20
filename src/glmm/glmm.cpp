#include "glmm.h"
#include <iostream>
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
//#include "stats.hpp"
#include <unordered_map>
#include "utils/parameters.h"
#include "glmm_helper.h"
#include "glmmresults.h"
#include "utils/utils.h"

void fitNULLGLMM(const Parameters& params) {
    std::vector<std::string> filePaths = {
        params.bedFile,
        params.bimFile,
        params.famFile,
        params.phenoFile,
        params.sparseGRMSampleIDFile,
        params.sparseGRMS,
        params.SampleIDIncludeFile
    };
    checkFiles(filePaths);

    // determine the samples to include
    std::vector<std::string> sampleIDsinPheno, sampleIDsinfamFile, sampleIDsinsparseGRMFile, sampleIDsinSampleIncludeFile, Xname;
    arma::mat X, eMat;
    arma::vec pheno;
    readPhenoFile(params.phenoFile, params.phenoCol, params.covarColList, params.categoricCovarCol, params.envCovarCol, params.MaleOnly, params.FemaleOnly, params.sexCol, params.FemaleCode, params.MaleCode, sampleIDsinPheno, X, Xname, eMat, pheno);
    readSampleIDsFromFamFile(params.famFile, sampleIDsinfamFile);
    readSampleIDsFromSingleColFile(params.sparseGRMSampleIDFile, sampleIDsinsparseGRMFile, "sparseGRMFile");
    readSampleIDsFromSingleColFile(params.SampleIDIncludeFile, sampleIDsinSampleIncludeFile, "SampleIDIncludeFile");


    // Check for duplicates in sampleIDsinPheno and obtain unique sample IDs while maintaining order
    bool duplicatesFound = false;
    std::vector<std::string> uniqueSampleIDsinPheno;
    std::set<std::string> seen;
    for (const auto& id : sampleIDsinPheno) {
        if (seen.find(id) == seen.end()) {
            uniqueSampleIDsinPheno.push_back(id);
            seen.insert(id);
        } else {
            duplicatesFound = true;
        }
    }

    // Find the samples in all the non-empty vectors
    std::set<std::string> commonSamples;
    std::vector<std::vector<std::string>> sampleVectors = {uniqueSampleIDsinPheno, sampleIDsinfamFile, sampleIDsinsparseGRMFile, sampleIDsinSampleIncludeFile};

    // Initialize commonSamples with the first non-empty vector
    for (const auto& vec : sampleVectors) {
        if (!vec.empty()) {
            commonSamples.insert(vec.begin(), vec.end());
            break;
        }
    }

    // Intersect with other non-empty vectors
    for (const auto& vec : sampleVectors) {
        if (!vec.empty()) {
            std::set<std::string> tempSet;
            std::set_intersection(commonSamples.begin(), commonSamples.end(), vec.begin(), vec.end(), std::inserter(tempSet, tempSet.begin()));
            commonSamples = std::move(tempSet);
        }
    }

    // Convert set to vector
    std::vector<std::string> overlappedSamples;
    overlappedSamples.assign(commonSamples.begin(), commonSamples.end());
    std::cout << overlappedSamples.size() << " samples with non-missing data were found" << std::endl;

    // Print if duplicates were found
    if (duplicatesFound) {
        std::cout << "Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements" << std::endl;
    }

    // Filter the data based on the overlapped samples
    arma::uvec sampleIndices = findSampleIndices(sampleIDsinPheno, overlappedSamples);
    X = X.rows(sampleIndices);
    eMat = eMat.rows(sampleIndices);
    pheno = pheno.elem(sampleIndices);
    std::vector<std::string> sampleIDsinPhenonew = subsetVector(sampleIDsinPheno, sampleIndices);

    // Check if the trait type is quantitative and inverse normalize the phenotype
    if (traitType == "quantitative" && params.invNormalize) {
        std::cout << "Perform the inverse normalization for " << params.phenoCol << std::endl;

        arma::vec ranks = arma::conv_to<arma::vec>::from(arma::sort_index(pheno) + 1);
        arma::vec invPheno = arma::norminv((ranks - 0.5) / ranks.n_elem);

        for (size_t i = 0; i < pheno.n_elem; ++i) {
            pheno(i) = invPheno(i);
        }
    }

    if ((traitType == "binary" || traitType == "tte" )& (length(covarColList) > 0)) {
        checkPerfectSep(X1, pheno, Xname, minCovariateCount);
    }   

    // Call the Covariate_Transform function
    CovariateTransformResult result = Covariate_Transform(X1, pheno, Xname);

    // Subset the sparse GRM file
    if (!params.sparseGRMFile.empty()) {
        arma::uvec sampleIndicesSparseGRM = findSampleIndices(sampleIDsinsparseGRMFile, overlappedSamples);
        arma::sp_mat sparseGRM = readSparseGRM(params.sparseGRMFile, sampleIndicesSparseGRM);
    }

    arma::uvec sampleIndicesinFam = findSampleIndices(sampleIDsinfamFile, overlappedSamples);

    NullGENO::NullGenoClass gNULLGENOobj;
    NullGENO::NullGenoClass* ptr_gNULLGENOobj = nullptr;

    if (!sampleIndicesinFam.is_empty()) {
        std::vector<bool> indicatorGenoSamplesWithPheno(sampleIDsinfamFile.size(), false);
        for (size_t i = 0; i < sampleIndicesinFam.n_elem; ++i) {
            indicatorGenoSamplesWithPheno[sampleIndicesinFam[i]] = true;
        }
        gNULLGENOobj.setGenoObj(params.bedFile, params.bimFile, params.famFile, sampleIndicesinFam, indicatorGenoSamplesWithPheno, params.memoryChunk, params.isDiagofKinSetAsOne);
        ptr_gNULLGENOobj = &gNULLGENOobj;
    }

    // Check for duplicates in overlappedSamples
    std::set<std::string> uniqueSamples;
    bool hasDuplicates = false;
    for (const auto& sample : overlappedSamples) {
        if (!uniqueSamples.insert(sample).second) {
            hasDuplicates = true;
            break;
        }
    }
    arma::uvec Ivec_start_indices;
    arma::vec REEt_diagVec, EEt_sqrtEigenMat, EEt_eigenvalVec. EEt_diagVec;
    arma::mat EEt_eigenvecMat;
    if(eMat.n_rows != 0){
        arma::mat U;
        arma::vec S;
        arma::mat V;
        arma::svd_econ(U, S, V, eMat, "left");

        EEt_eigenvalVec = arma::square(S);
        EEt_eigenvalVec.print("EEt_eigenvalVec");
        EEt_eigenvecMat = U;
        EEt_diagVec.set_size(eMat.n_rows);
        EEt_diagVec = arma::sum(arma::square(eMat), 1);
        EEt_sqrtEigenMat.set_size(arma::size(EEt_eigenvecMat));
        for(unsigned int i = 0; i < EEt_eigenvalVec.n_elem; i++){
                arma::fcolvec sqrtEigenVec = std::sqrt(EEt_eigenvalVec(i)) * (EEt_eigenvalVec.col(i));
                EEt_sqrtEigenMat.col(i) = sqrtEigenVec;
        }

        REEt_diagVec.set_size(eMat.n_rows);
        REEt_diagVec.zeros();
        for(unsigned int i = 0; i < EEt_eigenvalVec.n_elem; i++){
                 arma::vec eigenvec = EEt_eigenvecMat.col(i);
                 double eigenval = g_EEt_eigenvalVec(i);
                 REEt_diagVec = REEt_diagVec + eigenval*arma::square(eigenvec);
        }

    }
    if (hasDuplicates) {
        std::cout << "Duplicated IDs were found " << std::endl;
        std::cout << overlappedSamples.size() << " observations will be used for analysis" << std::endl;
        set_Ivec_start_indices(sampleIDsinPhenonew, Ivec_start_indices);
    }


    GLMMResults modglmm = glmmkin_ai_PCG(
        ptr_gNULLGENOobj,
        X,
        pheno,
        traitType,
        offset,
        var_weights,
        sampleIDsinPhenonew,
        isCovariateOffset,
        tau,
        fixtau,
        params.isCovariateOffset,
        params.maxiter,
        params.tol,
        params.verbose,
        params.nrun,
        params.tolPCG,
        params.maxiterPCG,
        result.subPheno,
        result.indicatorGenoSamplesWithPheno,
        result.out_transform,
        result.tauInit,
        params.memoryChunk,
        params.LOCO,
        result.chromosomeStartIndexVec,
        result.chromosomeEndIndexVec,
        params.traceCVcutoff,
        result.isCovariateTransform,
        result.isDiagofKinSetAsOne,
        params.isLowMemLOCO,
        result.covarianceIdxMat,
        params.isStoreSigma,
        params.useSparseGRMtoFitNULL,
        params.useGRMtoFitNULL,
        params.isSparseGRMIdentity,
        arma::vec & REEt_diagVec
    );

}


// Function to fit the null generalized linear mixed model
GLMMResults glmmkin_ai_PCG(
    const NullGENO::NullGenoClass* ptr_gNULLGENOobj,
    const arma::mat& X,
    const arma::vec& y,
    const std::string& family,
    const arma::vec& offset,
    const arma::vec& var_weights,
    const std::vector<std::string>& sampleIDsinPheno,
    bool isCovariateOffset,
    arma::vec tau,
    arma::vec fixtau,
    int maxiter,
    double tol,
    bool verbose,
    int nrun,
    double tolPCG,
    int maxiterPCG,
    const arma::mat& subPheno,
    const arma::vec& indicatorGenoSamplesWithPheno,
    const arma::mat& out_transform,
    const arma::vec& tauInit,
    double memoryChunk,
    bool LOCO,
    const arma::vec& chromosomeStartIndexVec,
    const arma::vec& chromosomeEndIndexVec,
    double traceCVcutoff,
    bool isCovariateTransform,
    bool isDiagofKinSetAsOne,
    bool isLowMemLOCO,
    const arma::mat& covarianceIdxMat,
    bool isStoreSigma,
    bool useSparseGRMtoFitNULL,
    bool useGRMtoFitNULL,
    bool isSparseGRMIdentity
) {

    GLMMResults glmmResult;
    int n = y.n_elem;
    if (offset.is_empty()) {
        offset = arma::zeros(n);
    }
    if (var_weights.is_empty()) {
        var_weights = arma::ones(n);
    }

    GLMResults fit0 = fitGLM(X, y, traitType);
    arma::vec eta = fit0.eta
    arma::vec mu = fit0.mu
    arma::vec mu_eta = fit0.mu_eta(eta); // family["mu.eta"](eta)
    arma::vec Y = eta - offset + (y - mu) / mu_eta;
    arma::vec sqrtW = mu_eta / arma::sqrt(1 / var_weights % fit0.variance(mu)); // family["variance"](mu)
    arma::vec W = arma::square(sqrtW);
    arma::vec alpha = fit0.beta;
    arma::vec eta = eta;


    tau.zeros();

    if (fit0.family == "poisson" || fit0.family == "binomial") { 
        tau[0] = 1;
        fixtau[0] = 1;
        tauInit[0] = 1;
        arma::uvec idxtau = arma::find(fixtau == 0);
        if (arma::accu(tauInit(idxtau)) == 0) {
            tau(idxtau).fill(0.1);
        } else {
            tau(idxtau) = tauInit(idxtau);
        }
    } else {
        arma::uvec idxtau = arma::find(fixtau == 0);
        if (arma::sum(tauInit(idxtau)) == 0) {
            tau[0] = 1;
            tau(idxtau).fill(arma::var(Y) / tau.n_elem);
            if (std::abs(arma::var(Y)) < 0.1) {
                throw std::runtime_error("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T");
            }
        } else {
            tau(arma::find(fixtau == 0)) = tauInit(arma::find(fixtau == 0));
        }
    }
    arma::mat cov;
    const GLMResults* fit0ptr = &fit0; 
    Get_Coef(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);

    arma::mat AI;
    arma::vec YPAPY, Trace, PY;

    getAIScore(Y, X, W, tau, fixtau, Sigma_iY, Sigma_iX, cov, AI, YPAPY, Trace, PY, k1, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO, Ivec_start_indices, eMat);

    arma::vec tau0 = tau;
    arma::vec tau0_q2 = tau(arma::find(fixtau == 0));
    arma::vec tau_q2 = arma::max(arma::zeros(tau0_q2.n_elem), tau0_q2 + tau0_q2 % tau0_q2 % (re_YPAPY - re_Trace) / n);
    tau(arma::find(fixtau == 0)) = tau_q2;


    if (verbose) {
        std::cout << "Variance component estimates:" << std::endl;
        tau.print();
    }

    for (int i = 0; i < maxiter; ++i) {
        if (verbose) {
            std::cout << "\nIteration " << i << " " << tau << ":\n";
        }
        alpha0 = alpha;
        tau0 = tau;
        eta0 = eta;

        Get_Coef(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);

        fitglmmaiRPCG_multiV_updateTau(Y, X, W, tau, fixtau, Sigma_iY, Sigma_iX, cov, alpha, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff, LOCO);

        mu_eta = fit0.mu_eta(eta);
        sqrtW = mu_eta / arma::sqrt(1 / var_weights % fit0.variance(mu));
        W = arma::square(sqrtW);

        if (arma::max(arma::abs(tau - tau0) / (arma::abs(tau) + arma::abs(tau0) + tol)) < tol) {
            break;
        }

        if (arma::max(tau) > std::pow(tol, -2)) {
            std::cerr << "Large variance estimate observed in the iterations, model not converged..." << std::endl;
            i = maxiter;
            break;
        }
    }

    std::cout << "\nFinal " << tau << ":\n";
    Get_Coef(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);
    bool converged = (i <= maxiter);

    arma::vec res = y - mu;
    std::string traitType;

    arma::vec mu2 = (1/tau[0])*fit0.variance(mu);

    arma::vec mu2_rescaled = mu2 % var_weights;
    arma::vec y_rescaled = y % var_weights;
    arma::vec mu_rescaled = mu % var_weights;

    // Placeholder for ScoreTest_NULL_Model function
    SA_NULL scoreTestResult;
    if (!isCovariateOffset) {
        scoreTestResult = ScoreTest_NULL_Model(mu, mu2, y, X);
    } else {
        scoreTestResult = ScoreTest_NULL_Model(mu, mu2, y, Xorig);
    }
    
    glmmResult.theta = tau;
    glmmResult.coefficients = alpha;
    glmmResult.linear_predictors = eta;
    glmmResult.fitted_values = mu;
    glmmResult.Y = Y;
    glmmResult.residuals = res;
    glmmResult.cov = cov;
    glmmResult.converged = converged;
    glmmResult.sampleID = sampleIDsinPheno;
    glmmResult.obj_noK = fit0;
    glmmResult.y = y;
    glmmResult.X = X;
    glmmResult.traitType = traitType;
    glmmResult.isCovariateOffset = isCovariateOffset;
    glmmResult.varWeights = var_weights;
    glmmResult.obj_noK = scoreTestResult;
    glmmResult.LOCO = LOCO;



    if (!isLowMemLOCO && LOCO) {
        // Placeholder for set_Diagof_StdGeno_LOCO function
        ptr_gNULLGENOobj->set_Diagof_StdGeno_LOCO();
        arma::mat LOCOResult(22, 1);
        GLMMResults glmmResultLOCO;

        for (int j = 0; j < 22; ++j) {
            int startIndex = chromosomeStartIndexVec[j];
            int endIndex = chromosomeEndIndexVec[j];
            if (startIndex != -1 && endIndex != -1) {
                setStartEndIndex(startIndex, endIndex, j - 1)
                Get_Coef(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);
                glmmResultLOCO.coefficients = alpha;
                glmmResultLOCO.linear_predictors = eta;
                glmmResultLOCO.fitted_values = mu;
                glmmResultLOCO.Y = Y;
                glmmResultLOCO.residuals = res;
                glmmResultLOCO.cov = cov;
                arma::vec mu2LOCO = (1/tau[0])*fit0.variance(mu);
                SA_NULL scoreTestResultLOCO;
                if (!isCovariateOffset) {
                    scoreTestResultLOCO = ScoreTest_NULL_Model(mu, mu2LOCO, y, X);
                } else {
                    scoreTestResultLOCO = ScoreTest_NULL_Model(mu, mu2LOCO, y, Xorig);
                }
                glmmResultLOCO.obj_noK = scoreTestResultLOCO;
                glmmResultLOCO.LOCO = true;
                glmmResult.LOCOResult.push_back(glmmResultLOCO);

            } else {
                glmmResultLOCO.LOCO = false;   
            }
            glmmResult.LOCOResult.push_back(glmmResultLOCO);

        }
        glmmResult.chromosomeStartIndexVec = chromosomeStartIndexVec;
        glmmResult.chromosomeEndIndexVec =  chromosomeEndIndexVec;
        
    }
    return glmmResult;
}  


