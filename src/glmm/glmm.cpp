#include "glmm.h"
#include <iostream>
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include "stats.hpp"
#include <unordered_map>



void checkFiles(const std::vector<std::string>& filePaths) {
    for (const auto& filePath : filePaths) {
        if (!filePath.empty()) {
            if (!std::filesystem::exists(filePath)) {
                throw std::runtime_error("ERROR! File does not exist: " + filePath);
            }
        }
    }
}
void preprocessData(const Parameters& params) {
    
    std::cout << "Preprocessing data" << std::endl;
    // Check if the bed, bim, and fam files exist
    if (!std::filesystem::exists(params.bedFile)) {
        throw std::runtime_error("ERROR! bed file does not exist");
    }
    if (!std::filesystem::exists(params.bimFile)) {
        throw std::runtime_error("ERROR! bim file does not exist");
    }
    if (!std::filesystem::exists(params.famFile)) {
        throw std::runtime_error("ERROR! fam file does not exist");
    }


    if (!fileExists(famFile)) {
        stop("ERROR! fam file does not exist\n");
      } else {
        if (famFile.empty()) {
        stop("ERROR! fam file is not specified\n");
        }
        std::vector<std::string> sampleListwithGenov0 = readFile(famFile);
        std::vector<std::string> IIDgeno;
        for (const auto& line : sampleListwithGenov0) {
        std::vector<std::string> tokens = split(line, ' ');
        IIDgeno.push_back(tokens[1]);
        }
        std::vector<int> IndexGeno(IIDgeno.size());
        std::iota(IndexGeno.begin(), IndexGeno.end(), 1);
        std::cout << IIDgeno.size() << " samples have genotypes\n";
      }
      } else {
      if (useSparseGRMtoFitNULL || useSparseGRMforVarRatio) {
        std::vector<std::string> sampleListwithGenov0 = readFile(sparseGRMSampleIDFile);
        std::vector<std::string> IIDgeno;
        for (const auto& line : sampleListwithGenov0) {
        IIDgeno.push_back(line);
        }
        std::vector<int> IndexGeno(IIDgeno.size());
        std::iota(IndexGeno.begin(), IndexGeno.end(), 1);
        std::cout << IIDgeno.size() << " samples are in the sparse GRM\n";
      }
      }



    // Read only the second column (sample IDs) from the fam file and store in an arma::mat
    std::ifstream famFile(params.famFile);
    if (!famFile.is_open()) {
        throw std::runtime_error("ERROR! Unable to open fam file");
    }

    std::vector<std::string> sampleIDs;
    std::string line;
    while (std::getline(famFile, line)) {
        std::istringstream iss(line);
        std::string id;
        for (int i = 0; i < 2; ++i) {
            iss >> id;
        }
        sampleIDs.push_back(id);
    }
    famFile.close();

    // Convert sampleIDs to arma::mat
    arma::mat famData(sampleIDs.size(), 1);
    for (size_t i = 0; i < sampleIDs.size(); ++i) {
        famData(i, 0) = std::stod(sampleIDs[i]);
    }

    arma::mat famData;
    famData.load(params.famFile, arma::csv_ascii);

    // Extract sample IDs
    std::vector<std::string> sampleIDs;
    for (size_t i = 0; i < famData.n_rows; ++i) {
        sampleIDs.push_back(famData(i, 1));
    }

    // Print the number of samples with genotypes
    std::cout << sampleIDs.size() << " samples have genotypes" << std::endl;





}


void fitNULLGLMM_multiV(const Parameters& params) {
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
    std::vector<std::string> sampleIDsinPheno, sampleIDsinfamFile, sampleIDsinsparseGRMFile, sampleIDsinSampleIncludeFile;
    arma::mat X, eMat;
    arma::vec pheno;
    std::string Xname;
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
    if (hasDuplicates) {
        std::cout << "Duplicated IDs were found " << std::endl;
        std::cout << overlappedSamples.size() << " observations will be used for analysis" << std::endl;
        set_Ivec_start_indices(sampleIDsinPhenonew, Ivec_start_indices);

    }

    // Output the results
    //std::cout << "Transformed Y: " << result.Y.t() << std::endl;
    //std::cout << "Transformed X1_new: " << result.X1_new << std::endl;
    //std::cout << "QR decomposition matrix qrr: " << result.qrr << std::endl;
    //std::cout << "Indices of non-finite coefficients: " << result.idx_na.t() << std::endl;
 

}


// Function to fit the null generalized linear mixed model
GLMMResults glmmkin_ai_PCG_Rcpp_multiV(
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
    bool isSparseGRMIdentity,
    NullGENO::NullGenoClass* ptr_gNULLGENOobj
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
    Get_Coef_multiV(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);
    bool converged = (i <= maxiter);

    arma::vec res = y - mu;
    arma::vec mu2;
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
    glmmResult.LOCO = LOCO;



    if (!isLowMemLOCO && LOCO) {
        // Placeholder for set_Diagof_StdGeno_LOCO function
        ptr_gNULLGENOobj->set_Diagof_StdGeno_LOCO();
        arma::mat LOCOResult(22, 1);
        for (int j = 0; j < 22; ++j) {
            int startIndex = chromosomeStartIndexVec[j];
            int endIndex = chromosomeEndIndexVec[j];
            if (startIndex != -1 && endIndex != -1) {
                setStartEndIndex(startIndex, endIndex, j - 1)
                Get_Coef_multiV(y, X, tau, offset, var_weights, fit0ptr, alpha, eta, Sigma_iY, Sigma_iX, cov, Y, maxiterPCG, tolPCG, maxiter, verbose, LOCO);

        
                re_coef_W = W;
                re_coef_Sigma_iY = Y; = X;
                re_coef_Sigma_iX = X;
                re_coef_cov = X;pha0;
                re_coef_alpha = alpha0;a;
                re_coef_eta = eta;                
                re_coef_mu = mu;
                re_coef_mu = mu;
oreTest_NULL_Model function
                // Placeholder for ScoreTest_NULL_Model function                obj_noK_rescaled = X;
                obj_noK_rescaled = X;
Result.row(j) = arma::join_horiz(arma::conv_to<arma::vec>::from({1}), re_coef_alpha, re_coef_eta, re_coef_mu, re_coef_Y, y - re_coef_mu, re_coef_cov, obj_noK_rescaled);
                LOCOResult.row(j) = arma::join_horiz(arma::conv_to<arma::vec>::from({1}), re_coef_alpha, re_coef_eta, re_coef_mu, re_coef_Y, y - re_coef_mu, re_coef_cov, obj_noK_rescaled);
            } else {   LOCOResult.row(j) = arma::conv_to<arma::vec>::from({0});
                LOCOResult.row(j) = arma::conv_to<arma::vec>::from({0});   }
            }
        }   glmmResult = arma::join_horiz(glmmResult, LOCOResult);
        glmmResult.result = arma::join_horiz(glmmResult.result, LOCOResult);    }
    }

    if (isLowMemLOCO && LOCO) {   glmmResult = arma::join_horiz(glmmResult, chromosomeStartIndexVec, chromosomeEndIndexVec);
        glmmResult.result = arma::join_horiz(glmmResult.result, chromosomeStartIndexVec, chromosomeEndIndexVec);    }
    }

    return glmmResult;
}   return glmmResult;
}
    glmmResult["chromosomeEndIndexVec"] = chromosomeEndIndexVec;
  }
  
  return glmmResult;
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

