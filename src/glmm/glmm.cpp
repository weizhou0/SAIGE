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
    readPhenoFile(params.phenoFile, params.phenoCol, params.covarColList, params.categoricCovarCol, params.envCovarCol, params.MaleOnly, params.FemaleOnly, params.sexCol, params.FemaleCode, params.MaleCode, sampleIDsinPheno, X, eMat, pheno);
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

    std::string sampleIDsinPhenonew = sampleIDsinPheno

    // Subset the sparse GRM file
    if (!params.sparseGRMFile.empty()) {
        arma::uvec sampleIndicesSparseGRM = findSampleIndices(sampleIDsinsparseGRMFile, overlappedSamples);
        arma::sp_mat sparseGRM = readSparseGRM(params.sparseGRMFile, sampleIndicesSparseGRM);
    }





}



void fitNullGLMM(const Parameters& params) {
    std::cout << "Fitting null GLMM" << std::endl;
    // 
}
struct GLMResults {
    arma::vec beta;
    arma::vec eta;
    arma::vec mu;
    arma::vec W;
    arma::vec z;
};

GLMResults fitGLM(const arma::mat& X, const arma::vec& y, const std::string& family) {
    GLMResults results;
    results.beta = stats::glm(X, y, 
        (family == "logistic") ? stats::glm_family::logistic :
        (family == "poisson") ? stats::glm_family::poisson :
        (family == "gaussian") ? stats::glm_family::gaussian :
        throw std::invalid_argument("Invalid family type. Use 'logistic', 'poisson', or 'gaussian'.")
    );

    results.eta = X * results.beta;  // Linear predictor
    
    if (family == "logistic") {
        results.mu = 1 / (1 + arma::exp(-results.eta)); // Sigmoid function
        results.W = results.mu % (1 - results.mu);      // IRLS weights
    } else if (family == "poisson") {
        results.mu = arma::exp(results.eta); // Poisson mean
        results.W = results.mu;              // IRLS weights
    } else if (family == "gaussian") { 
        results.mu = results.eta; 
        results.W = arma::ones(results.mu.n_elem); // Constant weights for OLS
    }

    results.z = results.eta + (y - results.mu) / results.W; // Adjusted response

    return results;
}


//GLMResults logistic_results = fitGLM(X, y_logistic, "logistic");




// Function to fit the null generalized linear mixed model
arma::mat model_matrix(const arma::mat& fit0) {
  // Extract the model matrix from the fit0 object
  return fit0;
}

// Function to fit the null generalized linear mixed model
arma::mat glmmkin_ai_PCG_Rcpp_multiV_eMat(
    const std::string& bedFile,
    const std::string& bimFile,
    const std::string& famFile,
    const arma::mat& Xorig,
    bool isCovariateOffset,
    const arma::mat& fit0,
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
    const arma::mat& obj_noK,
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
    // Initialize variables
    arma::vec subSampleInGeno = subPheno.col(0);
    if (subSampleInGeno.is_empty()) {
        subSampleInGeno = subPheno.col(1);
    }

    if (verbose) {
        std::cout << "Start reading genotype plink file here" << std::endl;
    }

    // Set duplicate sample index
    arma::vec dup_sample_index = arma::conv_to<arma::vec>::from(arma::unique(subPheno.col(2)));

    if (verbose) {
        std::cout << "Genotype reading is done" << std::endl;
    }

    if (LOCO) {
        // Update chromosome start and end indices
        // ... (code to update chromosome indices)
    }

    arma::vec y = fit0.col(0);
    int n = y.n_elem;
    arma::mat X = model_matrix(fit0);
    arma::vec offset = fit0.col(1);
    if (offset.is_empty()) {
        offset = arma::zeros(n);
    }

    arma::vec var_weights = fit0.col(2);
    arma::vec eta = fit0.col(3);
    arma::vec mu = fit0.col(4);
    arma::vec mu_eta = eta; // Placeholder for family["mu.eta"](eta)
    arma::vec Y = eta - offset + (y - mu) / mu_eta;

    if (var_weights.is_empty()) {
        var_weights = arma::ones(mu_eta.n_elem);
    }

    arma::vec sqrtW = mu_eta / arma::sqrt(1 / var_weights % mu); // Placeholder for family["variance"](mu)
    arma::vec W = sqrtW % sqrtW;

    arma::vec alpha0 = fit0.col(5);
    arma::vec eta0 = eta;

    tau.zeros();
    if (true) { // Placeholder for family["family"] == "poisson" || family["family"] == "binomial"
        tau[0] = 1;
        fixtau[0] = 1;
        tauInit[0] = 1;
        arma::uvec idxtau = arma::find(fixtau == 0);
        if (arma::sum(tauInit(idxtau)) == 0) {
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

    if (!covarianceIdxMat.is_empty()) {
        arma::uvec idxtau2 = arma::intersect(covarianceIdxMat.col(0), arma::find(fixtau == 0));
        if (!idxtau2.is_empty()) {
            tau(idxtau2).zeros();
        }
        arma::vec Kmatdiag = arma::ones(tau.n_elem - 1); // Placeholder for getMeanDiagofKmat(LOCO)
        tau.subvec(1, tau.n_elem - 1) /= Kmatdiag;
    }

    if (isSparseGRMIdentity) {
        tau[1] = 0;
    }

    // Placeholder for Get_Coef_multiV function
    arma::mat re_coef_Y = Y;
    arma::vec re_coef_W = W;
    arma::vec re_coef_Sigma_iY = Y;
    arma::mat re_coef_Sigma_iX = X;
    arma::mat re_coef_cov = X;
    arma::vec re_coef_alpha = alpha0;
    arma::vec re_coef_eta = eta;
    arma::vec re_coef_mu = mu;

    // Placeholder for getAIScore_multiV_eMat function
    arma::vec re_YPAPY = arma::ones(n);
    arma::vec re_Trace = arma::ones(n);

    arma::vec tau0 = tau;
    arma::vec tau0_q2 = tau(arma::find(fixtau == 0));
    arma::vec tau_q2 = arma::max(arma::zeros(tau0_q2.n_elem), tau0_q2 + tau0_q2 % tau0_q2 % (re_YPAPY - re_Trace) / n);
    tau(arma::find(fixtau == 0)) = tau_q2;

    if (!covarianceIdxMat.is_empty()) {
        tau(arma::find(fixtau == 0) % arma::intersect(covarianceIdxMat.col(0), arma::find(fixtau == 0))).zeros();
    }

    if (verbose) {
        std::cout << "Variance component estimates:" << std::endl;
        tau.print();
    }

    int maxiter_in = maxiter;
    if (isSparseGRMIdentity) {
        tau[1] = 0;
        maxiter_in = 0;
        alpha0 = re_coef_alpha;
        tau0 = tau;
        eta0 = eta;
    }

    for (int i = 0; i < maxiter_in; ++i) {
        if (verbose) {
            std::cout << "\nIteration " << i << " " << tau << ":\n";
        }
        alpha0 = re_coef_alpha;
        tau0 = tau;
        eta0 = eta;

        // Placeholder for Get_Coef_multiV function
        re_coef_Y = Y;
        re_coef_W = W;
        re_coef_Sigma_iY = Y;
        re_coef_Sigma_iX = X;
        re_coef_cov = X;
        re_coef_alpha = alpha0;
        re_coef_eta = eta;
        re_coef_mu = mu;

        // Placeholder for fitglmmaiRPCG_multiV_eMat function
        arma::vec fit_tau = tau;
        arma::mat fit_cov = X;
        arma::vec fit_alpha = alpha0;
        arma::vec fit_eta = eta;
        arma::vec fit_Y = Y;
        arma::vec fit_mu = mu;

        tau = fit_tau;
        re_coef_cov = fit_cov;
        re_coef_alpha = fit_alpha;
        re_coef_eta = fit_eta;
        re_coef_Y = fit_Y;
        re_coef_mu = fit_mu;

        mu_eta = eta; // Placeholder for family["mu.eta"](eta)
        sqrtW = mu_eta / arma::sqrt(1 / var_weights % mu); // Placeholder for family["variance"](mu)
        W = sqrtW % sqrtW;

        if (arma::max(arma::abs(tau - tau0) / (arma::abs(tau) + arma::abs(tau0) + tol)) < tol) {
            break;
        }

        if (arma::max(tau) > std::pow(tol, -2)) {
            std::cerr << "Large variance estimate observed in the iterations, model not converged..." << std::endl;
            i = maxiter;
            break;
        }
    }

    if (verbose) {
        std::cout << "\nFinal " << tau << ":\n";
    }

    // Placeholder for Get_Coef_multiV function
    re_coef_Y = Y;
    re_coef_W = W;
    re_coef_Sigma_iY = Y;
    re_coef_Sigma_iX = X;
    re_coef_cov = X;
    re_coef_alpha = alpha0;
    re_coef_eta = eta;
    re_coef_mu = mu;

    arma::vec res = y - re_coef_mu;
    arma::vec mu2;
    std::string traitType;
    if (true) { // Placeholder for family["family"] == "binomial"
        mu2 = re_coef_mu % (1 - re_coef_mu);
        traitType = "binary";
    } else if (false) { // Placeholder for family["family"] == "poisson"
        mu2 = re_coef_mu;
        traitType = "count";
    } else if (false) { // Placeholder for family["family"] == "gaussian"
        mu2 = arma::ones(res.n_elem) / tau[0];
        traitType = "quantitative";
    }

    arma::vec mu2_rescaled = mu2 % var_weights;
    arma::vec y_rescaled = y % var_weights;
    arma::vec mu_rescaled = re_coef_mu % var_weights;

    // Placeholder for ScoreTest_NULL_Model function
    arma::mat obj_noK_rescaled = X;

    arma::mat glmmResult = arma::join_horiz(tau, re_coef_alpha, re_coef_eta, re_coef_mu, re_coef_Y, res, re_coef_cov, arma::conv_to<arma::vec>::from({(int)(maxiter_in < maxiter)}), subPheno.col(2), obj_noK_rescaled, y, X, arma::conv_to<arma::vec>::from({(int)isCovariateOffset}), var_weights, arma::conv_to<arma::vec>::from({(int)LOCO}));

    if (!isLowMemLOCO && LOCO) {
        // Placeholder for set_Diagof_StdGeno_LOCO function
        arma::mat LOCOResult(22, 1);
        for (int j = 0; j < 22; ++j) {
            int startIndex = chromosomeStartIndexVec[j];
            int endIndex = chromosomeEndIndexVec[j];
            if (startIndex != -1 && endIndex != -1) {
                // Placeholder for setStartEndIndex function
                // Placeholder for Get_Coef_multiV function
                re_coef_Y = Y;
                re_coef_W = W;
                re_coef_Sigma_iY = Y;
                re_coef_Sigma_iX = X;
                re_coef_cov = X;
                re_coef_alpha = alpha0;
                re_coef_eta = eta;
                re_coef_mu = mu;

                // Placeholder for ScoreTest_NULL_Model function
                obj_noK_rescaled = X;

                LOCOResult.row(j) = arma::join_horiz(arma::conv_to<arma::vec>::from({1}), re_coef_alpha, re_coef_eta, re_coef_mu, re_coef_Y, y - re_coef_mu, re_coef_cov, obj_noK_rescaled);
            } else {
                LOCOResult.row(j) = arma::conv_to<arma::vec>::from({0});
            }
        }
        glmmResult = arma::join_horiz(glmmResult, LOCOResult);
    }

    if (isLowMemLOCO && LOCO) {
        glmmResult = arma::join_horiz(glmmResult, chromosomeStartIndexVec, chromosomeEndIndexVec);
    }

    return glmmResult;
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

    // Add factor covariates to combinedCovariates if not empty
    if (!factorCovariates.empty()) {
        for (size_t i = 0; i < factorCovariates.size(); ++i) {
            for (size_t j = 0; j < factorCovariates[i].size(); ++j) {
                combinedCovariates[j].push_back(static_cast<double>(factorCovariates[i][j]));
            }
        }
    }

    if (!covariateData_numeric.empty()) {
        combinedCovariates.resize(covariateData_numeric.size());
        for (size_t i = 0; i < covariateData_numeric.size(); ++i) {
            combinedCovariates[i] = covariateData_numeric[i];
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

