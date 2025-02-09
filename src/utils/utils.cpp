#include <iostream>
#include <cstring>
#include <getopt.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <armadillo>

// Helper function to split a string by a delimiter and return a vector of strings
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Helper function to split a string by a delimiter and return a vector of doubles
std::vector<double> splitToDouble(const std::string& str, char delimiter) {
    std::vector<double> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(std::stod(token));
    }
    return tokens;
}

// Helper function to split a string by a delimiter and return a vector of ints
std::vector<int> splitToInt(const std::string& str, char delimiter) {
    std::vector<int> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(std::stoi(token));
    }
    return tokens;
}

// Helper function to convert a vector of strings to a vector of integers representing factor levels
std::vector<int> asFactor(const std::vector<std::string>& vec) {
    std::map<std::string, int> levels;
    std::vector<int> factorVec(vec.size());
    int level = 0;

    for (size_t i = 0; i < vec.size(); ++i) {
        if (levels.find(vec[i]) == levels.end()) {
            levels[vec[i]] = level++;
        }
        factorVec[i] = levels[vec[i]];
    }

    return factorVec;
}

// Function to find indices of elements that are in the first vector and also in the second vector
arma::uvec findSampleIndices(const std::vector<std::string>& firstVector, const std::vector<std::string>& secondVector) {
    std::unordered_map<std::string, arma::uword> elementIndexMap;
    for (arma::uword i = 0; i < secondVector.size(); ++i) {
        elementIndexMap[secondVector[i]] = i;
    }

    std::vector<arma::uword> indices;
    for (arma::uword i = 0; i < firstVector.size(); ++i) {
        if (elementIndexMap.find(firstVector[i]) != elementIndexMap.end()) {
            indices.push_back(i);
        }
    }

    return arma::uvec(indices);
}

// Function to subset a vector based on indices
template <typename T>
std::vector<T> subsetVector(const std::vector<T>& vec, const arma::uvec& indices) {
    std::vector<T> subset;
    for (arma::uword i : indices) {
        subset.push_back(vec[i]);
    }
    return subset;
}

struct CovariateTransformResult {
    arma::vec Y;
    arma::mat X1_new;
    arma::mat qrr;
    std::vector<std::string> X_name;
    arma::uvec idx_na;
};

CovariateTransformResult Covariate_Transform(const arma::mat& X1, const arma::vec& Y, const std::vector<std::string>& X_name) {
    arma::mat X1_copy = X1;
    arma::vec coef = arma::solve(X1_copy, Y);
    arma::uvec idx_na = arma::find_nonfinite(coef);

    if (!idx_na.is_empty()) {
        X1_copy.shed_cols(idx_na);
        std::cerr << "Warning: multi collinearity is detected in covariates! Columns will be excluded in the model\n";
    }

    arma::mat qrr;
    arma::mat X1_Q;
    arma::qr(X1_Q, qrr, X1_copy);

    int N = X1_copy.n_rows;
    arma::mat X1_new = X1_Q * std::sqrt(N);

    CovariateTransformResult result;
    result.Y = Y;
    result.X1_new = X1_new;
    result.qrr = qrr;
    result.X_name = {}; // Placeholder for column names
    result.idx_na = idx_na;

    return result;


void set_Ivec_start_indices(const std::vector<std::string>& modelID, arma::uvec & Ivec_start_indices) {
    std::vector<int> b = asFactor(modelID);
    arma::uvec unique_vals = arma::unique(b); // Automatically sorted
    Ivec_start_indices.set_size(unique_vals.n_elem); // Storage for start indices

    // Find start indices for each unique value
    for (size_t i = 0; i < unique_vals.n_elem; ++i) {
            Ivec_start_indices(i) = arma::as_scalar(arma::find(modelID == unique_vals(i), 1, "first"));
    }
    Ivec_start_indices.insert_rows(Ivec_start_indices.n_elem, 1);
    Ivec_start_indices[Ivec_start_indices.n_elem-1] = Ivec_start_indices.n_rows;
    /*int n = b.size();
    int max_b = *std::max_element(b.begin(), b.end()) + 1;
    
    arma::sp_mat I_mat(n, max_b);
    for (int i = 0; i < n; ++i) {
        I_mat(i, b[i]) = 1.0;
    }

    set_I_longl_mat(I_mat, b);
    */
}

// Function to check for perfect separation in a dataset
void checkPerfectSep(arma::mat& X1, const arma::vec& Y, std::vector<std::string>& X_name, int minCovariateCount) {
    std::vector<std::string> colnamesDelete;
    int q = X1.n_cols;

    for (int i = 0; i < q; ++i) {
        arma::vec column = X1.col(i);
        if (arma::unique(column).n_elem == 2) {
            arma::umat sumTable = arma::hist(Y, column);
            if (arma::any(sumTable < minCovariateCount)) {
                colnamesDelete.push_back(X_name[i]);
                std::cout << "less than " << minCovariateCount << " samples in a covariate detected! " << X_name[i] << " will be excluded in the model\n";
            }
            X1.shed_col(i);
            X_name.erase(X_name.begin() + i);
            --i; // Adjust index after erasing element
            --q; // Adjust total column count
        }
    }
    return colnamesDelete;
}