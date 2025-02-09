#include "glmm.h"
#include <iostream>
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include "stats.hpp"
#include <unordered_map>

struct GLMResults {
    arma::vec beta;
    arma::vec eta;
    arma::vec mu;
    arma::vec W;
    arma::vec z;
    std::string familyGLM};
    ;


GLMResults fitGLM(const arma::mat& X, const arma::vec& y, const std::string& traitType) {
    GLMResults results;
    if (traitType == "binary") {
        results.beta = stats::glm(X, y, stats::glm_family::binomial);
        results.family = "binomial";
    } else if (traitType == "quantitative") {
        results.beta = stats::glm(X, y, stats::glm_family::gaussian);
        results.family = "gaussian";
    } else if (traitType == "tte") {
        results.beta = stats::glm(X, y, stats::glm_family::poisson);
        results.family = "poisson";
    } else if (traitType == "count") {
        results.beta = stats::glm(X, y, stats::glm_family::poisson);
        results.family = "poisson";
    } else {
        throw std::runtime_error("ERROR! Invalid trait type");
    }

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
    // Define the linkinv function for the respective family
    if (family == "logistic") {
        results.linkinv = [](const arma::vec& eta) { return 1 / (1 + arma::exp(-eta)); }; // Inverse of logit link
    } else if (family == "poisson") {
        results.linkinv = [](const arma::vec& eta) { return arma::exp(eta); }; // Inverse of log link
    } else if (family == "gaussian") {
        results.linkinv = [](const arma::vec& eta) { return eta; }; // Identity function for Gaussian
    }

    // Add the mu.eta function for the respective family
    if (family == "logistic") {
        results.mu_eta = [](const arma::vec& eta) { 
            arma::vec mu = linkinv(eta);  // First, compute mu
            return variance(mu); }; // Sigmoid function
        } else if (family == "poisson") {
        results.mu_eta = [](const arma::vec& eta) { return arma::exp(eta); }; // Poisson mean
        } else if (family == "gaussian") {
        results.mu_eta = [](const arma::vec& eta) { return arma::ones<arma::vec>(eta.n_elem); }; // Identity function for Gaussian
    }


    // Define family["variance"] as a function of mu
    if (family == "logistic") {
        results.variance = results.mu % (1 - results.mu); // Variance for logistic
    } else if (family == "poisson") {
        results.variance = results.mu; // Variance for Poisson
    } else if (family == "gaussian") {
        results.variance = arma::ones(results.mu.n_elem); // Constant variance for Gaussian
    }

    results.z = results.eta + (y - results.mu) / results.W; // Adjusted response
    results.familyGLM = family;
    return results;
}


//GLMResults logistic_results = fitGLM(X, y_logistic, "logistic");