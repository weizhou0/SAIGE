#ifndef GLM_H
#define GLM_H

#include <armadillo>
#include <string>
#include <vector>
#include <functional>

// Structure to hold the result of the GLM fitting
struct GLMResults {
    arma::vec beta;
    arma::vec eta;
    arma::vec mu;
    arma::vec W;
    arma::vec z;
    std::string familyGLM;
    std::function<arma::vec(const arma::vec&)> linkinv;
    std::function<arma::vec(const arma::vec&)> mu_eta;
    arma::vec variance;
};

// Function to fit a Generalized Linear Model (GLM)
GLMResults fitGLM(const arma::mat& X, const arma::vec& y, const std::string& traitType);

#endif // GLM_H