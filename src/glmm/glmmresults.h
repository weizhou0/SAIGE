#ifndef GLMMRESULTS_H
#define GLMMRESULTS_H

#include <armadillo>
#include <string>
#include <vector>
#include <functional>

// Structure to hold the result of the GLMM fitting
struct GLMMResults {
    arma::vec theta;
    arma::vec coefficients;
    arma::vec linear_predictors;
    arma::vec fitted_values;
    arma::vec Y;
    arma::vec residuals;
    arma::mat cov;
    bool converged;
    std::vector<std::string> sampleID;
    // Assuming obj_noK is a serializable type
    // List obj_noK;
    arma::vec y;
    arma::mat X;
    std::string traitType;
    bool isCovariateOffset;
    arma::vec varWeights;
    bool LOCO;
    // Assuming LOCOResult is a serializable type
    // List LOCOResult;
    arma::vec chromosomeStartIndexVec;
    arma::vec chromosomeEndIndexVec;

    // Default constructor
    GLMMResults() {}

    // Serialization function
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & theta;
        ar & coefficients;
        ar & linear_predictors;
        ar & fitted_values;
        ar & Y;
        ar & residuals;
        ar & cov;
        ar & converged;
        ar & sampleID;
        // ar & obj_noK;
        ar & y;
        ar & X;
        ar & traitType;
        ar & isCovariateOffset;
        ar & varWeights;
        ar & LOCO;
        // ar & LOCOResult;
        ar & chromosomeStartIndexVec;
        ar & chromosomeEndIndexVec;
    }
};

#endif // GLMMRESULTS_H