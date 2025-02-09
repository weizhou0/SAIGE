#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <fstream>


class GLMMResult {
    public:
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
        GLMMResult() {}
    
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

    //Serialize the object to a file:
    void saveGLMMResult(const GLMMResult& result, const std::string& filename) {
        std::ofstream ofs(filename);
        boost::archive::text_oarchive oa(ofs);
        oa << result;
    }
    //Deserialize the object from a file:
    GLMMResult loadGLMMResult(const std::string& filename) {
        GLMMResult result;
        std::ifstream ifs(filename);
        boost::archive::text_iarchive ia(ifs);
        ia >> result;
        return result;
    }

    int main() {
        // Create an example GLMMResult object
        GLMMResult result;
        result.theta = arma::vec({1.0, 2.0, 3.0});
        result.coefficients = arma::vec({0.1, 0.2, 0.3});
        result.linear_predictors = arma::vec({0.5, 0.6, 0.7});
        result.fitted_values = arma::vec({0.8, 0.9, 1.0});
        result.Y = arma::vec({1.1, 1.2, 1.3});
        result.residuals = arma::vec({1.4, 1.5, 1.6});
        result.cov = arma::mat(3, 3, arma::fill::eye);
        result.converged = true;
        result.sampleID = {"sample1", "sample2", "sample3"};
        result.y = arma::vec({2.1, 2.2, 2.3});
        result.X = arma::mat(3, 3, arma::fill::randu);
        result.traitType = "binary";
        result.isCovariateOffset = false;
        result.varWeights = arma::vec({1.0, 1.0, 1.0});
        result.LOCO = false;
        result.chromosomeStartIndexVec = arma::vec({0, 1, 2});
        result.chromosomeEndIndexVec = arma::vec({3, 4, 5});
    
        // Save the object to a file
        saveGLMMResult(result, "glmm_result.txt");
    
        // Load the object from the file
        GLMMResult loadedResult = loadGLMMResult("glmm_result.txt");
    
        // Print some values to verify
        std::cout << "Theta: " << loadedResult.theta.t() << std::endl;
        std::cout << "Coefficients: " << loadedResult.coefficients.t() << std::endl;
    
        return 0;
    }