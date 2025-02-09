#include "spgrm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <stdexcept>
#include <vector>
#include <string>

// Function to read a sparse GRM file
arma::sp_mat readSparseGRM(const std::string& sparseGRMFile, const arma::uvec& sampleIndices) {
    std::ifstream file(sparseGRMFile);
    if (!file.is_open()) {
        throw std::runtime_error("ERROR! Unable to open sparse GRM file");
    }

    std::vector<arma::uword> x, y;
    std::vector<double> z;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        arma::uword xi, yi;
        double zi;
        if (!(iss >> xi >> yi >> zi)) {
            throw std::runtime_error("ERROR! Invalid format in sparse GRM file");
        }
        x.push_back(xi);
        y.push_back(yi);
        z.push_back(zi);
    }
    file.close();

    // Create a sparse matrix
    arma::sp_mat sparseGRM(sampleIndices.n_elem, sampleIndices.n_elem);
    for (size_t i = 0; i < x.size(); ++i) {
        if (sampleIndices.is_elem(x[i]) && sampleIndices.is_elem(y[i])) {
            sparseGRM(sampleIndices(x[i]), sampleIndices(y[i])) = z[i];
        }
    }

    return sparseGRM;
}