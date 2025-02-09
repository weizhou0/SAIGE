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