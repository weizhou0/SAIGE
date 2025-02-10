#ifndef VALIDATION_H
#define VALIDATION_H

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

// Function to check if a string can be converted to a boolean
bool isBool(const std::string& str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return (lowerStr == "true" || lowerStr == "false");
}

// Function to check if a string can be converted to an integer
bool isInt(const std::string& str) {
    std::istringstream iss(str);
    int val;
    return (iss >> val >> std::ws).eof();
}

// Function to check if a string can be converted to a double
bool isDouble(const std::string& str) {
    std::istringstream iss(str);
    double val;
    return (iss >> val >> std::ws).eof();
}

// Function to check if a string can be split into a vector of doubles
bool isDoubleVector(const std::string& str, char delimiter) {
    std::istringstream tokenStream(str);
    std::string token;
    while (std::getline(tokenStream, token, delimiter)) {
        if (!isDouble(token)) {
            return false;
        }
    }
    return true;
}

// Function to check if a string can be split into a vector of integers
bool isIntVector(const std::string& str, char delimiter) {
    std::istringstream tokenStream(str);
    std::string token;
    while (std::getline(tokenStream, token, delimiter)) {
        if (!isInt(token)) {
            return false;
        }
    }
    return true;
}

#endif // VALIDATION_H