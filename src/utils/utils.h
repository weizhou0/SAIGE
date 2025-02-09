#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <sstream>

// Helper function to split a string by a delimiter and return a vector of strings
std::vector<std::string> split(const std::string& str, char delimiter);

// Helper function to split a string by a delimiter and return a vector of doubles
std::vector<double> splitToDouble(const std::string& str, char delimiter);

// Helper function to split a string by a delimiter and return a vector of ints
std::vector<int> splitToInt(const std::string& str, char delimiter);

#endif // UTILS_H