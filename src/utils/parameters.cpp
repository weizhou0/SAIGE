#include "parameters.h"
#include <iostream>
#include <cstring>

bool Parameters::parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--step") == 0 && i + 1 < argc) {
            step = argv[++i];
        } else if (strcmp(argv[i], "--input") == 0 && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (strcmp(argv[i], "--test") == 0 && i + 1 < argc) {
            testType = argv[++i];
        } else {
            std::cerr << "Unknown parameter: " << argv[i] << std::endl;
            return false;
        }
    }
    return true;
}
