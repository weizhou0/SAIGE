#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

struct Parameters {
    std::string step;
    std::string inputFile;
    std::string outputFile;
    std::string testType;
    // Add other parameters as needed

    bool parse(int argc, char* argv[]);
};

#endif // PARAMETERS_H
