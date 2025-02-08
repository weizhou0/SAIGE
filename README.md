# SAIGE
C++ version of SAIGE

# SAIGE

## Overview
The `SAIGE` project is a C++ application designed to perform multiple tasks related to genetic data analysis, including creating a sparse GRM, fitting a null generalized linear mixed model (GLMM), conducting genetic association tests, and writing out an LD matrix. This program provides researchers with tools to analyze genetic data while accounting for random effects, making it suitable for various genetic studies.

## Features
- Create a sparse Genetic Relationship Matrix (GRM).
- Fit null generalized linear mixed models.
- Conduct genetic association tests (single-variant and set-based).
- Write out Linkage Disequilibrium (LD) matrix.
- Output model diagnostics and results.
- Easy integration with CMake for building and managing dependencies.

## Requirements
- C++11 or higher
- CMake 3.10 or higher

## Building the Project
To build the project, follow these steps:

1. Clone the repository:

git clone <repository-url> cd SAIGE


2. Create a build directory:

mkdir build
cd build

3. Run CMake to configure the project:
cmake ..

4. Compile the project:
make

## Running the Application
After building the project, you can run the application using the following command:

./SAIGE <step> [options]

## Usage
The application supports the following steps:

- **step0**: Create a sparse GRM
- **step1**: Fit a null GLMM
- **step2**: Conduct genetic association tests (single-variant or set-based)
- **step3**: Write out LD matrix

### Example Commands
- Create a sparse GRM:
  ```sh
  ./SAIGE step0 --input <input-file> --output <output-file> --param1 <value1> --param2 <value2> ...

GitHub Copilot
Sure, here's the updated README.md file with the correct project name and the steps renamed to step0, step1, step2, and step3:

git clone <repository-url> cd SAIGE

2. Create a build directory:
mkdir build cd build

3. Run CMake to configure the project:
cmake ..


4. Compile the project:
make

## Running the Application
After building the project, you can run the application using the following command:
./SAIGE <step> [options]

## Usage
The application supports the following steps:

- **step0**: Create a sparse GRM
- **step1**: Fit a null GLMM
- **step2**: Conduct genetic association tests (single-variant or set-based)
- **step3**: Write out LD matrix

### Example Commands
- Create a sparse GRM:
  ```sh
  ./SAIGE step0 --input <input-file> --output <output-file> --param1 <value1> --param2 <value2> ...
- Fit a null GLMM:
./SAIGE step1 --input <input-file> --output <output-file> --param1 <value1> --param2 <value2> ...

- Conduct genetic association tests:
./SAIGE step2 --input <input-file> --output <output-file> --test <single|set> --param1 <value1> --param2 <value2> ...

- Write out LD matrix:
./SAIGE step3 --input <input-file> --output <output-file> --param1 <value1> --param2 <value2> ...

Example
To fit a model and conduct an association test, provide the necessary input files and parameters when prompted. The results will be displayed in the console and can also be saved to a specified output file.

Contributing
Contributions are welcome! Please submit a pull request or open an issue for any enhancements or bug fixes.

License
This project is licensed under the MIT License. See the LICENSE file for more details. ```









