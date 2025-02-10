#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include "utils.h"

// Structure to hold the parameters for Step 1
struct Step1Parameters {
    std::string plinkFile;
    std::string bedFile;
    std::string bimFile;
    std::string famFile;
    std::string phenoFile;
    std::string phenoCol;
    bool isRemoveZerosinPheno;
    std::string traitType;
    bool invNormalize;
    std::vector<std::string> covarColList;
    std::vector<std::string> qCovarCol;
    std::vector<std::string> eCovarCol;
    std::vector<std::string> sampleCovarCol;
    std::vector<std::string> offsetCol;
    std::vector<std::string> varWeightsCol;
    std::string longlCol;
    std::string sampleIDColinphenoFile;
    double tol;
    int maxiter;
    double tolPCG;
    int maxiterPCG;
    int nThreads;
    double SPAcutoff;
    int numMarkersForVarRatio;
    bool skipModelFitting;
    double memoryChunk;
    std::vector<double> tauInit;
    bool LOCO;
    bool isLowMemLOCO;
    double traceCVcutoff;
    double ratioCVcutoff;
    std::string outputPrefix;
    std::string outputPrefix_varRatio;
    bool IsOverwriteVarianceRatioFile;
    std::string sparseGRMFile;
    std::string sparseGRMSampleIDFile;
    int numRandomMarkerforSparseKin;
    double relatednessCutoff;
    bool isCateVarianceRatio;
    std::vector<int> cateVarRatioIndexVec;
    std::vector<double> cateVarRatioMinMACVecExclude;
    std::vector<double> cateVarRatioMaxMACVecInclude;
    bool isCovariateTransform;
    bool isDiagofKinSetAsOne;
    int minCovariateCount;
    double minMAFforGRM;
    double maxMissingRateforGRM;
    bool useSparseGRMtoFitNULL;
    bool useSparseGRMforVarRatio;
    bool includeNonautoMarkersforVarRatio;
    std::string sexCol;
    int FemaleCode;
    bool FemaleOnly;
    int MaleCode;
    bool MaleOnly;
    std::string SampleIDIncludeFile;
    bool isCovariateOffset;
    bool skipVarianceRatioEstimation;
    int nrun;
    std::string VmatFilelist;
    std::string VmatSampleFilelist;
    std::string VcellmatFilelist;
    std::string VcellmatSampleFilelist;
    bool useGRMtoFitNULL;
    bool isStoreSigma;
    bool isShrinkModelOutput;
};

// Structure to hold the parameters for Step 2
struct Step2Parameters {
    std::string vcfFile;
    std::string vcfFileIndex;
    std::string vcfField;
    std::string savFile;
    std::string savFileIndex;
    std::string bgenFile;
    std::string bgenFileIndex;
    std::string sampleFile;
    std::string bedFile;
    std::string bimFile;
    std::string famFile;
    std::string AlleleOrder;
    std::string idstoIncludeFile;
    std::string rangestoIncludeFile;
    std::string chrom;
    bool is_imputed_data;
    double minMAF;
    double minMAC;
    double minGroupMAC_in_BurdenTest;
    double minInfo;
    double maxMissing;
    std::string impute_method;
    bool LOCO;
    std::string GMMATmodelFile;
    std::string varianceRatioFile;
    std::string SAIGEOutputFile;
    int markers_per_chunk;
    int groups_per_chunk;
    bool is_output_moreDetails;
    bool is_overwrite_output;
    std::string maxMAF_in_groupTest;
    std::string maxMAC_in_groupTest;
    std::string annotation_in_groupTest;
    std::string groupFile;
    std::string sparseGRMFile;
    std::string sparseGRMSampleIDFile;
    double relatednessCutoff;
    double MACCutoff_to_CollapseUltraRare;
    std::string cateVarRatioMinMACVecExclude;
    std::string cateVarRatioMaxMACVecInclude;
    std::string weights_beta;
    double r_corr;
    int markers_per_chunk_in_groupTest;
    std::string condition;
    std::string weights_for_condition;
    double SPAcutoff;
    double dosage_zerod_cutoff;
    double dosage_zerod_MAC_cutoff;
    bool is_single_in_groupTest;
    bool is_no_weight_in_groupTest;
    bool is_output_markerList_in_groupTest;
    bool is_Firth_beta;
    double pCutoffforFirth;
    bool is_fastTest;
    double max_MAC_for_ER;
};

// Structure to hold the overall parameters
struct Parameters {
    std::string step;
    Step1Parameters step1Params;
    Step2Parameters step2Params;

    bool parse(int argc, char* argv[]);
};

// Function to parse Step 1 parameters
bool parseStep1Parameters(int argc, char* argv[], Step1Parameters& params);

// Function to parse Step 2 parameters
bool parseStep2Parameters(int argc, char* argv[], Step2Parameters& params);

// Function to check if a file exists
bool fileExists(const std::string& filename);

// Function to print the parameters
void printParameters(const Parameters& params);

// Function to parse command-line arguments
void parseArguments(int argc, char* argv[], Parameters& params);

// Function to check if files exist
void checkFiles(const std::vector<std::string>& filePaths);

// Function to preprocess data
void preprocessData(const Parameters& params);

#endif // PARAMETERS_H