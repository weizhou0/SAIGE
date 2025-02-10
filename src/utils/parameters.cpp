#include "parameters.h"
#include "utils.h"
#include "validation.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <numeric>

bool Parameters::parse(int argc, char* argv[]) {
    if (strcmp(argv[1], "--step") == 0 && argc > 2) {
        step = argv[2];
    } else {
        std::cerr << "Step parameter is required" << std::endl;
        return false;
    }

    if (step == "1") {
        return parseStep1Parameters(argc, argv, step1Params);
    } else if (step == "2") {
        return parseStep2Parameters(argc, argv, step2Params);
    } else {
        std::cerr << "Unknown step: " << step << std::endl;
        return false;
    }
}

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
};;

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

struct Parameters {
    std::string step;
    Step1Parameters step1Params;
    Step2Parameters step2Params;

    bool parse(int argc, char* argv[]);
};


bool parseStep1Parameters(int argc, char* argv[], Step1Parameters& params) {
    static struct option long_options[] = {
        {"plinkFile", required_argument, 0, 0},
        {"bedFile", required_argument, 0, 0},
        {"bimFile", required_argument, 0, 0},
        {"famFile", required_argument, 0, 0},
        {"phenoFile", required_argument, 0, 0},
        {"phenoCol", required_argument, 0, 0},
        {"isRemoveZerosinPheno", no_argument, 0, 0},
        {"traitType", required_argument, 0, 0},
        {"invNormalize", no_argument, 0, 0},
        {"covarColList", required_argument, 0, 0},
        {"qCovarCol", required_argument, 0, 0},
        {"eCovarCol", required_argument, 0, 0},
        {"sampleCovarCol", required_argument, 0, 0},
        {"offsetCol", required_argument, 0, 0},
        {"varWeightsCol", required_argument, 0, 0},
        {"longlCol", required_argument, 0, 0},
        {"sampleIDColinphenoFile", required_argument, 0, 0},
        {"tol", required_argument, 0, 0},
        {"maxiter", required_argument, 0, 0},
        {"tolPCG", required_argument, 0, 0},
        {"maxiterPCG", required_argument, 0, 0},
        {"nThreads", required_argument, 0, 0},
        {"SPAcutoff", required_argument, 0, 0},
        {"numMarkersForVarRatio", required_argument, 0, 0},
        {"skipModelFitting", no_argument, 0, 0},
        {"memoryChunk", required_argument, 0, 0},
        {"tauInit", required_argument, 0, 0},
        {"LOCO", no_argument, 0, 0},
        {"isLowMemLOCO", no_argument, 0, 0},
        {"traceCVcutoff", required_argument, 0, 0},
        {"ratioCVcutoff", required_argument, 0, 0},
        {"outputPrefix", required_argument, 0, 0},
        {"outputPrefix_varRatio", required_argument, 0, 0},
        {"IsOverwriteVarianceRatioFile", no_argument, 0, 0},
        {"sparseGRMFile", required_argument, 0, 0},
        {"sparseGRMSampleIDFile", required_argument, 0, 0},
        {"numRandomMarkerforSparseKin", required_argument, 0, 0},
        {"relatednessCutoff", required_argument, 0, 0},
        {"isCateVarianceRatio", no_argument, 0, 0},
        {"cateVarRatioIndexVec", required_argument, 0, 0},
        {"cateVarRatioMinMACVecExclude", required_argument, 0, 0},
        {"cateVarRatioMaxMACVecInclude", required_argument, 0, 0},
        {"isCovariateTransform", no_argument, 0, 0},
        {"isDiagofKinSetAsOne", no_argument, 0, 0},
        {"minCovariateCount", required_argument, 0, 0},
        {"minMAFforGRM", required_argument, 0, 0},
        {"maxMissingRateforGRM", required_argument, 0, 0},
        {"useSparseGRMtoFitNULL", no_argument, 0, 0},
        {"useSparseGRMforVarRatio", no_argument, 0, 0},
        {"includeNonautoMarkersforVarRatio", no_argument, 0, 0},
        {"sexCol", required_argument, 0, 0},
        {"FemaleCode", required_argument, 0, 0},
        {"FemaleOnly", no_argument, 0, 0},
        {"MaleCode", required_argument, 0, 0},
        {"MaleOnly", no_argument, 0, 0},
        {"SampleIDIncludeFile", required_argument, 0, 0},
        {"isCovariateOffset", no_argument, 0, 0},
        {"skipVarianceRatioEstimation", no_argument, 0, 0},
        {"nrun", required_argument, 0, 0},
        {"VmatFilelist", required_argument, 0, 0},
        {"VmatSampleFilelist", required_argument, 0, 0},
        {"VcellmatFilelist", required_argument, 0, 0},
        {"VcellmatSampleFilelist", required_argument, 0, 0},
        {"useGRMtoFitNULL", no_argument, 0, 0},
        {"isStoreSigma", no_argument, 0, 0},
        {"isShrinkModelOutput", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    while (true) {
        int c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                switch (option_index) {
                    case 0: params.plinkFile = optarg; break;
                    case 1: params.bedFile = optarg; break;
                    case 2: params.bimFile = optarg; break;
                    case 3: params.famFile = optarg; break;
                    case 4: params.phenoFile = optarg; break;
                    case 5: params.phenoCol = optarg; break;
                    case 6: params.isRemoveZerosinPheno = true; break;
                    case 7: params.traitType = optarg; break;
                    case 8: params.invNormalize = true; break;
                    case 9: params.covarColList = split(optarg, ','); break;
                    case 10: params.qCovarCol = split(optarg, ','); break;
                    case 11: params.eCovarCol = split(optarg, ','); break;
                    case 12: params.sampleCovarCol = split(optarg, ','); break;
                    case 13: params.offsetCol = split(optarg, ','); break;
                    case 14: params.varWeightsCol = split(optarg, ','); break;
                    case 15: params.longlCol = optarg; break;
                    case 16: params.sampleIDColinphenoFile = optarg; break;
                    case 17: params.tol = std::stod(optarg); break;
                        if (isDouble(optarg)) {
                            params.tol = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for tol: " << optarg << std::endl;
                            return false;
                        }
                        break;
                        case 18: 
                        if (isInt(optarg)) {
                            params.maxiter = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for maxiter: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 19: 
                        if (isDouble(optarg)) {
                            params.tolPCG = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for tolPCG: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 20: 
                        if (isInt(optarg)) {
                            params.maxiterPCG = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for maxiterPCG: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 21: 
                        if (isInt(optarg)) {
                            params.nThreads = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for nThreads: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 22: 
                        if (isDouble(optarg)) {
                            params.SPAcutoff = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for SPAcutoff: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 23: 
                        if (isInt(optarg)) {
                            params.numMarkersForVarRatio = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for numMarkersForVarRatio: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 24: params.skipModelFitting = true; break;
                    case 25: 
                        if (isDouble(optarg)) {
                            params.memoryChunk = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for memoryChunk: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 26: 
                        if (isDoubleVector(optarg, ',')) {
                            params.tauInit = splitToDouble(optarg, ',');
                        } else {
                            std::cerr << "Invalid value for tauInit: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 27: params.LOCO = true; break;
                    case 28: params.isLowMemLOCO = true; break;
                    case 29: 
                        if (isDouble(optarg)) {
                            params.traceCVcutoff = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for traceCVcutoff: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 30: 
                        if (isDouble(optarg)) {
                            params.ratioCVcutoff = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for ratioCVcutoff: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 31: params.outputPrefix = optarg; break;
                    case 32: params.outputPrefix_varRatio = optarg; break;
                    case 33: params.IsOverwriteVarianceRatioFile = true; break;
                    case 34: params.sparseGRMFile = optarg; break;
                    case 35: params.sparseGRMSampleIDFile = optarg; break;
                    case 36: 
                        if (isInt(optarg)) {
                            params.numRandomMarkerforSparseKin = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for numRandomMarkerforSparseKin: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 37: 
                        if (isDouble(optarg)) {
                            params.relatednessCutoff = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for relatednessCutoff: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 38: params.isCateVarianceRatio = true; break;
                    case 39: 
                        if (isIntVector(optarg, ',')) {
                            params.cateVarRatioIndexVec = splitToInt(optarg, ',');
                        } else {
                            std::cerr << "Invalid value for cateVarRatioIndexVec: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 40: 
                        if (isDoubleVector(optarg, ',')) {
                            params.cateVarRatioMinMACVecExclude = splitToDouble(optarg, ',');
                        } else {
                            std::cerr << "Invalid value for cateVarRatioMinMACVecExclude: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 41: 
                        if (isDoubleVector(optarg, ',')) {
                            params.cateVarRatioMaxMACVecInclude = splitToDouble(optarg, ',');
                        } else {
                            std::cerr << "Invalid value for cateVarRatioMaxMACVecInclude: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 42: params.isCovariateTransform = true; break;
                    case 43: params.isDiagofKinSetAsOne = true; break;
                    case 44: 
                        if (isInt(optarg)) {
                            params.minCovariateCount = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for minCovariateCount: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 45: 
                        if (isDouble(optarg)) {
                            params.minMAFforGRM = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for minMAFforGRM: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 46: 
                        if (isDouble(optarg)) {
                            params.maxMissingRateforGRM = std::stod(optarg);
                        } else {
                            std::cerr << "Invalid value for maxMissingRateforGRM: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 47: params.useSparseGRMtoFitNULL = true; break;
                    case 48: params.useSparseGRMforVarRatio = true; break;
                    case 49: params.includeNonautoMarkersforVarRatio = true; break;
                    case 50: params.sexCol = optarg; break;
                    case 51: 
                        if (isInt(optarg)) {
                            params.FemaleCode = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for FemaleCode: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 52: params.FemaleOnly = true; break;
                    case 53: 
                        if (isInt(optarg)) {
                            params.MaleCode = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for MaleCode: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 54: params.MaleOnly = true; break;
                    case 55: params.SampleIDIncludeFile = optarg; break;
                    case 56: params.isCovariateOffset = true; break;
                    case 57: params.skipVarianceRatioEstimation = true; break;
                    case 58: 
                        if (isInt(optarg)) {
                            params.nrun = std::stoi(optarg);
                        } else {
                            std::cerr << "Invalid value for nrun: " << optarg << std::endl;
                            return false;
                        }
                        break;
                    case 59: params.VmatFilelist = optarg; break;
                    case 60: params.VmatSampleFilelist = optarg; break;
                    case 61: params.VcellmatFilelist = optarg; break;
                    case 62: params.VcellmatSampleFilelist = optarg; break;
                    case 63: params.useGRMtoFitNULL = true; break;
                    case 64: params.isStoreSigma = true; break;
                    case 65: params.isShrinkModelOutput = true; break;
                    default:
                        std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
                        return false;
                }
                break;
            default:
                std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
                return false;
        }
    }

    // Check if file parameters exist
    if (!params.plinkFile.empty()) {
        params.bedFile = params.plinkFile + ".bed";
        params.bimFile = params.plinkFile + ".bim";
        params.famFile = params.plinkFile + ".fam";
        if (!fileExists(params.bedFile)) {
            std::cerr << "File not found: " << params.bedFile << std::endl;
            return false;
        }
        if (!fileExists(params.bimFile)) {
            std::cerr << "File not found: " << params.bimFile << std::endl;
            return false;
        }
        if (!fileExists(params.famFile)) {
            std::cerr << "File not found: " << params.famFile << std::endl;
            return false;
        }
    }

    if (!params.phenoFile.empty() && !fileExists(params.phenoFile)) {
        std::cerr << "File not found: " << params.phenoFile << std::endl;
        return false;
    }
    if (!params.sparseGRMFile.empty()) {
        if (params.sparseGRMSampleIDFile.empty() || !fileExists(params.sparseGRMSampleIDFile)) {
            std::cerr << "sparseGRMSampleIDFile must be specified and exist when sparseGRMFile is specified" << std::endl;
            return false;
        }
        if (!fileExists(params.sparseGRMFile)) {
            std::cerr << "File not found: " << params.sparseGRMFile << std::endl;
            return false;
        }
    }
    if (params.FemaleOnly && params.MaleOnly) {
        std::cerr << "Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job" << std::endl;
        return false;
    }

    if (params.FemaleOnly) {
        params.outputPrefix += "_FemaleOnly";
        std::cout << "Female-specific model will be fitted. Samples coded as " << params.FemaleCode << " in the column " << params.sexCol << " in the phenotype file will be included" << std::endl;
    } else if (params.MaleOnly) {
        params.outputPrefix += "_MaleOnly";
        std::cout << "Male-specific model will be fitted. Samples coded as " << params.MaleCode << " in the column " << params.sexCol << " in the phenotype file will be included" << std::endl;
    }

    if ((!params.useSparseGRMtoFitNULL && params.useGRMtoFitNULL) || !params.skipVarianceRatioEstimation) {
        if (!fileExists(params.bedFile)) {
            std::cerr << "ERROR! bed file does not exist" << std::endl;
            return false;
        }
        if (!fileExists(params.bimFile)) {
            std::cerr << "ERROR! bim file does not exist" << std::endl;
            return false;
        }
    }

    return true;
}

bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}
    }
    return true;
}


bool parseStep2Parameters(int argc, char* argv[], Step2Parameters& params) {
    static struct option long_options[] = {
        {"vcfFile", required_argument, 0, 0},
        {"vcfFileIndex", required_argument, 0, 0},
        {"vcfField", required_argument, 0, 0},
        {"savFile", required_argument, 0, 0},
        {"savFileIndex", required_argument, 0, 0},
        {"bgenFile", required_argument, 0, 0},
        {"bgenFileIndex", required_argument, 0, 0},
        {"sampleFile", required_argument, 0, 0},
        {"bedFile", required_argument, 0, 0},
        {"bimFile", required_argument, 0, 0},
        {"famFile", required_argument, 0, 0},
        {"AlleleOrder", required_argument, 0, 0},
        {"idstoIncludeFile", required_argument, 0, 0},
        {"rangestoIncludeFile", required_argument, 0, 0},
        {"chrom", required_argument, 0, 0},
        {"is_imputed_data", no_argument, 0, 0},
        {"minMAF", required_argument, 0, 0},
        {"minMAC", required_argument, 0, 0},
        {"minGroupMAC_in_BurdenTest", required_argument, 0, 0},
        {"minInfo", required_argument, 0, 0},
        {"maxMissing", required_argument, 0, 0},
        {"impute_method", required_argument, 0, 0},
        {"LOCO", no_argument, 0, 0},
        {"GMMATmodelFile", required_argument, 0, 0},
        {"varianceRatioFile", required_argument, 0, 0},
        {"SAIGEOutputFile", required_argument, 0, 0},
        {"markers_per_chunk", required_argument, 0, 0},
        {"groups_per_chunk", required_argument, 0, 0},
        {"is_output_moreDetails", no_argument, 0, 0},
        {"is_overwrite_output", no_argument, 0, 0},
        {"maxMAF_in_groupTest", required_argument, 0, 0},
        {"maxMAC_in_groupTest", required_argument, 0, 0},
        {"annotation_in_groupTest", required_argument, 0, 0},
        {"groupFile", required_argument, 0, 0},
        {"sparseGRMFile", required_argument, 0, 0},
        {"sparseGRMSampleIDFile", required_argument, 0, 0},
        {"relatednessCutoff", required_argument, 0, 0},
        {"MACCutoff_to_CollapseUltraRare", required_argument, 0, 0},
        {"cateVarRatioMinMACVecExclude", required_argument, 0, 0},
        {"cateVarRatioMaxMACVecInclude", required_argument, 0, 0},
        {"weights.beta", required_argument, 0, 0},
        {"r.corr", required_argument, 0, 0},
        {"markers_per_chunk_in_groupTest", required_argument, 0, 0},
        {"condition", required_argument, 0, 0},
        {"weights_for_condition", required_argument, 0, 0},
        {"SPAcutoff", required_argument, 0, 0},
        {"dosage_zerod_cutoff", required_argument, 0, 0},
        {"dosage_zerod_MAC_cutoff", required_argument, 0, 0},
        {"is_single_in_groupTest", no_argument, 0, 0},
        {"is_no_weight_in_groupTest", no_argument, 0, 0},
        {"is_output_markerList_in_groupTest", no_argument, 0, 0},
        {"is_Firth_beta", no_argument, 0, 0},
        {"pCutoffforFirth", required_argument, 0, 0},
        {"is_fastTest", no_argument, 0, 0},
        {"max_MAC_for_ER", required_argument, 0, 0},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    while (true) {
        int c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
            switch (option_index) {
            case 0: params.vcfFile = optarg; break;
            case 1: params.vcfFileIndex = optarg; break;
            case 2: params.vcfField = optarg; break;
            case 3: params.savFile = optarg; break;
            case 4: params.savFileIndex = optarg; break;
            case 5: params.bgenFile = optarg; break;
            case 6: params.bgenFileIndex = optarg; break;
            case 7: params.sampleFile = optarg; break;
            case 8: params.bedFile = optarg; break;
            case 9: params.bimFile = optarg; break;
            case 10: params.famFile = optarg; break;
            case 11: params.AlleleOrder = optarg; break;
            case 12: params.idstoIncludeFile = optarg; break;
            case 13: params.rangestoIncludeFile = optarg; break;
            case 14: params.chrom = optarg; break;
            case 15: params.is_imputed_data = true; break;
            case 16: 
            if (isDouble(optarg)) {
                params.minMAF = std::stod(optarg);
            } else {
                std::cerr << "Invalid value for minMAF: " << optarg << std::endl;
                return false;
            }
            break;
            case 17: 
            if (isDouble(optarg)) {
                params.minMAC = std::stod(optarg);
            } else {
                std::cerr << "Invalid value for minMAC: " << optarg << std::endl;
                return false;
            }
            break;
            case 18: 
            if (isDouble(optarg)) {
                params.minGroupMAC_in_BurdenTest = std::stod(optarg);
            } else {
                std::cerr << "Invalid value for minGroupMAC_in_BurdenTest: " << optarg << std::endl;
                return false;
            }
            break;
            case 19: 
            if (isDouble(optarg)) {
                params.minInfo = std::stod(optarg);
            } else {
                std::cerr << "Invalid value for minInfo: " << optarg << std::endl;
                return false;
            }
            break;
            case 20: 
            if (isDouble(optarg)) {
                params.maxMissing = std::stod(optarg);
            } else {
                std::cerr << "Invalid value for maxMissing: " << optarg << std::endl;
                return false;
            }
            break;
            case 21: params.impute_method = optarg; break;
            case
                case 23: params.GMMATmodelFile = optarg; break;
                case 24: params.varianceRatioFile = optarg; break;
                case 25: params.SAIGEOutputFile = optarg; break;
                case 26: 
                if (isInt(optarg)) {
                    params.markers_per_chunk = std::stoi(optarg);
                } else {
                    std::cerr << "Invalid value for markers_per_chunk: " << optarg << std::endl;
                    return false;
                }
                break;
                case 27: 
                if (isInt(optarg)) {
                    params.groups_per_chunk = std::stoi(optarg);
                } else {
                    std::cerr << "Invalid value for groups_per_chunk: " << optarg << std::endl;
                    return false;
                }
                break;
                case 28: params.is_output_moreDetails = true; break;
                case 29: params.is_overwrite_output = true; break;
                case 30: params.maxMAF_in_groupTest = optarg; break;
                case 31: params.maxMAC_in_groupTest = optarg; break;
                case 32: params.annotation_in_groupTest = optarg; break;
                case 33: params.groupFile = optarg; break;
                case 34: params.sparseGRMFile = optarg; break;
                case 35: params.sparseGRMSampleIDFile = optarg; break;
                case 36: 
                if (isDouble(optarg)) {
                    params.relatednessCutoff = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for relatednessCutoff: " << optarg << std::endl;
                    return false;
                }
                break;
                case 37: 
                if (isDouble(optarg)) {
                    params.MACCutoff_to_CollapseUltraRare = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for MACCutoff_to_CollapseUltraRare: " << optarg << std::endl;
                    return false;
                }
                break;
                case 38: params.cateVarRatioMinMACVecExclude = optarg; break;
                case 39: params.cateVarRatioMaxMACVecInclude = optarg; break;
                case 40: params.weights_beta = optarg; break;
                case 41: 
                if (isDouble(optarg)) {
                    params.r_corr = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for r_corr: " << optarg << std::endl;
                    return false;
                }
                break;
                case 42: 
                if (isInt(optarg)) {
                    params.markers_per_chunk_in_groupTest = std::stoi(optarg);
                } else {
                    std::cerr << "Invalid value for markers_per_chunk_in_groupTest: " << optarg << std::endl;
                    return false;
                }
                break;
                case 43: params.condition = optarg; break;
                case 44: params.weights_for_condition = optarg; break;
                case 45: 
                if (isDouble(optarg)) {
                    params.SPAcutoff = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for SPAcutoff: " << optarg << std::endl;
                    return false;
                }
                break;
                case 46: 
                if (isDouble(optarg)) {
                    params.dosage_zerod_cutoff = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for dosage_zerod_cutoff: " << optarg << std::endl;
                    return false;
                }
                break;
                case 47: 
                if (isDouble(optarg)) {
                    params.dosage_zerod_MAC_cutoff = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for dosage_zerod_MAC_cutoff: " << optarg << std::endl;
                    return false;
                }
                break;
                case 48: params.is_single_in_groupTest = true; break;
                case 49: params.is_no_weight_in_groupTest = true; break;
                case 50: params.is_output_markerList_in_groupTest = true; break;
                case 51: params.is_Firth_beta = true; break;
                case 52: 
                if (isDouble(optarg)) {
                    params.pCutoffforFirth = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for pCutoffforFirth: " << optarg << std::endl;
                    return false;
                }
                break;
                case 53: params.is_fastTest = true; break;
                case 54: 
                if (isDouble(optarg)) {
                    params.max_MAC_for_ER = std::stod(optarg);
                } else {
                    std::cerr << "Invalid value for max_MAC_for_ER: " << optarg << std::endl;
                    return false;
                }
                break;
                default:
                std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
                return false;
            }
            break;
            default:
            std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
            return false;
        }break;
                    case 1: params.vcfFileIndex = optarg; break;
                    case 2: params.vcfField = optarg; break;
                    case 3: params.savFile = optarg; break;
                    case 4: params.savFileIndex = optarg; break;
                    case 5: params.bgenFile = optarg; break;
                    case 6: params.bgenFileIndex = optarg; break;
                    case 7: params.sampleFile = optarg; break;
                    case 8: params.bedFile = optarg; break;
                    case 9: params.bimFile = optarg; break;
                    case 10: params.famFile = optarg; break;
                    case 11: params.AlleleOrder = optarg; break;
                    case 12: params.idstoIncludeFile = optarg; break;
                    case 13: params.rangestoIncludeFile = optarg; break;
                    case 14: params.chrom = optarg; break;
                    case 15: params.is_imputed_data = true; break;
                    case 16: params.minMAF = std::stod(optarg); break;
                    case 17: params.minMAC = std::stod(optarg); break;
                    case 18: params.minGroupMAC_in_BurdenTest = std::stod(optarg); break;
                    case 19: params.minInfo = std::stod(optarg); break;
                    case 20: params.maxMissing = std::stod(optarg); break;
                    case 21: params.impute_method = optarg; break;
                    case 22: params.LOCO = true; break;
                    case 23: params.GMMATmodelFile = optarg; break;
                    case 24: params.varianceRatioFile = optarg; break;
                    case 25: params.SAIGEOutputFile = optarg; break;
                    case 26: params.markers_per_chunk = std::stoi(optarg); break;
                    case 27: params.groups_per_chunk = std::stoi(optarg); break;
                    case 28: params.is_output_moreDetails = true; break;
                    case 29: params.is_overwrite_output = true; break;
                    case 30: params.maxMAF_in_groupTest = optarg; break;
                    case 31: params.maxMAC_in_groupTest = optarg; break;
                    case 32: params.annotation_in_groupTest = optarg; break;
                    case 33: params.groupFile = optarg; break;
                    case 34: params.sparseGRMFile = optarg; break;
                    case 35: params.sparseGRMSampleIDFile = optarg; break;
                    case 36: params.relatednessCutoff = std::stod(optarg); break;
                    case 37: params.MACCutoff_to_CollapseUltraRare = std::stod(optarg); break;
                    case 38: params.cateVarRatioMinMACVecExclude = optarg; break;
                    case 39: params.cateVarRatioMaxMACVecInclude = optarg; break;
                    case 40: params.weights_beta = optarg; break;
                    case 41: params.r_corr = std::stod(optarg); break;
                    case 42: params.markers_per_chunk_in_groupTest = std::stoi(optarg); break;
                    case 43: params.condition = optarg; break;
                    case 44: params.weights_for_condition = optarg; break;
                    case 45: params.SPAcutoff = std::stod(optarg); break;
                    case 46: params.dosage_zerod_cutoff = std::stod(optarg); break;
                    case 47: params.dosage_zerod_MAC_cutoff = std::stod(optarg); break;
                    case 48: params.is_single_in_groupTest = true; break;
                    case 49: params.is_no_weight_in_groupTest = true; break;
                    case 50: params.is_output_markerList_in_groupTest = true; break;
                    case 51: params.is_Firth_beta = true; break;
                    case 52: params.pCutoffforFirth = std::stod(optarg); break;
                    case 53: params.is_fastTest = true; break;
                    case 54: params.max_MAC_for_ER = std::stod(optarg); break;
                }
                break;
            default:
                std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
                return false;
        }
    }
    return true;
}



void printParameters(const Parameters& params) {
    std::cout << "Parsed options:" << std::endl;
    if (params.step == "1") {
        std::cout << "plinkFile: " << params.plinkFile << std::endl;
        std::cout << "phenoFile: " << params.phenoFile << std::endl;
        std::cout << "phenoCol: " << params.phenoCol << std::endl;
        std::cout << "traitType: " << params.traitType << std::endl;
        std::cout << "invNormalize: " << params.invNormalize << std::endl;
        std::cout << "covarColList: ";
        for (const auto& col : params.covarColList) std::cout << col << " ";
        std::cout << std::endl;
        std::cout << "qCovarCol: ";
        for (const auto& col : params.qCovarCol) std::cout << col << " ";
        std::cout << std::endl;
        std::cout << "sampleIDColinphenoFile: " << params.sampleIDColinphenoFile << std::endl;
        std::cout << "tol: " << params.tol << std::endl;
        std::cout << "maxiter: " << params.maxiter << std::endl;
        std::cout << "tolPCG: " << params.tolPCG << std::endl;
        std::cout << "maxiterPCG: " << params.maxiterPCG << std::endl;
        std::cout << "nThreads: " << params.nThreads << std::endl;
        std::cout << "SPAcutoff: " << params.SPAcutoff << std::endl;
        std::cout << "numMarkersForVarRatio: " << params.numMarkersForVarRatio << std::endl;
        std::cout << "skipModelFitting: " << params.skipModelFitting << std::endl;
        std::cout << "memoryChunk: " << params.memoryChunk << std::endl;
        std::cout << "tauInit: ";
        for (const auto& tau : params.tauInit) std::cout << tau << " ";
        std::cout << std::endl;
        std::cout << "LOCO: " << params.LOCO << std::endl;
        std::cout << "traceCVcutoff: " << params.traceCVcutoff << std::endl;
        std::cout << "ratioCVcutoff: " << params.ratioCVcutoff << std::endl;
        std::cout << "outputPrefix: " << params.outputPrefix << std::endl;
        std::cout << "outputPrefix_varRatio: " << params.outputPrefix_varRatio << std::endl;
        std::cout << "IsOverwriteVarianceRatioFile: " << params.IsOverwriteVarianceRatioFile << std::endl;
        std::cout << "sparseGRMFile: " << params.sparseGRMFile << std::endl;
        std::cout << "sparseGRMSampleIDFile: " << params.sparseGRMSampleIDFile << std::endl;
        std::cout << "numRandomMarkerforSparseKin: " << params.numRandomMarkerforSparseKin << std::endl;
        std::cout << "relatednessCutoff: " << params.relatednessCutoff << std::endl;
        std::cout << "isCateVarianceRatio: " << params.isCateVarianceRatio << std::endl;
        std::cout << "cateVarRatioIndexVec: ";
        for (const auto& idx : params.cateVarRatioIndexVec) std::cout << idx << " ";
        std::cout << std::endl;
        std::cout << "cateVarRatioMinMACVecExclude: ";
        for (const auto& mac : params.cateVarRatioMinMACVecExclude) std::cout << mac << " ";
        std::cout << std::endl;
        std::cout << "cateVarRatioMaxMACVecInclude: ";
        for (const auto& mac : params.cateVarRatioMaxMACVecInclude) std::cout << mac << " ";
        std::cout << std::endl;
        std::cout << "isCovariateTransform: " << params.isCovariateTransform << std::endl;
        std::cout << "isDiagofKinSetAsOne: " << params.isDiagofKinSetAsOne << std::endl;
        std::cout << "useSparseGRMtoFitNULL: " << params.useSparseGRMtoFitNULL << std::endl;
        std::cout << "useSparseGRMforVarRatio: " << params.useSparseGRMforVarRatio << std::endl;
        std::cout << "includeNonautoMarkersforVarRatio: " << params.includeNonautoMarkersforVarRatio << std::endl;
        std::cout << "isStoreSigma: " << params.isStoreSigma << std::endl;
        std::cout << "isShrinkModelOutput: " << params.isShrinkModelOutput << std::endl;
    } else if (params.step == "2") {
        std::cout << "vcfFile: " << params.vcfFile << std::endl;
        std::cout << "vcfFileIndex: " << params.vcfFileIndex << std::endl;
        std::cout << "vcfField: " << params.vcfField << std::endl;
        std::cout << "savFile: " << params.savFile << std::endl;
        std::cout << "savFileIndex: " << params.savFileIndex << std::endl;
        std::cout << "bgenFile: " << params.bgenFile << std::endl;
        std::cout << "bgenFileIndex: " << params.bgenFileIndex << std::endl;
        std::cout << "sampleFile: " << params.sampleFile << std::endl;
        std::cout << "bedFile: " << params.bedFile << std::endl;
        std::cout << "bimFile: " << params.bimFile << std::endl;
        std::cout << "famFile: " << params.famFile << std::endl;
        std::cout << "AlleleOrder: " << params.AlleleOrder << std::endl;
        std::cout << "idstoIncludeFile: " << params.idstoIncludeFile << std::endl;
        std::cout << "rangestoIncludeFile: " << params.rangestoIncludeFile << std::endl;
        std::cout << "chrom: " << params.chrom << std::endl;
        std::cout << "is_imputed_data: " << params.is_imputed_data << std::endl;
        std::cout << "minMAF: " << params.minMAF << std::endl;
        std::cout << "minMAC: " << params.minMAC << std::endl;
        std::cout << "minGroupMAC_in_BurdenTest: " << params.minGroupMAC_in_BurdenTest << std::endl;
        std::cout << "minInfo: " << params.minInfo << std::endl;
        std::cout << "maxMissing: " << params.maxMissing << std::endl;
        std::cout << "impute_method: " << params.impute_method << std::endl;
        std::cout << "LOCO: " << params.LOCO << std::endl;
        std::cout << "GMMATmodelFile: " << params.GMMATmodelFile << std::endl;
        std::cout << "varianceRatioFile: " << params.varianceRatioFile << std::endl;
        std::cout << "SAIGEOutputFile: " << params.SAIGEOutputFile << std::endl;
        std::cout << "markers_per_chunk: " << params.markers_per_chunk << std::endl;
        std::cout << "groups_per_chunk: " << params.groups_per_chunk << std::endl;
        std::cout << "is_output_moreDetails: " << params.is_output_moreDetails << std::endl;
        std::cout << "is_overwrite_output: " << params.is_overwrite_output << std::endl;
        std::cout << "maxMAF_in_groupTest: " << params.maxMAF_in_groupTest << std::endl;
        std::cout << "maxMAC_in_groupTest: " << params.maxMAC_in_groupTest << std::endl;
        std::cout << "annotation_in_groupTest: " << params.annotation_in_groupTest << std::endl;
        std::cout << "groupFile: " << params.groupFile << std::endl;
        std::cout << "sparseGRMFile: " << params.sparseGRMFile << std::endl;
        std::cout << "sparseGRMSampleIDFile: " << params.sparseGRMSampleIDFile << std::endl;
        std::cout << "relatednessCutoff: " << params.relatednessCutoff << std::endl;
        std::cout << "MACCutoff_to_CollapseUltraRare: " << params.MACCutoff_to_CollapseUltraRare << std::endl;
        std::cout << "cateVarRatioMinMACVecExclude: " << params.cateVarRatioMinMACVecExclude << std::endl;
        std::cout << "cateVarRatioMaxMACVecInclude: " << params.cateVarRatioMaxMACVecInclude << std::endl;
        std::cout << "weights.beta: " << params.weights_beta << std::endl;
        std::cout << "r.corr: " << params.r_corr << std::endl;
        std::cout << "markers_per_chunk_in_groupTest: " << params.markers_per_chunk_in_groupTest << std::endl;
        std::cout << "condition: " << params.condition << std::endl;
        std::cout << "weights_for_condition: " << params.weights_for_condition << std::endl;
        std::cout << "SPAcutoff: " << params.SPAcutoff << std::endl;
        std::cout << "dosage_zerod_cutoff: " << params.dosage_zerod_cutoff << std::endl;
        std::cout << "dosage_zerod_MAC_cutoff: " << params.dosage_zerod_MAC_cutoff << std::endl;
        std::cout << "is_single_in_groupTest: " << params.is_single_in_groupTest << std::endl;
        std::cout << "is_no_weight_in_groupTest: " << params.is_no_weight_in_groupTest << std::endl;
        std::cout << "is_output_markerList_in_groupTest: " << params.is_output_markerList_in_groupTest << std::endl;
        std::cout << "is_Firth_beta: " << params.is_Firth_beta << std::endl;
        std::cout << "pCutoffforFirth: " << params.pCutoffforFirth << std::endl;
        std::cout << "is_fastTest: " << params.is_fastTest << std::endl;
        std::cout << "max_MAC_for_ER: " << params.max_MAC_for_ER << std::endl;
    } else {
        std::cerr << "Unknown step: " << params.step << std::endl;
    }
}

void parseArguments(int argc, char* argv[], Parameters& params) {
    static struct option long_options[] = {
        {"plinkFile", required_argument, 0, 0},
        {"phenoFile", required_argument, 0, 0},
        {"phenoCol", required_argument, 0, 0},
        {"traitType", required_argument, 0, 0},
        {"invNormalize", no_argument, 0, 0},
        {"covarColList", required_argument, 0, 0},
        {"qCovarCol", required_argument, 0, 0},
        {"sampleIDColinphenoFile", required_argument, 0, 0},
        {"tol", required_argument, 0, 0},
        {"maxiter", required_argument, 0, 0},
        {"tolPCG", required_argument, 0, 0},
        {"maxiterPCG", required_argument, 0, 0},
        {"nThreads", required_argument, 0, 0},
        {"SPAcutoff", required_argument, 0, 0},
        {"numMarkersForVarRatio", required_argument, 0, 0},
        {"skipModelFitting", no_argument, 0, 0},
        {"memoryChunk", required_argument, 0, 0},
        {"tauInit", required_argument, 0, 0},
        {"LOCO", no_argument, 0, 0},
        {"traceCVcutoff", required_argument, 0, 0},
        {"ratioCVcutoff", required_argument, 0, 0},
        {"outputPrefix", required_argument, 0, 0},
        {"outputPrefix_varRatio", required_argument, 0, 0},
        {"IsOverwriteVarianceRatioFile", no_argument, 0, 0},
        {"sparseGRMFile", required_argument, 0, 0},
        {"sparseGRMSampleIDFile", required_argument, 0, 0},
        {"numRandomMarkerforSparseKin", required_argument, 0, 0},
        {"relatednessCutoff", required_argument, 0, 0},
        {"isCateVarianceRatio", no_argument, 0, 0},
        {"cateVarRatioIndexVec", required_argument, 0, 0},
        {"cateVarRatioMinMACVecExclude", required_argument, 0, 0},
        {"cateVarRatioMaxMACVecInclude", required_argument, 0, 0},
        {"isCovariateTransform", no_argument, 0, 0},
        {"isDiagofKinSetAsOne", no_argument, 0, 0},
        {"useSparseGRMtoFitNULL", no_argument, 0, 0},
        {"useSparseGRMforVarRatio", no_argument, 0, 0},
        {"includeNonautoMarkersforVarRatio", no_argument, 0, 0},
        {"isStoreSigma", no_argument, 0, 0},
        {"isShrinkModelOutput", no_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    while (true) {
        int c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0:
                switch (option_index) {
                    case 0: params.plinkFile = optarg; break;
                    case 1: params.phenoFile = optarg; break;
                    case 2: params.phenoCol = optarg; break;
                    case 3: params.traitType = optarg; break;
                    case 4: params.invNormalize = true; break;
                    case 5: params.covarColList = split(optarg, ','); break;
                    case 6: params.qCovarCol = split(optarg, ','); break;
                    case 7: params.sampleIDColinphenoFile = optarg; break;
                    case 8: params.tol = std::stod(optarg); break;
                    case 9: params.maxiter = std::stoi(optarg); break;
                    case 10: params.tolPCG = std::stod(optarg); break;
                    case 11: params.maxiterPCG = std::stoi(optarg); break;
                    case 12: params.nThreads = std::stoi(optarg); break;
                    case 13: params.SPAcutoff = std::stod(optarg); break;
                    case 14: params.numMarkersForVarRatio = std::stoi(optarg); break;
                    case 15: params.skipModelFitting = true; break;
                    case 16: params.memoryChunk = std::stod(optarg); break;
                    case 17: params.tauInit = splitToDouble(optarg, ','); break;
                    case 18: params.LOCO = true; break;
                    case 19: params.traceCVcutoff = std::stod(optarg); break;
                    case 20: params.ratioCVcutoff = std::stod(optarg); break;
                    case 21: params.outputPrefix = optarg; break;
                    case 22: params.outputPrefix_varRatio = optarg; break;
                    case 23: params.IsOverwriteVarianceRatioFile = true; break;
                    case 24: params.sparseGRMFile = optarg; break;
                    case 25: params.sparseGRMSampleIDFile = optarg; break;
                    case 26: params.numRandomMarkerforSparseKin = std::stoi(optarg); break;
                    case 27: params.relatednessCutoff = std::stod(optarg); break;
                    case 28: params.isCateVarianceRatio = true; break;
                    case 29: params.cateVarRatioIndexVec = splitToInt(optarg, ','); break;
                    case 30: params.cateVarRatioMinMACVecExclude = splitToDouble(optarg, ','); break;
                    case 31: params.cateVarRatioMaxMACVecInclude = splitToDouble(optarg, ','); break;
                    case 32: params.isCovariateTransform = true; break;
                    case 33: params.isDiagofKinSetAsOne = true; break;
                    case 34: params.useSparseGRMtoFitNULL = true; break;
                    case 35: params.useSparseGRMforVarRatio = true; break;
                    case 36: params.includeNonautoMarkersforVarRatio = true; break;
                    case 37: params.isStoreSigma = true; break;
                    case 38: params.isShrinkModelOutput = true; break;
                }
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
            default:
                printHelp();
                exit(EXIT_FAILURE);
        }
    }

    // Print out the parsed options
    printParameters(params);
}

void checkFiles(const std::vector<std::string>& filePaths) {
    for (const auto& filePath : filePaths) {
        if (!filePath.empty()) {
            if (!std::filesystem::exists(filePath)) {
                throw std::runtime_error("ERROR! File does not exist: " + filePath);
            }
        }
    }
}
void preprocessData(const Parameters& params) {
    std::cout << "Preprocessing data" << std::endl;

    // Check if the bed, bim, and fam files exist
    checkFiles({params.bedFile, params.bimFile, params.famFile});

    // Read sample IDs from the fam file
    std::ifstream famFile(params.famFile);
    if (!famFile.is_open()) {
        throw std::runtime_error("ERROR! Unable to open fam file");
    }

    std::vector<std::string> sampleIDs;
    std::string line;
    while (std::getline(famFile, line)) {
        std::istringstream iss(line);
        std::string id;
        for (int i = 0; i < 2; ++i) {
            iss >> id;
        }
        sampleIDs.push_back(id);
    }
    famFile.close();

    std::cout << sampleIDs.size() << " samples have genotypes" << std::endl;

    // Check if sparse GRM files exist and read sample IDs if required
    if (params.useSparseGRMtoFitNULL || params.useSparseGRMforVarRatio) {
        checkFiles({params.sparseGRMSampleIDFile});

        std::ifstream sparseGRMFile(params.sparseGRMSampleIDFile);
        if (!sparseGRMFile.is_open()) {
            throw std::runtime_error("ERROR! Unable to open sparse GRM sample ID file");
        }

        std::vector<std::string> sparseSampleIDs;
        while (std::getline(sparseGRMFile, line)) {
            sparseSampleIDs.push_back(line);
        }
        sparseGRMFile.close();

        std::cout << sparseSampleIDs.size() << " samples are in the sparse GRM" << std::endl;
    }
}


