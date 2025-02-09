#include "parameters.h"
#include "utils.h"
#include <iostream>
#include <cstring>

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
    std::string phenoFile;
    std::string phenoCol;
    std::string traitType;
    bool invNormalize;
    std::vector<std::string> covarColList;
    std::vector<std::string> qCovarCol;
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
    bool useSparseGRMtoFitNULL;
    bool useSparseGRMforVarRatio;
    bool includeNonautoMarkersforVarRatio;
    bool isStoreSigma;
    bool isShrinkModelOutput;
};

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
                    params.isShrinkModelOutput = true; break;
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
                    if (FemaleOnly && MaleOnly) {
                        std::cerr << "Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job" << std::endl;
                        return false;
                    }

                    if (FemaleOnly) {
                        params.outputPrefix += "_FemaleOnly";
                        std::cout << "Female-specific model will be fitted. Samples coded as " << FemaleCode << " in the column " << sexCol << " in the phenotype file will be included" << std::endl;
                    } else if (MaleOnly) {
                        params.outputPrefix += "_MaleOnly";
                        std::cout << "Male-specific model will be fitted. Samples coded as " << MaleCode << " in the column " << sexCol << " in the phenotype file will be included" << std::endl;
                    }
        if ((!params.useSparseGRMtoFitNULL && params.useGRMtoFitNULL) || !params.skipVarianceRatioEstimation) {
            if (!fileExists(params.bedFile)) {
                std::cerr << "ERROR! bed file does not exist" << std::endl;
                return false;
            }
            if (!fileExists(params.bimFile)) {
                std::cerr << "ERROR! bim file does not exist" << std::endl;
                return false;
            } else {
                if (params.LOCO) {
                    std::vector<std::string> chrVec = readBimFile(params.bimFile);
                    auto updatechrList = updateChrStartEndIndexVec(chrVec);
                    params.LOCO = updatechrList.LOCO;
                    params.chromosomeStartIndexVec = updatechrList.chromosomeStartIndexVec;
                    params.chromosomeEndIndexVec = updatechrList.chromosomeEndIndexVec;
                }
                if (!params.LOCO) {
                    params.chromosomeStartIndexVec.assign(22, NA);
                    params.chromosomeEndIndexVec.assign(22, NA);
                }
            }
        }
                    return true;
                    }

                    bool fileExists(const std::string& filename) {
                        std::ifstream file(filename);
                        return file.good();
                    }
                std::cerr << "Unknown parameter: " << argv[optind - 1] << std::endl;
                return false;
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
    if (!std::filesystem::exists(params.bedFile)) {
        throw std::runtime_error("ERROR! bed file does not exist");
    }
    if (!std::filesystem::exists(params.bimFile)) {
        throw std::runtime_error("ERROR! bim file does not exist");
    }
    if (!std::filesystem::exists(params.famFile)) {
        throw std::runtime_error("ERROR! fam file does not exist");
    }


    if (!fileExists(famFile)) {
        stop("ERROR! fam file does not exist\n");
      } else {
        if (famFile.empty()) {
        stop("ERROR! fam file is not specified\n");
        }
        std::vector<std::string> sampleListwithGenov0 = readFile(famFile);
        std::vector<std::string> IIDgeno;
        for (const auto& line : sampleListwithGenov0) {
        std::vector<std::string> tokens = split(line, ' ');
        IIDgeno.push_back(tokens[1]);
        }
        std::vector<int> IndexGeno(IIDgeno.size());
        std::iota(IndexGeno.begin(), IndexGeno.end(), 1);
        std::cout << IIDgeno.size() << " samples have genotypes\n";
      }

      
      } else {
      if (useSparseGRMtoFitNULL || useSparseGRMforVarRatio) {
        std::vector<std::string> sampleListwithGenov0 = readFile(sparseGRMSampleIDFile);
        std::vector<std::string> IIDgeno;
        for (const auto& line : sampleListwithGenov0) {
        IIDgeno.push_back(line);
        }
        std::vector<int> IndexGeno(IIDgeno.size());
        std::iota(IndexGeno.begin(), IndexGeno.end(), 1);
        std::cout << IIDgeno.size() << " samples are in the sparse GRM\n";
      }
      }



    // Read only the second column (sample IDs) from the fam file and store in an arma::mat
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

    // Convert sampleIDs to arma::mat
    arma::mat famData(sampleIDs.size(), 1);
    for (size_t i = 0; i < sampleIDs.size(); ++i) {
        famData(i, 0) = std::stod(sampleIDs[i]);
    }

    arma::mat famData;
    famData.load(params.famFile, arma::csv_ascii);

    // Extract sample IDs
    std::vector<std::string> sampleIDs;
    for (size_t i = 0; i < famData.n_rows; ++i) {
        sampleIDs.push_back(famData(i, 1));
    }

    // Print the number of samples with genotypes
    std::cout << sampleIDs.size() << " samples have genotypes" << std::endl;
}