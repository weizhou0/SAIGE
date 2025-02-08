#include <iostream>
#include "utils/parameters.h"
#include "grm/grm.h"
#include "glmm/glmm.h"
#include "association/single_variant.h"
#include "association/set_based.h"
#include "ld/ld_matrix.h"

int main(int argc, char* argv[]) {
    Parameters params;
    if (!params.parse(argc, argv)) {
        return 1;
    }

    if (params.step == "step0") {
        createSparseGRM(params);
    } else if (params.step == "step1") {
        fitNullGLMM(params);
    } else if (params.step == "step2") {
        if (params.testType == "single") {
            conductSingleVariantTest(params);
        } else if (params.testType == "set") {
            conductSetBasedTest(params);
        }
    } else if (params.step == "step3") {
        writeLDMatrix(params);
    } else {
        std::cerr << "Invalid step specified." << std::endl;
        return 1;
    }

    return 0;
}
