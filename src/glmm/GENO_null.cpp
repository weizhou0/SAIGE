#define ARMA_USE_SUPERLU 1

#include <armadillo>
#include <omp.h>
#include <tbb.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>// include this header for calculating execution time
#include <cassert>
#include <boost/date_time.hpp> // for gettimeofday and timeval
#include "getMem.hpp"
#include "GENO_null.hpp"


namespace NullGENO {

        NullGenoClass::NullGenoClass(){	
        }	

        void NullGenoClass::setStdGenoLookUpArr(double mafVal, double invsdVal, arma::vec & stdGenoLookUpArr) {
                double mafVal2 = 2 * mafVal;
                stdGenoLookUpArr(0) = (0 - mafVal2) * invsdVal;
                stdGenoLookUpArr(1) = (1 - mafVal2) * invsdVal;
                stdGenoLookUpArr(2) = (2 - mafVal2) * invsdVal;
        }



        void NullGenoClass::setSparseKinLookUpArr(double mafVal, double invsdVal) {
                double mafVal2 = 2 * mafVal;
                double a0 = (0 - mafVal2) * invsdVal;
                double a1 = (1 - mafVal2) * invsdVal;
                double a2 = (2 - mafVal2) * invsdVal;

                sKinLookUpArr[0][0] = a0 * a0;
                sKinLookUpArr[0][1] = a0 * a1;
                sKinLookUpArr[0][2] = a0 * a2;
                sKinLookUpArr[1][0] = sKinLookUpArr[0][1];
                sKinLookUpArr[1][1] = a1 * a1;
                sKinLookUpArr[1][2] = a1 * a2;
                sKinLookUpArr[2][0] = sKinLookUpArr[0][2];
                sKinLookUpArr[2][1] = sKinLookUpArr[1][2];
                sKinLookUpArr[2][2] = a2 * a2;
        }

        void NullGenoClass::setBit(unsigned char & ch, int ii, int aVal, int bVal) {
                if (bVal == 1 && aVal == 1) {
                        ch ^= char(1 << ((ii * 2) + 1)); // set a to be 1
                } else if (bVal == 0) {
                        ch ^= char(1 << (ii * 2)); // change b to 0
                        if (aVal == 1) {
                                ch ^= char(1 << ((ii * 2) + 1)); // change a to 1
                        }
                }
        }


        void NullGenoClass::setGenotype(unsigned char* c, const int pos, const int geno) {
                (*c) |= (geno << (pos << 1));
        }

        void NullGenoClass::getGenotype(unsigned char* c, const int pos, int& geno) {
                geno = ((*c) >> (pos << 1)) & 0x3;  // 0b11 = 0x3
        }

        void NullGenoClass::Init_OneSNP_Geno() {
                m_size_of_esi = (Nnomissing + 3) / 4;
                int k = 8;
                while (k > 0) {
                        --k;
                        m_bits_val[k] = 1 << k;
                }
        }

        arma::ivec* NullGenoClass::Get_OneSNP_Geno(size_t SNPIdx) {
                m_OneSNP_Geno.zeros(Nnomissing);

                // Avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx / numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind = 0;
                unsigned char geno1;
                int bufferGeno;

                for (size_t i = Start_idx; i < Start_idx + m_size_of_esi - 1; i++) {
                        geno1 = genoVecofPointers[indexOfVectorPointer]->at(i); // Avoid large continuous memory usage
                        for (int j = 0; j < 4; j++) {
                                int b = geno1 & 1;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1;
                                bufferGeno = 2 - (a + b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                        }
                }

                size_t i = Start_idx + m_size_of_esi - 1;
                geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);
                for (int j = 0; j < 4; j++) {
                        int b = geno1 & 1;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1;
                        bufferGeno = 2 - (a + b);
                        m_OneSNP_Geno[ind] = bufferGeno;
                        ind++;
                        geno1 = geno1 >> 1;
                        if (ind >= Nnomissing) {
                                return &m_OneSNP_Geno;
                        }
                }

                return &m_OneSNP_Geno;
        }



        arma::ivec* NullGenoClass::Get_OneSNP_Geno_forVarRatio(size_t SNPIdx) {
                m_OneSNP_Geno.zeros(Nnomissing);
                // Avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx / numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind = 0;
                unsigned char geno1;
                int bufferGeno;

                for (size_t i = Start_idx; i < Start_idx + m_size_of_esi - 1; i++) {
                        geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i); // Avoid large continuous memory usage
                        for (int j = 0; j < 4; j++) {
                                int b = geno1 & 1;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1;
                                bufferGeno = 2 - (a + b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                        }
                }

                size_t i = Start_idx + m_size_of_esi - 1;
                geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i);
                for (int j = 0; j < 4; j++) {
                        int b = geno1 & 1;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1;
                        bufferGeno = 2 - (a + b);
                        m_OneSNP_Geno[ind] = bufferGeno;
                        ind++;
                        geno1 = geno1 >> 1;
                        if (ind >= Nnomissing) {
                                return &m_OneSNP_Geno;
                        }
                }

                return &m_OneSNP_Geno;
        }



        void NullGenoClass::Get_OneSNP_Geno_atBeginning(
                size_t SNPIdx, 
                std::vector<int> & indexNA, 
                std::vector<unsigned char> & genoVecOneMarkerOld, 
                double & altFreq, 
                double & missingRate, 
                int & mac,  
                int & alleleCount, 
                bool & passQC, 
                size_t SNPIdx_new, 
                bool & passVarRatio, 
                size_t SNPIdx_vr
        ) {
                arma::ivec m_OneSNP_GenoTemp;
                m_OneSNP_GenoTemp.zeros(N);
                m_OneSNP_Geno.zeros(Nnomissing);
                int m_size_of_esi_temp = (N+3)/4;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                int u;
                alleleCount = 0;
                int numMissing = 0;
                for(int i=0; i< (m_size_of_esi_temp - 1); i++){
                        geno1 = genoVecOneMarkerOld[i];
                        for(int j=0; j<4; j++){
                                u = j & 3;
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                if (b == 1 && a == 0){
                                        bufferGeno = 3;
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;
                                }else{
                                        std::cout << "Error GENO!!\n";
                                        break;
                                }
                                m_OneSNP_GenoTemp[ind] = bufferGeno;
                                if(indicatorGenoSamplesWithPheno_in[ind]){
                                        if(bufferGeno == 3){
                                                numMissing++;
                                        }else{
                                                alleleCount = alleleCount + bufferGeno;
                                        }
                                }
                                ind++;
                                geno1 = geno1 >> 1;
                        }
                }
                int i = m_size_of_esi_temp - 1;
                geno1 = genoVecOneMarkerOld[i];
                for(int j=0; j<4; j++){
                        u = j & 3;
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        if (b == 1 && a == 0){
                                bufferGeno = 3;
                        }else if(b == 0 && a == 0){
                                bufferGeno = 2;
                        }else if(b == 0 && a == 1){
                                bufferGeno = 1;
                        }else if(b == 1 && a == 1){
                                bufferGeno = 0;
                        }else{
                                std::cout << "Error GENO!!\n";
                                break;
                        }
                        m_OneSNP_GenoTemp[ind] = bufferGeno;
                        if(indicatorGenoSamplesWithPheno_in[ind]){
                                if(bufferGeno == 3){
                                        numMissing++;
                                }else{
                                        alleleCount = alleleCount + bufferGeno;
                                }
                        }
                        ind++;
                        geno1 = geno1 >> 1;
                        if(ind == (N)){
                                break;
                        }
                }
                altFreq = alleleCount/((Nnomissing-numMissing) * 2.0);
                missingRate = numMissing/static_cast<double>(Nnomissing);
                int fillinMissingGeno = static_cast<int>(round(2*altFreq));
                if(numMissing > 0){
                        alleleCount = alleleCount + fillinMissingGeno*numMissing;
                }
                altFreq = alleleCount/(Nnomissing * 2.0);
                unsigned char geno2;
                passQC = false;
                passVarRatio = false;
                double maf = std::min(altFreq, 1-altFreq);
                mac = std::min(alleleCount, int(Nnomissing) * 2 - alleleCount);
                if(maf >= minMAFtoConstructGRM && missingRate <= maxMissingRate){
                        passQC = true;
                }
                if(isVarRatio){
                        if(g_maxMACVarRatio != -1){
                                if(mac >= g_minMACVarRatio && mac < g_maxMACVarRatio){
                                        passVarRatio = true;
                                        genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
                                        genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(double(Nnomissing)/4));
                                }else if(mac >= g_maxMACVarRatio){
                                        passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        if(passVarRatio){
                                                genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(double(Nnomissing)/4));
                                                numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
                                        }
                                }
                        }else{
                                if(mac >= g_minMACVarRatio){
                                        passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        if(passVarRatio){
                                                genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(double(Nnomissing)/4));
                                                numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
                                        }
                                }
                        }
                        if(passVarRatio){
                                passQC = false;
                        }
                }
                if(passQC || passVarRatio){
                        for(int indx=0; indx < Nnomissing; indx++){
                                u = indx & 3;
                                bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
                                if(bufferGeno == 3){
                                        bufferGeno = fillinMissingGeno;
                                }
                                if(bufferGeno == 0){
                                        setGenotype(&geno2, u, HOM_ALT);
                                }else if(bufferGeno == 1){
                                        setGenotype(&geno2, u, HET);
                                }else if(bufferGeno == 2){
                                        setGenotype(&geno2, u, HOM_REF);
                                }
                                if(u == 3 || indx == (Nnomissing-1)){
                                        if(passVarRatio){
                                                genoVecofPointers_forVarRatio[SNPIdx_vr/numMarkersofEachArray]->push_back(geno2);
                                        }
                                        if(passQC){
                                                genoVecofPointers[SNPIdx_new/numMarkersofEachArray]->push_back(geno2);
                                        }
                                        geno2 = 0;
                                }
                        }
                }
        }

        int NullGenoClass::Get_OneSNP_StdGeno(size_t SNPIdx, arma::vec * out) {
                // Avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx / numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;

                out->zeros(Nnomissing);

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind = 0;
                unsigned char geno1;

                double freq = alleleFreqVec[SNPIdx];
                double invStd = invstdvVec[SNPIdx];

                arma::vec stdGenoLookUpArr(3);
                setStdGenoLookUpArr(freq, invStd, stdGenoLookUpArr);

                for (size_t i = Start_idx; i < Start_idx + m_size_of_esi - 1; i++) {
                        geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

                        for (int j = 0; j < 4; j++) {
                                int b = geno1 & 1;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1;
                                (*out)[ind] = stdGenoLookUpArr(2 - (a + b));
                                ind++;
                                geno1 = geno1 >> 1;
                        }
                }

                size_t i = Start_idx + m_size_of_esi - 1;
                geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

                for (int j = 0; j < 4; j++) {
                        int b = geno1 & 1;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1;
                        (*out)[ind] = stdGenoLookUpArr(2 - (a + b));
                        ind++;
                        geno1 = geno1 >> 1;

                        if (ind >= Nnomissing) {
                                stdGenoLookUpArr.clear();
                                return 1;
                        }
                }

                stdGenoLookUpArr.clear();
                return 1;
        }


        arma::vec* NullGenoClass::Get_Diagof_StdGeno() {
                arma::vec* temp = &m_OneSNP_StdGeno;
                // Not yet calculated
                if (size(m_DiagStd)[0] != Nnomissing) {
                        m_DiagStd.zeros(Nnomissing);
                        for (size_t i = 0; i < numberofMarkerswithMAFge_minMAFtoConstructGRM; i++) {
                                Get_OneSNP_StdGeno(i, temp);
                                m_DiagStd = m_DiagStd + (*temp) % (*temp);
                        }
                }
                return &m_DiagStd;
        }


        arma::vec* NullGenoClass::Get_Diagof_StdGeno_LOCO() {

                m_DiagStd_LOCO = mtx_DiagStd_LOCO.col(chromIndex);
                Msub_MAFge_minMAFtoConstructGRM_singleChr = Msub_MAFge_minMAFtoConstructGRM_byChr(chromIndex);
                // }

                return &m_DiagStd_LOCO;
        }


        void NullGenoClass::setGenoObj(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, double memoryChunk, bool isDiagofKinSetAsOne) {
                setKinDiagtoOne = isDiagofKinSetAsOne;
                ptrsubSampleInGeno = subSampleInGeno;
                indicatorGenoSamplesWithPheno_in = indicatorGenoSamplesWithPheno;
                Nnomissing = subSampleInGeno.size();

                alleleFreqVec.clear();
                MACVec.clear();
                invstdvVec.clear();

                M = 0;
                N = 0;

                std::string junk;

                // Count the number of individuals
                ifstream test_famfile(famfile.c_str());
                if (!test_famfile.is_open()) {
                        printf("Error! fam file not open!");
                        return;
                }
                int indexRow = 0;
                while (std::getline(test_famfile, junk)) {
                        indexRow++;
                        junk.clear();
                }
                N = indexRow;
                test_famfile.clear();

                // Count the number of markers
                ifstream test_bimfile(bimfile.c_str());
                if (!test_bimfile.is_open()) {
                        printf("Error! bim file not open!");
                        return;
                }
                indexRow = 0;
                while (std::getline(test_bimfile, junk)) {
                        indexRow++;
                        junk.clear();
                }
                M = indexRow;
                test_bimfile.clear();

                // Init OneSNP Geno
                Init_OneSNP_Geno();

                int buffer;
                int TotalRead = 0;

                std::vector<unsigned char> genoVecOneMarkerOld;
                std::vector<unsigned char> genoVecOneMarkerNew;

                size_t nbyteOld = ceil(double(N) / 4);
                size_t nbyteNew = ceil(double(Nnomissing) / 4);
                size_t reserve = ceil(double(Nnomissing) / 4) * M + M * 2;


                genoVecOneMarkerOld.reserve(nbyteOld);
                genoVecOneMarkerOld.resize(nbyteOld);

                ifstream test_bedfile(bedfile.c_str(), ios::binary);
                if (!test_bedfile.is_open()) {
                        printf("Error! file open!");
                        return;
                }

                numMarkersofEachArray = 1;
                numofGenoArray = M;
                genoVecofPointers.resize(numofGenoArray);
                genoVecofPointers_forVarRatio.resize(numofGenoArray);

                for (int i = 0; i < numofGenoArray; i++) {
                        genoVecofPointers[i] = new vector<unsigned char>;
                        genoVecofPointers[i]->reserve(numMarkersofEachArray * ceil(double(Nnomissing) / 4));
                }

                arma::ivec g_randMarkerIndforVR_temp;
                if (isVarRatio) {
                        g_randMarkerIndforVR_temp = arma::randi(1000, arma::distr_param(0, M - 1));
                        g_randMarkerIndforVR = arma::unique(g_randMarkerIndforVR_temp);
                }

                double freq, Std, invStd, missingRate;
                int alleleCount, mac;
                std::vector<int> indexNA;
                int lengthIndexNA;
                int indexGeno;
                int indexBit;
                int fillinMissingGeno;
                int b2;
                int a2;

                size_t ind = 0;
                unsigned char geno1 = 0;
                int bufferGeno;
                int u;

                bool isPassQC = false;
                bool isPass_vr = false;
                size_t SNPIdx_new = 0;
                size_t SNPIdx_vr = 0;

                for (int i = 0; i < M; i++) {
                        genoVecOneMarkerOld.clear();
                        genoVecOneMarkerOld.reserve(nbyteOld);
                        genoVecOneMarkerOld.resize(nbyteOld);

                        test_bedfile.seekg(3 + nbyteOld * i);
                        test_bedfile.read((char*)(&genoVecOneMarkerOld[0]), nbyteOld);

                        indexNA.clear();

                        Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld, freq, missingRate, mac, alleleCount, isPassQC, SNPIdx_new, isPass_vr, SNPIdx_vr);

                        if (isPassQC) {
                                Std = std::sqrt(2 * freq * (1 - freq));
                                invStd = (Std == 0) ? 0 : 1 / Std;
                                invstdvVec0.push_back(invStd);
                                alleleFreqVec0.push_back(freq);
                                numberofMarkerswithMAFge_minMAFtoConstructGRM++;
                                MACVec0.push_back(mac);
                                MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(true);
                                SNPIdx_new++;
                        } else {
                                MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(false);
                        }

                        if (isVarRatio && isPass_vr) {
                                Std = std::sqrt(2 * freq * (1 - freq));
                                invStd = (Std == 0) ? 0 : 1 / Std;
                                invstdvVec0_forVarRatio.push_back(invStd);
                                alleleFreqVec0_forVarRatio.push_back(freq);
                                MACVec0_forVarRatio.push_back(mac);
                                markerIndexVec0_forVarRatio.push_back(i);
                                SNPIdx_vr++;
                                numberofMarkers_varRatio++;
                        }
                }//end for(int i = 0; i < M; i++){
               

                if (minMAFtoConstructGRM > 0 || maxMissingRate < 1) {
                        cout << numberofMarkerswithMAFge_minMAFtoConstructGRM << " markers with MAF >= " << minMAFtoConstructGRM << " and missing rate <= " << maxMissingRate << endl;
                }

                int numofGenoArray_old = numofGenoArray;
                if (numberofMarkerswithMAFge_minMAFtoConstructGRM % numMarkersofEachArray == 0) {
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM / numMarkersofEachArray;
                } else {
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM / numMarkersofEachArray + 1;
                }

                if (numofGenoArray > numofGenoArray_old) {
                        for (int i = numofGenoArray; i < numofGenoArray_old; i++) {
                                delete genoVecofPointers[i];
                        }
                }

                invstdvVec.clear();
                invstdvVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
                alleleFreqVec.clear();
                alleleFreqVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
                MACVec.clear();
                MACVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);

                for (int i = 0; i < numberofMarkerswithMAFge_minMAFtoConstructGRM; i++) {
                        invstdvVec[i] = invstdvVec0.at(i);
                        alleleFreqVec[i] = alleleFreqVec0.at(i);
                        MACVec[i] = MACVec0.at(i);
                }

                if (isVarRatio) {
                        invstdvVec_forVarRatio.clear();
                        invstdvVec_forVarRatio.set_size(numberofMarkers_varRatio);
                        alleleFreqVec_forVarRatio.clear();
                        alleleFreqVec_forVarRatio.set_size(numberofMarkers_varRatio);
                        MACVec_forVarRatio.clear();
                        MACVec_forVarRatio.set_size(numberofMarkers_varRatio);
                        markerIndexVec_forVarRatio.clear();
                        markerIndexVec_forVarRatio.set_size(numberofMarkers_varRatio);

                        for (int i = 0; i < numberofMarkers_varRatio; i++) {
                                invstdvVec_forVarRatio[i] = invstdvVec0_forVarRatio.at(i);
                                alleleFreqVec_forVarRatio[i] = alleleFreqVec0_forVarRatio.at(i);
                                MACVec_forVarRatio[i] = MACVec0_forVarRatio.at(i);
                                markerIndexVec_forVarRatio[i] = markerIndexVec0_forVarRatio.at(i);
                        }
                }

                test_bedfile.close();
        }



        void NullGenoClass::printFromgenoVec(unsigned char genoBinary0) {
                unsigned char genoBinary = genoBinary0;
                for (int j = 0; j < 4; j++) {
                        int b = genoBinary & 1;
                        genoBinary = genoBinary >> 1;
                        int a = genoBinary & 1;
                        genoBinary = genoBinary >> 1;
                        cout << 2 - (a + b) << " " << endl;
                }
                cout << endl;
        }

        int NullGenoClass::getM() const {
                return M;
        }

        int NullGenoClass::getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const {
                return numberofMarkerswithMAFge_minMAFtoConstructGRM;
        }

        int NullGenoClass::getMsub() const {
                return Msub;
        }

        int NullGenoClass::getStartIndex() const {
                return startIndex;
        }

        int NullGenoClass::getEndIndex() const {
                return endIndex;
        }

        int NullGenoClass::getN() const {
                return N;
        }

        int NullGenoClass::getNnomissing() const {
                return Nnomissing;
        }

        double NullGenoClass::getAC(int m) {
                return alleleFreqVec[m] * 2 * Nnomissing;
        }

        double NullGenoClass::getMAC(int m) {
                if (alleleFreqVec[m] > 0.5) {
                        return (1 - alleleFreqVec[m]) * 2 * Nnomissing;
                } else {
                        return alleleFreqVec[m] * 2 * Nnomissing;
                }
        }

        int NullGenoClass::getMsub_MAFge_minMAFtoConstructGRM_in() const {
                return numberofMarkerswithMAFge_minMAFtoConstructGRM;
        }

        int NullGenoClass::getMsub_MAFge_minMAFtoConstructGRM_singleChr_in() const {
                return Msub_MAFge_minMAFtoConstructGRM_singleChr;
        }

        void NullGenoClass::Get_Samples_StdGeno(arma::ivec SampleIdsVec) {
                int indexOfVectorPointer;
                int SNPIdxinVec;

                int numSamples = SampleIdsVec.n_elem;
                stdGenoforSamples.clear();
                stdGenoforSamples.resize(M * numSamples);

                arma::ivec sampleGenoIdxVec(numSamples);
                arma::ivec sampleGenoIdxSubVec(numSamples);

                for (int j = 0; j < numSamples; j++) {
                        sampleGenoIdxVec[j] = SampleIdsVec[j] / 4;
                        sampleGenoIdxSubVec[j] = SampleIdsVec[j] % 4;
                }

                int startidx;
                unsigned char geno1;

                for (int i = 0; i < M; i++) {
                        indexOfVectorPointer = i / numMarkersofEachArray;
                        SNPIdxinVec = i % numMarkersofEachArray;
                        startidx = m_size_of_esi * SNPIdxinVec;

                        double freq = alleleFreqVec[i];
                        double invStd = invstdvVec[i];

                        for (int j = 0; j < numSamples; j++) {
                                int k = startidx + sampleGenoIdxVec[j];
                                geno1 = genoVecofPointers[indexOfVectorPointer]->at(k);
                                for (int q = 0; q < 4; q++) {
                                        if (q == sampleGenoIdxSubVec[j]) {
                                                int b = geno1 & 1;
                                                
                                                geno1 = geno1 >> 1;
                                                int a = geno1 & 1;
                                                stdGenoforSamples[i * numSamples + j] = ((2 - (a + b)) - 2 * freq) * invStd;
                                                geno1 = geno1 >> 1;
                                        } else {
                                                geno1 = geno1 >> 1;
                                                geno1 = geno1 >> 1;
                                        }
                                }
                        }
                }
        }

        void NullGenoClass::setstartendIndexVec(arma::ivec & t_startIndexVec, arma::ivec & t_endIndexVec) {
                startIndexVec = t_startIndexVec;
                endIndexVec = t_endIndexVec;
        }

        void NullGenoClass::setStartEndIndex(int startIndex, int endIndex, int chromIndex){
                startIndex = startIndex;
                endIndex = endIndex;
                Msub = 0;
                chromIndex = chromIndex;

                for(size_t i=0; i< M; i++){
                        if(i < startIndex || i > endIndex){
                                if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){

                                        Msub = Msub + 1;
                                }
                        }
                }
        }

} // end of namespace
