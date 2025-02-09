#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <armadillo>
#include <stdexcept>
#include <Rcpp.h>


namespace glmm {
struct GLMMResult {
    arma::vec theta;
    arma::vec coefficients;
    arma::vec linear_predictors;
    arma::vec fitted_values;
    arma::vec Y;
    arma::vec residuals;
    arma::mat cov;
    bool converged;
    std::vector<std::string> sampleID;
    GLMResults obj_noK;
    arma::vec y;
    arma::mat X;
    std::string traitType;
    bool isCovariateOffset;
    arma::vec varWeights;
    bool LOCO;
    arma::vec chromosomeStartIndexVec;
    arma::vec chromosomeEndIndexVec;
};

GLMMResult glmmkin_ai_PCG_Rcpp(const string& bedFile, const string& bimFile, const string& famFile, const mat& Xorig, bool isCovariateOffset, const List& fit0, vec tau, vec fixtau, int maxiter, double tol, bool verbose, int nrun, double tolPCG, int maxiterPCG, const DataFrame& subPheno, const vec& indicatorGenoSamplesWithPheno, const List& obj_noK, const List& out_transform, vec tauInit, double memoryChunk, bool LOCO, const vec& chromosomeStartIndexVec, const vec& chromosomeEndIndexVec, double traceCVcutoff, bool isCovariateTransform, bool isDiagofKinSetAsOne, bool isLowMemLOCO, const mat& covarianceIdxMat, bool isStoreSigma, bool useSparseGRMtoFitNULL, bool useGRMtoFitNULL, bool isSparseGRMIdentity) {
    clock_t t_begin = clock();
    std::cout << "Start time: " << t_begin << std::endl;

    arma::vec subSampleInGeno = Rcpp::as<arma::vec>(subPheno["IndexGeno"]);
    if (subSampleInGeno.is_empty()) {
        subSampleInGeno = Rcpp::as<arma::vec>(subPheno["IndexPheno"]);
    }

    if (verbose) {
        std::cout << "Start reading genotype plink file here" << std::endl;
    }

    // Read genotype plink file
    if (!useSparseGRMtoFitNULL && useGRMtoFitNULL && !bedFile.empty()) {
        std::cout << "Reading genotype plink file" << std::endl;
        // Implement setgeno function to read genotype plink file
        // setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne);
    }

    setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)


    if (verbose) {
        std::cout << "Genotype reading is done" << std::endl;
    }

    if (LOCO) {
        // Implement getQCdMarkerIndex and updateChrStartEndIndexVec functions
        // arma::vec MsubIndVec = getQCdMarkerIndex();
        // arma::vec chrVec = Rcpp::as<arma::vec>(data.table::fread(bimFile, header = false)[, 1]);
        // chrVec = chrVec.elem(arma::find(MsubIndVec == true));
        // Rcpp::List updatechrList = updateChrStartEndIndexVec(chrVec);
        // LOCO = Rcpp::as<bool>(updatechrList["LOCO"]);
        // chromosomeStartIndexVec = Rcpp::as<arma::vec>(updatechrList["chromosomeStartIndexVec"]);
        // chromosomeEndIndexVec = Rcpp::as<arma::vec>(updatechrList["chromosomeEndIndexVec"]);
    }

    arma::vec y = Rcpp::as<arma::vec>(fit0["y"]);
    int n = y.n_elem;
    arma::mat X = Rcpp::as<arma::mat>(fit0["X"]);
    arma::vec offset = Rcpp::as<arma::vec>(fit0["offset"]);
    if (offset.is_empty()) {
        offset = arma::vec(n, arma::fill::zeros);
    }

    arma::vec var_weights = Rcpp::as<arma::vec>(fit0["weights"]);
    Rcpp::List family = fit0["family"];
    arma::vec eta = Rcpp::as<arma::vec>(fit0["linear.predictors"]);
    arma::vec mu = Rcpp::as<arma::vec>(fit0["fitted.values"]);
    arma::vec mu_eta = Rcpp::as<arma::vec>(family["mu.eta"](eta));
    arma::vec Y = eta - offset + (y - mu) / mu_eta;

    if (var_weights.is_empty()) {
        var_weights = arma::vec(mu_eta.n_elem, arma::fill::ones);
    }

    arma::vec sqrtW = mu_eta / sqrt(1 / var_weights % family["variance"](mu));
    arma::vec W = arma::square(sqrtW);

    arma::vec alpha0 = Rcpp::as<arma::vec>(fit0["coef"]);
    arma::vec eta0 = eta;

    std::cout << "tauInit" << std::endl;
    tau.zeros();
    if (family["family"] == "poisson" || family["family"] == "binomial") {
        tau[0] = 1;
        fixtau[0] = 1;
        tauInit[0] = 1;
        arma::uvec idxtau = arma::find(fixtau == 0);
        std::cout << "fixtau " << fixtau << std::endl;
        std::cout << "tauInit " << tauInit << std::endl;
        std::cout << "idxtau " << idxtau << std::endl;
        if (arma::sum(tauInit(idxtau)) == 0) {
            tau(idxtau).fill(0.1);
        } else {
            tau(idxtau) = tauInit(idxtau);
        }
    } else {
        arma::uvec idxtau = arma::find(fixtau == 0);
        if (arma::sum(tauInit(idxtau)) == 0) {
            tau[0] = 1;
            tau(idxtau) = arma::var(Y) / tau.n_elem;
            if (std::abs(arma::var(Y)) < 0.1) {
                throw std::runtime_error("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T");
            }
        } else {
            tau(arma::find(fixtau == 0)) = tauInit(arma::find(fixtau == 0));
        }
    }

    std::cout << "initial tau is " << tau << std::endl;

    if (!covarianceIdxMat.is_empty()) {
        arma::uvec idxtau2 = arma::intersect(covarianceIdxMat.col(0), arma::find(fixtau == 0));
        std::cout << "covarianceIdxMat" << std::endl;
        std::cout << covarianceIdxMat << std::endl;
        std::cout << "idxtau2" << std::endl;
        std::cout << idxtau2 << std::endl;
        if (!idxtau2.is_empty()) {
            tau(idxtau2).zeros();
        }
        arma::vec Kmatdiag = getMeanDiagofKmat(LOCO);
        std::cout << Kmatdiag << std::endl;
        tau(arma::span(1, tau.n_elem - 1)) /= Kmatdiag;
    }

    std::cout << "tau" << std::endl;
    std::cout << tau << std::endl;

    if (isSparseGRMIdentity) {
        tau[1] = 0;
    }

    Rcpp::List re_coef = Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose, maxiterPCG, tolPCG, maxiter, false, var_weights);
    Rcpp::List re = getAIScore_multiV_eMat(re_coef["Y"], X, re_coef["W"], tau, fixtau, re_coef["Sigma_iY"], re_coef["Sigma_iX"], re_coef["cov"], nrun, maxiterPCG, tolPCG, traceCVcutoff, false);
    arma::vec tau0 = tau;
    arma::vec tau0_q2 = tau(arma::find(fixtau == 0));
    std::cout << "tau0_q2 a" << std::endl;
    std::cout << tau0_q2 << std::endl;
    std::cout << "idxtau" << std::endl;
    std::cout << arma::find(fixtau == 0) << std::endl;

    arma::vec tau_q2 = arma::max(0, tau0_q2 + arma::square(tau0_q2) % (re["YPAPY"] - re["Trace"]) / n);
    tau(arma::find(fixtau == 0)) = tau_q2;

    if (!covarianceIdxMat.is_empty()) {
        tau(arma::find(fixtau == 0) % arma::find(covarianceIdxMat.col(0) == 1)).zeros();
    }
    std::cout << "re$YPAPY" << std::endl;
    std::cout << re["YPAPY"] << std::endl;
    std::cout << "re$Trace" << std::endl;
    std::cout << re["Trace"] << std::endl;

    if (verbose) {
        std::cout << "Variance component estimates:" << std::endl;
        std::cout << tau << std::endl;
    }

    int maxiter_in = maxiter;
    if (isSparseGRMIdentity) {
        tau[1] = 0;
        maxiter_in = 0;
        alpha0 = re_coef["alpha"];
        tau0 = tau;
        std::cout << "tau0_v1: " << tau0 << std::endl;
        eta0 = eta;
    }

    for (int i = 0; i < maxiter_in; ++i) {
        if (verbose) std::cout << "\nIteration " << i << " " << tau << ":\n";
        alpha0 = re_coef["alpha"];
        tau0 = tau;
        std::cout << "tau0_v1: " << tau0 << std::endl;
        eta0 = eta;

        clock_t t_begin_Get_Coef = clock();
        re_coef = glmm::Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose, maxiterPCG, tolPCG, maxiter, false, var_weights);
        clock_t t_end_Get_Coef = clock();
        std::cout << "t_end_Get_Coef - t_begin_Get_Coef\n";
        std::cout << (t_end_Get_Coef - t_begin_Get_Coef) / CLOCKS_PER_SEC << " seconds" << std::endl;

        Rcpp::List fit = glmm::fitglmmaiRPCG_multiV_eMat(re_coef["Y"], X, re_coef["W"], tau, fixtau, re_coef["Sigma_iY"], re_coef["Sigma_iX"], re_coef["cov"], nrun, maxiterPCG, tolPCG, tol, traceCVcutoff, false);

        clock_t t_end_fitglmmaiRPCG = clock();
        std::cout << "t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n";
        std::cout << (t_end_fitglmmaiRPCG - t_begin_Get_Coef) / CLOCKS_PER_SEC << " seconds" << std::endl;

        tau = Rcpp::as<arma::vec>(fit["tau"]);
        arma::mat cov = Rcpp::as<arma::mat>(re_coef["cov"]);
        alpha0 = Rcpp::as<arma::vec>(re_coef["alpha"]);
        eta = Rcpp::as<arma::vec>(re_coef["eta"]);
        Y = Rcpp::as<arma::vec>(re_coef["Y"]);
        mu = Rcpp::as<arma::vec>(re_coef["mu"]);

        mu_eta = Rcpp::as<arma::vec>(family["mu.eta"](eta));
        if (var_weights.is_empty()) {
            sqrtW = mu_eta / sqrt(Rcpp::as<arma::vec>(family["variance"](mu)));
        } else {
            sqrtW = mu_eta / sqrt(1 / var_weights % Rcpp::as<arma::vec>(family["variance"](mu)));
        }
        W = arma::square(sqrtW);

        std::cout << arma::abs(tau - tau0) / (arma::abs(tau) + arma::abs(tau0) + tol) << std::endl;
        std::cout << "tau: " << tau << std::endl;
        std::cout << "tau0: " << tau0 << std::endl;

        if (arma::sum(tau(arma::span(1, tau.n_elem - 1))) == 0) {
            break;
        } else {
            if (arma::max(arma::abs(tau - tau0) / (arma::abs(tau) + arma::abs(tau0) + tol)) < tol) break;
            if (arma::max(tau) > std::pow(tol, -2)) {
            Rcpp::warning("Large variance estimate observed in the iterations, model not converged...");
            i = maxiter;
            break;
            }
        }
        }

        if (verbose) std::cout << "\nFinal " << tau << ":\n";

        re_coef = glmm::Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose, maxiterPCG, tolPCG, maxiter, false, var_weights);
        arma::mat cov = Rcpp::as<arma::mat>(re_coef["cov"]);
        alpha0 = Rcpp::as<arma::vec>(re_coef["alpha"]);
        eta = Rcpp::as<arma::vec>(re_coef["eta"]);
        Y = Rcpp::as<arma::vec>(re_coef["Y"]);
        mu = Rcpp::as<arma::vec>(re_coef["mu"]);
        bool converged = (i < maxiter);

        arma::vec res = y - mu;

        std::string traitType;
        if (family["family"] == "binomial") {
        mu2 = mu % (1 - mu);
        traitType = "binary";
        } else if (family["family"] == "poisson") {
        mu2 = mu;
        traitType = "count";
        } else if (family["family"] == "gaussian") {
        mu2 = arma::vec(res.n_elem, arma::fill::ones) / tau[0];
        traitType = "quantitative";
        }

        arma::vec mu2_rescaled = mu2 % var_weights;
        arma::vec y_rescaled = y % var_weights;
        arma::vec mu_rescaled = mu % var_weights;

    Rcpp::List obj_noK;
    if (!isCovariateOffset) {
        obj_noK = ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X);
    } else {
        obj_noK = ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig);
    }

    glmm::GLMMResult glmmResult;
    glmmResult.theta = tau;
    glmmResult.coefficients = alpha0;
    glmmResult.linear_predictors = eta;
    glmmResult.fitted_values = mu;
    glmmResult.Y = Y;
    glmmResult.residuals = res;
    glmmResult.cov = cov;
    glmmResult.converged = converged;
    glmmResult.sampleID = Rcpp::as<std::vector<std::string>>(subPheno["IID"]);
    glmmResult.obj_noK = obj_noK;
    glmmResult.y = y;
    glmmResult.X = X;
    glmmResult.traitType = traitType;
    glmmResult.isCovariateOffset = isCovariateOffset;
    glmmResult.varWeights = var_weights;
    glmmResult.LOCO = LOCO;

    clock_t t_end_null = clock();
    std::cout << "t_end_null - t_begin, fitting the NULL model without LOCO took\n";
    std::cout << (t_end_null - t_begin) / CLOCKS_PER_SEC << " seconds" << std::endl;

    if (!isLowMemLOCO && LOCO) {
        set_Diagof_StdGeno_LOCO();
        glmmResult.LOCOResult = Rcpp::List::create();
        for (int j = 0; j < 22; ++j) {
            int startIndex = chromosomeStartIndexVec[j];
            int endIndex = chromosomeEndIndexVec[j];
            if (!std::isnan(startIndex) && !std::isnan(endIndex)) {
                std::cout << "leave chromosome " << j + 1 << " out\n";
                setStartEndIndex(startIndex, endIndex, j);
                clock_t t_begin_Get_Coef_LOCO = clock();
                Rcpp::List re_coef_LOCO = glmm::Get_Coef_multiV(y, X, tau, family, alpha0, eta0, offset, verbose, maxiterPCG, tolPCG, maxiter, true, var_weights);
                clock_t t_end_Get_Coef_LOCO = clock();
                std::cout << "t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n";
                std::cout << (t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO) / CLOCKS_PER_SEC << " seconds" << std::endl;
                arma::mat cov = Rcpp::as<arma::mat>(re_coef_LOCO["cov"]);
                arma::vec alpha = Rcpp::as<arma::vec>(re_coef_LOCO["alpha"]);
                arma::vec eta = Rcpp::as<arma::vec>(re_coef_LOCO["eta"]);
                arma::vec Y = Rcpp::as<arma::vec>(re_coef_LOCO["Y"]);
                arma::vec mu = Rcpp::as<arma::vec>(re_coef_LOCO["mu"]);
                arma::vec res = y - mu;

                if (family["family"] == "binomial") {
                    mu2 = mu % (1 - mu);
                } else if (family["family"] == "poisson") {
                    mu2 = mu;
                } else if (family["family"] == "gaussian") {
                    mu2 = arma::vec(res.n_elem, arma::fill::ones) / tau[0];
                }

                arma::vec mu2_rescaled = mu2 % var_weights;
                arma::vec mu_rescaled = mu % var_weights;

                Rcpp::List obj_noK;
                if (!isCovariateOffset) {
                    obj_noK = ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, X);
                } else {
                    obj_noK = ScoreTest_NULL_Model(mu_rescaled, mu2_rescaled, y_rescaled, Xorig);
                }

                glmmResult.LOCOResult[j] = Rcpp::List::create(
                    Rcpp::Named("isLOCO") = true,
                    Rcpp::Named("coefficients") = alpha,
                    Rcpp::Named("linear.predictors") = eta,
                    Rcpp::Named("fitted.values") = mu,
                    Rcpp::Named("Y") = Y,
                    Rcpp::Named("residuals") = res,
                    Rcpp::Named("cov") = cov,
                    Rcpp::Named("obj.noK") = obj_noK
                );
            } else {
                glmmResult.LOCOResult[j] = Rcpp::List::create(Rcpp::Named("isLOCO") = false);
            }
        }
    }

    if (isLowMemLOCO && LOCO) {
        glmmResult.chromosomeStartIndexVec = chromosomeStartIndexVec;
        glmmResult.chromosomeEndIndexVec = chromosomeEndIndexVec;
    }

    return glmmResult;
}