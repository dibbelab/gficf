// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RunModularityClusteringCpp
IntegerVector RunModularityClusteringCpp(Eigen::SparseMatrix<double> SNN, int modularityFunction, double resolution, int algorithm, int nRandomStarts, int nIterations, int randomSeed, bool printOutput, std::string edgefilename);
RcppExport SEXP _gficf_RunModularityClusteringCpp(SEXP SNNSEXP, SEXP modularityFunctionSEXP, SEXP resolutionSEXP, SEXP algorithmSEXP, SEXP nRandomStartsSEXP, SEXP nIterationsSEXP, SEXP randomSeedSEXP, SEXP printOutputSEXP, SEXP edgefilenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type SNN(SNNSEXP);
    Rcpp::traits::input_parameter< int >::type modularityFunction(modularityFunctionSEXP);
    Rcpp::traits::input_parameter< double >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< int >::type algorithm(algorithmSEXP);
    Rcpp::traits::input_parameter< int >::type nRandomStarts(nRandomStartsSEXP);
    Rcpp::traits::input_parameter< int >::type nIterations(nIterationsSEXP);
    Rcpp::traits::input_parameter< int >::type randomSeed(randomSeedSEXP);
    Rcpp::traits::input_parameter< bool >::type printOutput(printOutputSEXP);
    Rcpp::traits::input_parameter< std::string >::type edgefilename(edgefilenameSEXP);
    rcpp_result_gen = Rcpp::wrap(RunModularityClusteringCpp(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename));
    return rcpp_result_gen;
END_RCPP
}
// jaccard_coeff
NumericMatrix jaccard_coeff(NumericMatrix idx, bool printOutput);
RcppExport SEXP _gficf_jaccard_coeff(SEXP idxSEXP, SEXP printOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< bool >::type printOutput(printOutputSEXP);
    rcpp_result_gen = Rcpp::wrap(jaccard_coeff(idx, printOutput));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_WMU_test
Rcpp::NumericMatrix rcpp_WMU_test(Rcpp::NumericMatrix M, Rcpp::NumericVector idx1, Rcpp::NumericVector idx2);
RcppExport SEXP _gficf_rcpp_WMU_test(SEXP MSEXP, SEXP idx1SEXP, SEXP idx2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type idx1(idx1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type idx2(idx2SEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_WMU_test(M, idx1, idx2));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_parallel_jaccard_coef
Rcpp::NumericMatrix rcpp_parallel_jaccard_coef(Rcpp::NumericMatrix mat, bool printOutput);
RcppExport SEXP _gficf_rcpp_parallel_jaccard_coef(SEXP matSEXP, SEXP printOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< bool >::type printOutput(printOutputSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_parallel_jaccard_coef(mat, printOutput));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_parallel_WMU_test
Rcpp::NumericMatrix rcpp_parallel_WMU_test(Rcpp::NumericMatrix matX, Rcpp::NumericMatrix matY, bool printOutput);
RcppExport SEXP _gficf_rcpp_parallel_WMU_test(SEXP matXSEXP, SEXP matYSEXP, SEXP printOutputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type matX(matXSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type matY(matYSEXP);
    Rcpp::traits::input_parameter< bool >::type printOutput(printOutputSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_parallel_WMU_test(matX, matY, printOutput));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gficf_RunModularityClusteringCpp", (DL_FUNC) &_gficf_RunModularityClusteringCpp, 9},
    {"_gficf_jaccard_coeff", (DL_FUNC) &_gficf_jaccard_coeff, 2},
    {"_gficf_rcpp_WMU_test", (DL_FUNC) &_gficf_rcpp_WMU_test, 3},
    {"_gficf_rcpp_parallel_jaccard_coef", (DL_FUNC) &_gficf_rcpp_parallel_jaccard_coef, 2},
    {"_gficf_rcpp_parallel_WMU_test", (DL_FUNC) &_gficf_rcpp_parallel_WMU_test, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gficf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
