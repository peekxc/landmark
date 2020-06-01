// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// landmark_maxmin
IntegerVector landmark_maxmin(const NumericMatrix& x, const int n, const int seed_index);
RcppExport SEXP _landmark_landmark_maxmin(SEXP xSEXP, SEXP nSEXP, SEXP seed_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type seed_index(seed_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(landmark_maxmin(x, n, seed_index));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_landmark_landmark_maxmin", (DL_FUNC) &_landmark_landmark_maxmin, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_landmark(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
