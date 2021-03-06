// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sample_start_info
CharacterVector sample_start_info(NumericMatrix m1, NumericMatrix m2, NumericMatrix m3, NumericMatrix m4, NumericVector w);
RcppExport SEXP sequencesAnonymizer_sample_start_info(SEXP m1SEXP, SEXP m2SEXP, SEXP m3SEXP, SEXP m4SEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m3(m3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m4(m4SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(sample_start_info(m1, m2, m3, m4, w));
    return __result;
END_RCPP
}
// sample_starting_state
NumericVector sample_starting_state(NumericMatrix m1, NumericMatrix m2, NumericMatrix m3, NumericMatrix m4, NumericVector w);
RcppExport SEXP sequencesAnonymizer_sample_starting_state(SEXP m1SEXP, SEXP m2SEXP, SEXP m3SEXP, SEXP m4SEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m3(m3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m4(m4SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(sample_starting_state(m1, m2, m3, m4, w));
    return __result;
END_RCPP
}
// sample_transition
NumericVector sample_transition(NumericMatrix m1, NumericMatrix m2, NumericMatrix m3, NumericMatrix m4, NumericVector w);
RcppExport SEXP sequencesAnonymizer_sample_transition(SEXP m1SEXP, SEXP m2SEXP, SEXP m3SEXP, SEXP m4SEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m3(m3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m4(m4SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(sample_transition(m1, m2, m3, m4, w));
    return __result;
END_RCPP
}
