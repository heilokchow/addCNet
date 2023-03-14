// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SimSetC
List SimSetC(int n, double shift1, double shift2, NumericVector Zij);
RcppExport SEXP _addCNet_SimSetC(SEXP nSEXP, SEXP shift1SEXP, SEXP shift2SEXP, SEXP ZijSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shift1(shift1SEXP);
    Rcpp::traits::input_parameter< double >::type shift2(shift2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Zij(ZijSEXP);
    rcpp_result_gen = Rcpp::wrap(SimSetC(n, shift1, shift2, Zij));
    return rcpp_result_gen;
END_RCPP
}
// pvp
Rcpp::List pvp(Eigen::MatrixXd P);
RcppExport SEXP _addCNet_pvp(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(pvp(P));
    return rcpp_result_gen;
END_RCPP
}
// splinecalc
void splinecalc(NumericVector Patemp, NumericVector N_tall, NumericVector N_tallM, NumericVector b, int n, int nb, int z1, int z2, int z3, double h1, double t, NumericVector t_sep_t);
RcppExport SEXP _addCNet_splinecalc(SEXP PatempSEXP, SEXP N_tallSEXP, SEXP N_tallMSEXP, SEXP bSEXP, SEXP nSEXP, SEXP nbSEXP, SEXP z1SEXP, SEXP z2SEXP, SEXP z3SEXP, SEXP h1SEXP, SEXP tSEXP, SEXP t_sep_tSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Patemp(PatempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N_tall(N_tallSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N_tallM(N_tallMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< int >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< int >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< int >::type z3(z3SEXP);
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t_sep_t(t_sep_tSEXP);
    splinecalc(Patemp, N_tall, N_tallM, b, n, nb, z1, z2, z3, h1, t, t_sep_t);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_addCNet_SimSetC", (DL_FUNC) &_addCNet_SimSetC, 4},
    {"_addCNet_pvp", (DL_FUNC) &_addCNet_pvp, 1},
    {"_addCNet_splinecalc", (DL_FUNC) &_addCNet_splinecalc, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_addCNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
