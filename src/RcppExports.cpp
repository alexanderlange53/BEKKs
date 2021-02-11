// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// elimination_mat
arma::mat elimination_mat(int n);
RcppExport SEXP _BEKKs_elimination_mat(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(elimination_mat(n));
    return rcpp_result_gen;
END_RCPP
}
// commutation_mat
arma::mat commutation_mat(int n);
RcppExport SEXP _BEKKs_commutation_mat(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(commutation_mat(n));
    return rcpp_result_gen;
END_RCPP
}
// duplication_mat
arma::mat duplication_mat(int n);
RcppExport SEXP _BEKKs_duplication_mat(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(duplication_mat(n));
    return rcpp_result_gen;
END_RCPP
}
// loglike_bekk
double loglike_bekk(arma::vec theta, arma::mat r);
RcppExport SEXP _BEKKs_loglike_bekk(SEXP thetaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike_bekk(theta, r));
    return rcpp_result_gen;
END_RCPP
}
// score_bekk
arma::mat score_bekk(arma::mat theta, arma::mat r);
RcppExport SEXP _BEKKs_score_bekk(SEXP thetaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(score_bekk(theta, r));
    return rcpp_result_gen;
END_RCPP
}
// SigmaLagCr
arma::mat SigmaLagCr(arma::mat y, int Tob, int q, int p, double ucvar);
RcppExport SEXP _BEKKs_SigmaLagCr(SEXP ySEXP, SEXP TobSEXP, SEXP qSEXP, SEXP pSEXP, SEXP ucvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type Tob(TobSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ucvar(ucvarSEXP);
    rcpp_result_gen = Rcpp::wrap(SigmaLagCr(y, Tob, q, p, ucvar));
    return rcpp_result_gen;
END_RCPP
}
// GarchVariance
arma::mat GarchVariance(int Tob, int q, int p, double ucvar, arma::mat theta, arma::mat Z);
RcppExport SEXP _BEKKs_GarchVariance(SEXP TobSEXP, SEXP qSEXP, SEXP pSEXP, SEXP ucvarSEXP, SEXP thetaSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Tob(TobSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ucvar(ucvarSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(GarchVariance(Tob, q, p, ucvar, theta, Z));
    return rcpp_result_gen;
END_RCPP
}
// ScoreGarch
arma::mat ScoreGarch(arma::mat epsilon2, arma::mat Z, int Tob, int q, int p, arma::mat theta, double ucvar);
RcppExport SEXP _BEKKs_ScoreGarch(SEXP epsilon2SEXP, SEXP ZSEXP, SEXP TobSEXP, SEXP qSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP ucvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type epsilon2(epsilon2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type Tob(TobSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type ucvar(ucvarSEXP);
    rcpp_result_gen = Rcpp::wrap(ScoreGarch(epsilon2, Z, Tob, q, p, theta, ucvar));
    return rcpp_result_gen;
END_RCPP
}
// LikelihoodGarch
double LikelihoodGarch(arma::mat Z, int Tob, int q, int p, arma::mat theta, arma::mat epsilon2, double ucvar);
RcppExport SEXP _BEKKs_LikelihoodGarch(SEXP ZSEXP, SEXP TobSEXP, SEXP qSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP epsilon2SEXP, SEXP ucvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type Tob(TobSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type epsilon2(epsilon2SEXP);
    Rcpp::traits::input_parameter< double >::type ucvar(ucvarSEXP);
    rcpp_result_gen = Rcpp::wrap(LikelihoodGarch(Z, Tob, q, p, theta, epsilon2, ucvar));
    return rcpp_result_gen;
END_RCPP
}
// BhhhGarch
arma::vec BhhhGarch(arma::mat r2, int q, int p, arma::mat theta, arma::mat epsilon2, arma::mat Z, int Tob, int max_iter, double crit, double ucvar);
RcppExport SEXP _BEKKs_BhhhGarch(SEXP r2SEXP, SEXP qSEXP, SEXP pSEXP, SEXP thetaSEXP, SEXP epsilon2SEXP, SEXP ZSEXP, SEXP TobSEXP, SEXP max_iterSEXP, SEXP critSEXP, SEXP ucvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type epsilon2(epsilon2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type Tob(TobSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type crit(critSEXP);
    Rcpp::traits::input_parameter< double >::type ucvar(ucvarSEXP);
    rcpp_result_gen = Rcpp::wrap(BhhhGarch(r2, q, p, theta, epsilon2, Z, Tob, max_iter, crit, ucvar));
    return rcpp_result_gen;
END_RCPP
}
// YLagCr
arma::mat YLagCr(arma::mat y, int p);
RcppExport SEXP _BEKKs_YLagCr(SEXP ySEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(YLagCr(y, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BEKKs_elimination_mat", (DL_FUNC) &_BEKKs_elimination_mat, 1},
    {"_BEKKs_commutation_mat", (DL_FUNC) &_BEKKs_commutation_mat, 1},
    {"_BEKKs_duplication_mat", (DL_FUNC) &_BEKKs_duplication_mat, 1},
    {"_BEKKs_loglike_bekk", (DL_FUNC) &_BEKKs_loglike_bekk, 2},
    {"_BEKKs_score_bekk", (DL_FUNC) &_BEKKs_score_bekk, 2},
    {"_BEKKs_SigmaLagCr", (DL_FUNC) &_BEKKs_SigmaLagCr, 5},
    {"_BEKKs_GarchVariance", (DL_FUNC) &_BEKKs_GarchVariance, 6},
    {"_BEKKs_ScoreGarch", (DL_FUNC) &_BEKKs_ScoreGarch, 7},
    {"_BEKKs_LikelihoodGarch", (DL_FUNC) &_BEKKs_LikelihoodGarch, 7},
    {"_BEKKs_BhhhGarch", (DL_FUNC) &_BEKKs_BhhhGarch, 10},
    {"_BEKKs_YLagCr", (DL_FUNC) &_BEKKs_YLagCr, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_BEKKs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
