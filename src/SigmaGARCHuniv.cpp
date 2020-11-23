#include <RcppArmadillo.h>
using namespace Rcpp;


// Univariate GARCH(1, 1) variances -------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec SigmaGARCHuniv(const arma::vec& param, int Tob, double& SigmaE, const arma::vec est){

  arma::vec SigmaOut(Tob);
  SigmaOut(0) = SigmaE;

  for (int i = 1; i < Tob; i++) {
    SigmaOut(i) =  param(0) + param(1) * pow(est(i-1), 2) + param(2) * SigmaOut(i-1);
  }

  return SigmaOut;
}
