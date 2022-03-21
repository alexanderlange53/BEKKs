

#ifndef IndicatorFunctions_H
#define IndicatorFunctions_H


arma::mat indicatorFunction(arma::mat r, arma::mat signs);


arma::mat expected_indicator_value(arma::mat r, arma::mat signs);


arma::mat elimination_mat(const int& n);


arma::mat commutation_mat(const int& n);


arma::mat duplication_mat(const int& n);


arma::mat inv_gen(const arma::mat& m);


#endif
