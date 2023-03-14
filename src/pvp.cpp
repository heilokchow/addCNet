
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;

// [[Rcpp::export]]
Rcpp::List pvp(Eigen::MatrixXd P) {
  MatrixXd Q = P.transpose();
  MatrixXd Va = Q * P;
  MatrixXd Pa = Va.inverse() * Q;

  Rcpp::List output = Rcpp::List::create(Rcpp::Named("Va") = Va,
                                         Rcpp::Named("Pa") = Pa);
  return(output);
}
