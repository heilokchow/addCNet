#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void cumcalc(NumericVector Patemp, NumericVector N_tC, NumericVector B, int n,
                int nb, int z3, double t, NumericVector t_sep_t) {
  int z33 = z3 - 1;
  double ktemp(0.0);
  int n1 = n*(n-1);

  for (int j = 0; j < 100; j++) {
    if (t < t_sep_t[j]) {
      for (int i = 0; i < nb; i++) {
        B[j*nb + i] += Patemp[i];
      }
    }
    N_tC[j*n1 + z33] += 1;
  }
}

