#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec getSuv(arma::vec Time, arma::vec censor, arma::vec wgt) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1) % wgt.elem(ind1));
    r.elem(arma::regspace<arma::uvec>(0, 1, i)) += sum(wgt.elem(ind1));
  }
	return(cumprod(1 - d / r));
}
