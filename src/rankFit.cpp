#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::rowvec log_ns_est(const arma::vec& a,
			const arma::mat& X,
			const arma::vec& D,
			const arma::vec& Y,
			const arma::vec& W) {
  int n = Y.n_elem;
  int p = X.n_cols;
  arma::vec yexa = Y % exp(-X * a);
  arma::uvec ind = stable_sort_index(yexa, "descend");
  arma::vec ordD = D(ind);
  arma::vec ordW = W(ind);
  arma::mat xz = X % repmat(W, 1, p);
  xz = cumsum(xz.rows(ind), 0);
  // Rcpp::Rcout << ind;
  arma::mat c1 = X.rows(ind);
  arma::vec tmp = cumsum(W(ind));
  // arma::vec tmp = cumsum(ebaxZ(ind));
  // arma::vec tmp = arma::regspace(1, n);
  arma::mat r = c1 - xz / repmat(tmp, 1, p);
  r.replace(arma::datum::nan, 0);
  return sum(repmat(ordW % ordD, 1, p) % r, 0);
}
