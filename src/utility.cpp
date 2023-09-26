#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

arma::mat matvec(arma::mat x, arma::vec y) {
	arma::mat out(x.n_rows, x.n_cols);
	for (size_t i = 0; i < x.n_cols; i++) {
		out.col(i) = x.col(i) % y;
	}
	return out;
}

arma::mat matvec2(arma::mat x, arma::vec y) {
	arma::mat out(x.n_rows, x.n_cols);
	for (size_t i = 0; i < x.n_rows; i++) {
		out.row(i) = x.row(i) % y.t();
	}
	return out;
}

bool iseye(const arma::mat& M) {
	int n = M.n_rows;
	arma::mat A(n, n, arma::fill::eye);
	return(arma::approx_equal(A, M, "absdiff", 0.001));
}

arma::vec genzi(arma::vec x) {
  arma::mat out = x * x.t();
  return out(arma::trimatl_ind(size(out), -1));
}
