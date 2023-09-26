#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace Rcpp;

#include "utility.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using cmp_par = std::pair<double, arma::uword>;

// log-rank type estimating function (non-smooth); old name = ulognsfun */
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec log_ns_est(const arma::vec& a,
		     const arma::mat& X,
		     const arma::vec& D,
		     const arma::vec& Y,
		     const arma::vec& W,
		     const arma::vec& gw) {
  arma::uword const n = Y.n_elem;
  arma::uword const p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec yexa = Y - X.t() * a;
  yexa.replace(-arma::datum::inf, arma::datum::nan);
  yexa.replace(arma::datum::nan, yexa.min());
  arma::uvec const idx = arma::sort_index(yexa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  arma::vec x_col_sum(p, arma::fill::zeros);
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    x_col_sum = sum(matvec2(X, W), 1);
    w_sum = sum(W);
    out = D(idx_i) * W(idx_i) * gw(idx_i) * (X.col(idx_i) - x_col_sum / w_sum);
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first) {
      while(yexa[idx_i] > indices_head->first) {
	x_col_sum -= W(indices_head->second) * X.col(indices_head->second);
	w_sum -= W(indices_head->second);
	++indices_head;
      }
    }
    else --indices_head;
    if (w_sum > 0)
      out += D(idx_i) * W(idx_i) * gw(idx_i) * (X.col(idx_i) - x_col_sum / w_sum);
  }
  return out;
}

// Gehan type estimating function (smooth); old name = ufun */
// @noRd
// [[Rcpp::export(rng = false)]]
arma::vec gehan_s_est(const arma::vec& a,
		      const arma::mat& X,
		      const arma::vec& D,
		      const arma::vec& Y,
		      const arma::vec& W,
		      const int& nc,
		      const arma::mat& sigma,
		      const arma::vec& gw) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  arma::vec yexa = Y - X * a;
  arma::vec Dgw = D % gw;
  arma::mat cSigma(p, p, arma::fill::eye);
  if (iseye(sigma) == false) cSigma = chol(sigma).t();
  for (int i = 0; i < n - 1; i++) {
    arma::mat xdif = repmat(X.row(i), n - i - 1, 1) - X.rows(i + 1, n - 1);
    arma::mat xs = xdif;
    if (iseye(sigma) == false) xs = xdif * cSigma;
    arma::vec rij = sqrt(sum(xs % xs, 1));
    arma::vec yexa2 = yexa.subvec(i + 1, n - 1);
    arma::vec H = arma::normcdf(sqrt(nc) * (yexa2 - yexa[i]) / rij); 
    H.replace(arma::datum::nan, 0);
    arma::vec Dgwj = Dgw.subvec(i + 1, n - 1);
    arma::vec Wj = W.subvec(i + 1, n - 1);
    out += W(i) * sum(matvec(xdif, (Dgw(i) + Dgwj) % Wj % H - Dgwj % Wj)).t();
  }
  return out;
}


// Faster version if no ties
// 
// // [[Rcpp::export(rng = false)]]
// arma::rowvec log_ns_est(const arma::vec& a,
// 			const arma::mat& X,
// 			const arma::vec& D,
// 			const arma::vec& Y,
// 			const arma::vec& W,
// 			const arma::vec& gw) {
//   int n = Y.n_elem;
//   int p = X.n_cols;
//   arma::vec yexa = Y % exp(-X * a);
//   arma::uvec ind = stable_sort_index(yexa, "descend");
//   arma::vec ordD = D(ind);
//   arma::vec ordW = W(ind);
//   arma::vec ordgw = gw(ind);
//   arma::mat xz = X % repmat(W, 1, p);
//   xz = cumsum(xz.rows(ind), 0);
//   arma::mat c1 = X.rows(ind);
//   arma::vec tmp = cumsum(W(ind));
//   arma::mat r = c1 - xz / repmat(tmp, 1, p);
//   r.replace(arma::datum::nan, 0);
//   return sum(repmat(ordgw % ordW % ordD, 1, p) % r, 0);
// }

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec gehan_ns_est(const arma::vec& a,
		       const arma::mat& X,
		       const arma::vec& D,
		       const arma::vec& Y,
		       const arma::vec& W,
		       const arma::vec& gw) {
  arma::uword const n = Y.n_elem;
  arma::uword const p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec yexa = Y - X.t() * a;
  yexa.replace(-arma::datum::inf, arma::datum::nan);
  yexa.replace(arma::datum::nan, yexa.min() - 0.01);
  arma::uvec const idx = arma::sort_index(yexa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  arma::vec x_col_sum(p, arma::fill::zeros);
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    x_col_sum = sum(matvec2(X, W), 1);
    w_sum = sum(W);
    out = D(idx_i) * W(idx_i) * gw(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first) {
      while(yexa[idx_i] > indices_head->first) {
	x_col_sum -= W(indices_head->second) * X.col(indices_head->second);
	w_sum -= W(indices_head->second);
	++indices_head;
      }
    }
    else --indices_head;
    if (w_sum > 0)
      out += D(idx_i) * W(idx_i) * gw(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
  }
  return out;
}

// Compute smooth gehan weight; used to prepare for method #3 and #4; old name = getgehan
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec gehan_s_wt(const arma::vec& a,
		     const arma::mat& X,
		     const arma::vec& Y,
		     const arma::vec& W,
		     const int& nc,
		     const arma::mat& sigma) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec out(n, arma::fill::zeros);
  arma::vec yexa = Y - X * a;
  arma::mat cSigma(p, p, arma::fill::eye);
  if (iseye(sigma) == false) cSigma = chol(sigma).t();
  for (int i = 0; i < n; i++) {
    arma::mat xdif = repmat(X.row(i), n, 1) - X;
    arma::mat xs = xdif;
    if (iseye(sigma) == false) xs = xdif * cSigma;
    arma::vec rij = sqrt(sum(xs % xs, 1));
    arma::vec H = arma::normcdf(sqrt(nc) * (yexa - yexa[i]) / rij);
    // std::cout << H << "\n";
    H.replace(arma::datum::nan, 0);
    out[i] = sum(H % W);
  }
  return out;
}

// Compute non-smooth gehan weight; used to prepare for method #3 and #4; old name = getnsgehan
// [[Rcpp::export(rng = false)]]
arma::vec gehan_ns_wt(const arma::vec& a,
		      const arma::mat& X,
		      const arma::vec& Y,
		      const arma::vec& W) {
  arma::uword const n = Y.n_elem;
  arma::vec out(n, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec yexa = Y - X.t() * a;
  yexa.replace(-arma::datum::inf, arma::datum::nan);
  yexa.replace(arma::datum::nan, yexa.min() - 0.01);
  arma::uvec const idx = arma::sort_index(yexa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    w_sum = sum(W);
    out[idx_i] = w_sum;
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first) {
      while(yexa[idx_i] > indices_head->first) {
	w_sum -= W(indices_head->second);
	++indices_head;
      }
    }
    else --indices_head;
    out[idx_i] = w_sum;
  }
  return out;
}

// log-rank type estimating function (smooth-equivent form); old name = ulogfun
// @noRd
// [[Rcpp::export(rng = false)]]
arma::vec log_s_est(const arma::vec& a,
		    const arma::mat& X,
		    const arma::vec& D,
		    const arma::vec& Y,
		    const arma::vec& W,
		    const int& nc,
		    const arma::mat& sigma,
		    const arma::vec& gw) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::rowvec out(p, arma::fill::zeros);
  arma::vec yexa = Y - X * a; 
  arma::mat cSigma(p, p, arma::fill::eye);
  if (iseye(sigma) == false) cSigma = chol(sigma).t();
  for (int i = 0; i < n; i++) {
    if (D(i) > 0) {
      arma::mat xdif = repmat(X.row(i), n, 1) - X;
      arma::mat xs = xdif;
      if (iseye(sigma) == false) xs = xdif * cSigma;
      arma::vec rij = sqrt(sum(xs % xs, 1));
      arma::vec H = arma::normcdf(sqrt(nc) * (yexa - yexa[i]) / rij);
      H.replace(arma::datum::nan, 0);
      if (sum(H) != 0) 
	out += gw(i) * W(i) * (X.row(i) - sum(matvec(X, H)) / sum(H));
    }
  }
  return out.t();
}
