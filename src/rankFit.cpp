#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using cmp_par = std::pair<double, arma::uword>;

arma::mat matvec2(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_rows; i++) {
    out.row(i) = x.row(i) % y.t();
  }
  return out;
}

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
    out = D(idx_i) * W(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
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
