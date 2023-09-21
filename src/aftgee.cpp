#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

arma::mat matvec(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_cols; i++) {
    out.col(i) = x.col(i) % y;
  }
  return out;
}

// Functions to calculate parameter in variance-covariance matrix
// Those functions assume equal cluster size and weights at cluster level

// double ahatEx(arma::vec a,
// 	      arma::vec nt,
// 	      arma::vec index,
// 	      arma::vec w) {
//   vec aa = a % sqrt(w);
//   double phi = sum(square(aa)) / sum(w); // this makes scale.fix = FALSE
//   int n = nt.n_elem;
//   double out = 0;
//   for (int i = 0; i < n; i++) {
//     arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
//     if (nt(i) > 1) 
//       out += (sum(a2) * sum(a2) - sum(a2 % a2)) * w(index(i)) / (nt(i) - 1) / sum(w) / phi;
//   }
//   return(out);
// }

double ahatEx(arma::vec a,
	      arma::vec nt,
	      arma::vec index,
	      arma::vec w) {
  vec aa = a % sqrt(w);
  double phi = sum(square(aa)) / sum(w); // this makes scale.fix = FALSE
  int n = nt.n_elem;
  double out = 0;
  double tmp = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
    if (nt(i) > 1) 
      out += (sum(a2) * sum(a2) - sum(a2 % a2)) * w(index(i));
  }
  for (int i = 0; i < n; i++) {
    if (nt(i) > 1)
      tmp += w(index(i)) * nt(i) * (nt(i) - 1); 
  }
  return(out / tmp / phi);
}

// double ahatAR1(arma::vec a,
// 	       arma::vec nt,
// 	       arma::vec index,
// 	       arma::vec w) {
//   int n = nt.n_elem;
//   vec aa = a % sqrt(w);
//   double phi = sum(square(aa)) / sum(w);
//   double out = 0;
//   for (int i = 0; i < n; i++) {
//     arma::vec a2 = aa(span(index(i), index(i) + nt(i) - 1));
//     arma::vec w2 = w(span(index(i), index(i) + nt(i) - 1));
//     if (nt(i) > 1) 
//       out += sum(a2(span(0, nt(i) - 2)) % a2(span(1, nt(i) - 1))) / sum(w2) / (nt(i) - 1) / phi;
//   }
//   return(out)
// }

double ahatAR1(arma::vec a,
	       arma::vec nt,
	       arma::vec index,
	       arma::vec w) {
  int n = nt.n_elem;
  vec aa = a % sqrt(w);
  double phi = sum(square(aa)) / sum(w);
  double out = 0;
  double tmp = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = aa(span(index(i), index(i) + nt(i) - 1));
    arma::vec w2 = w(span(index(i), index(i) + nt(i) - 1));
    if (nt(i) > 1) 
      out += sum(a2(span(0, nt(i) - 2)) % a2(span(1, nt(i) - 1)));
  }
  for (int i = 0; i < n; i++) {
        if (nt(i) > 1) 
	  tmp += w(index(i)) * (nt(i) - 1); 
  }
  return(out / tmp / phi);
}

// no penality version
// [[Rcpp::export(rng = false)]]
Rcpp::List gee(arma::vec y,
							 arma::mat X,
							 arma::vec b0,
							 arma::vec nt,
							 arma::vec w,
							 std::string corstr,
							 double tol,
							 int maxit){
  Rcpp::List out(7);
  int N = nt.n_elem;
  int nx = X.n_cols;
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  double ahat = 0;
  for (int j = 1; j <= maxit; j++) {
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    arma::vec mu = X * b0;
    arma::vec ym = y - mu;
    arma::mat bigD = matvec(X, sqrt(w));
    arma::vec ymw = ym % sqrt(w);
    for (int i = 0; i < N; i++) {
      // ////////////////////////////////////////////////////////////////////
      int k = nt(i);
      arma::mat Rhat(k, k, arma::fill::eye);
      // Rcpp::Rcout << "corstr: " << corstr << endl;
      if (corstr == "exchangeable" && k > 1) {
				ahat = ahatEx(ym / stddev(mu), nt, index, w);
				Rhat = Rhat * (1 - ahat) + ahat;
      }
      if (corstr == "ar1" && k > 1) {
				ahat = ahatAR1(ym / stddev(mu), nt, index, w);
				arma::vec tmp(k - 1, arma::fill::value(ahat));
				arma::mat Rhat2(k, k, arma::fill::zeros);
				tmp = cumprod(tmp);
				for (int i = 0; i < k - 1; i++) {
					Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
				}
				Rhat = Rhat + Rhat2 + Rhat2.t();
      }
      // ////////////////////////////////////////////////////////////////////         
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::vec ym2w = ymw(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2w;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    b1 = b0 + pinv(H) * S;
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  arma::mat E(nx, nx, arma::fill::zeros);
  out(3) = E;
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}

arma::vec getID(arma::vec nt) {
  arma::vec out(sum(nt));
  int offset = 0;
  int n = nt.n_elem;
  for (int i = 0; i < n; i++) {
    out(span(offset, nt(i) - 1 + offset)).fill(i);
    offset += nt(i);
  }
  return(out);
}

// [[Rcpp::export(rng = false)]]
arma::vec eResC(arma::vec const Time,
		arma::vec const censor,
		arma::vec const wgt) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1) % wgt.elem(ind1));
    r(span(0, i)) += sum(wgt.elem(ind1));
  }
  arma::vec surv = cumprod(1 - d / r);
  T0.resize(n + 1);
  T0(n) = T0(n - 1);
  arma::vec tmp = diff(T0) % surv;
  tmp.replace(datum::nan, 0);
  arma::vec eres(Time.n_elem, arma::fill::zeros);
  arma::uvec ind = sort_index(Time);
  eres(ind) = (reverse(cumsum(reverse(tmp)))) / surv + T0(span(0, n - 1));
  eres.replace(datum::nan, max(T0));
  return eres;
}


// aftgee algorithm in c, no penalty
// [[Rcpp::export(rng = false)]]
Rcpp::List aftgeeEst(arma::vec y,
		     arma::mat X,
		     arma::vec D,
		     arma::vec b0,
		     arma::vec nt,
		     arma::vec w,
		     std::string corstr,
		     double tol,
		     int maxit) {
  Rcpp::List out(3);
  arma::vec iter_gee(maxit, arma::fill::zeros);
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::vec e = y - X * b0;
    arma::vec eres = eResC(e, D, w);
    arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
    Rcpp::List fit = gee(Ey, X, b0, nt, w, corstr, tol, maxit);
    arma::vec b1 = fit(0);
    if(max(abs(b1 - b0)) < tol) break;
    b0 = b1;
    out(0) = b1;
    out(1) = j;
    iter_gee(j - 1) = fit(5);
  }
  out(2) = nonzeros(iter_gee);
  out.names() = Rcpp::CharacterVector::create("b", "iter", "iter_gee");
  return out;
}

