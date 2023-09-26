#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

#include "utility.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// only gives the 1st moment
// [[Rcpp::export(rng = false)]]
arma::vec eResC(arma::vec const Time,
		arma::vec const censor,
		arma::vec const wgt) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros);
  arma::uvec ind(Time.n_elem, arma::fill::zeros); 
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    ind(ind1).fill(i);
    d[i] = sum(censor.elem(ind1) % wgt.elem(ind1));
    r(span(0, i)) += sum(wgt.elem(ind1));
  }
  arma::vec surv = cumprod(1 - d / r);
  T0.resize(n + 1);
  T0(n) = T0(n - 1);
  arma::vec tmp = diff(T0) % surv;
  tmp.replace(datum::nan, 0);
  arma::vec eresT0 = (reverse(cumsum(reverse(tmp)))) / surv + T0(span(0, n - 1));
  eresT0.replace(datum::nan, max(T0));
  arma::vec eres(Time.n_elem, arma::fill::zeros);
  eres = eresT0(ind);
  return eres;
}

// Gives the 1st and 2nd moment
// [[Rcpp::export(rng = false)]]
Rcpp::List eResC2(arma::vec const Time,
		  arma::vec const censor,
		  arma::vec const wgt) {
  Rcpp::List out(2);
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros);
  arma::uvec ind(Time.n_elem, arma::fill::zeros); 
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    ind(ind1).fill(i);
    d[i] = sum(censor.elem(ind1) % wgt.elem(ind1));
    r(span(0, i)) += sum(wgt.elem(ind1));
  }
  arma::vec surv = cumprod(1 - d / r);
  T0.resize(n + 1);
  T0(n) = T0(n - 1);
  arma::vec tmp = diff(T0) % surv;
  tmp.replace(datum::nan, 0);
  arma::vec T00 = T0(span(0, n - 1));
  arma::vec eresT0 = (reverse(cumsum(reverse(tmp)))) / surv + T00;
  arma::vec eres2T0 = 2 * (reverse(cumsum(reverse(T00 % tmp)))) / surv + T00 % T00;
  eresT0.replace(datum::nan, max(T0));
  eres2T0.replace(datum::nan, max(T0 % T0));
  arma::vec eres(Time.n_elem, arma::fill::zeros);
  arma::vec eres2(Time.n_elem, arma::fill::zeros);
  eres = eresT0(ind);
  eres2 = eres2T0(ind);
  out(0) = eres;
  out(1) = eres2; 
  return out;
}

// Only can do indep, ex, and ar1. 
double get_alpha_delta(arma::vec mu,
		       arma::vec nt,
		       arma::vec index,
		       arma::vec w,
		       double a,
		       std::string corstr) {
  double phi = sum(w % square(mu)) / sum(w);
  // Rcpp::Rcout << "phi = " << phi << endl;		
  double H = 0;
  double G = 0; 
  int n = nt.n_elem;
  for (int i = 0; i < n; i++) {
    int s1 = nt(i);
    int crs = s1 * (s1 - 1) / 2;
    if (s1 > 1) {
      arma::vec PRi = mu(span(index(i), index(i) + s1 - 1));
      arma::vec Wi = w(span(index(i), index(i) + s1 - 1));
      arma::vec sPRi = PRi / sqrt(phi);
      arma::vec zi = genzi(sPRi);
      arma::mat Rhat(s1, s1, arma::fill::eye);
      arma::vec E(crs, arma::fill::ones);
      if (corstr == "exchangeable") {
	Rhat = Rhat * (1 - a) + a; 
      }
      if (corstr == "ar1") {
	arma::vec tmp(s1 - 1, arma::fill::value(a));
	arma::mat Rhat2(s1, s1, arma::fill::zeros);
	tmp = cumprod(tmp);
	int k = 0;
	for (int j = 0; j < s1 - 1; j++) {
	  Rhat2.submat(j + 1, j, s1 - 1, j) = tmp(span(0, s1 - j - 2));
	  E(span(k, k + s1 - 2 - j)) = regspace(1, s1 - 1 - j);
	  k += s1 - j - 1;
	}
	for (int j = 0; j < crs; j++) {
	  if (E(j) > 1) E(j) = E(j) * pow(a, E(j) - 1); 
	}
	Rhat += Rhat2 + Rhat2.t();
      }
      // Rcpp::Rcout << "E = " << E << endl;	
      // Rcpp::Rcout << "Rhat = " << Rhat << endl;
      // arma::vec rhoi = Rhat(arma::trimatu_ind(size(Rhat), 1));
      arma::vec rhoi = Rhat(arma::trimatl_ind(size(Rhat), -1));
      arma::vec WiV3inv = genzi(Wi);
      H += sum(E % E % WiV3inv);
      G += sum(E % WiV3inv % (zi - rhoi));
    } // end if s1 > 1
  } // end for i loop
  return G / H; // fabs(G / H);
}

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
  // Rcpp::Rcout << "index = " << index << endl;
  arma::vec b1 = b0;
  double ahat = 0;
  double ahatD = 0; 
  for (int j = 1; j <= maxit; j++) {
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    arma::vec mu = X * b0;
    arma::vec ym = y - mu;
    arma::mat bigD = matvec(X, sqrt(w));
    arma::vec ymw = ym % sqrt(w);
    if (max(nt) > 1) {
      ahatD = get_alpha_delta(ym, nt, index, w, ahat, corstr);
      ahat += ahatD;
    }
    // Rcpp::Rcout << "ahatD: " << ahatD << endl;    
    // Rcpp::Rcout << "ahat: " << ahat << endl;    
    for (int i = 0; i < N; i++) {
      int k = nt(i);
      arma::mat Rhat(k, k, arma::fill::eye);
      if (corstr == "exchangeable" && k > 1) {
	Rhat = Rhat * (1 - ahat) + ahat;
      }
      if (corstr == "ar1" && k > 1) {
	arma::vec tmp(k - 1, arma::fill::value(ahat));
	arma::mat Rhat2(k, k, arma::fill::zeros);
	tmp = cumprod(tmp);
	for (int i = 0; i < k - 1; i++) {
	  Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
	}
	Rhat = Rhat + Rhat2 + Rhat2.t();
      }
      // Rcpp::Rcout << "ahat: " << ahat << endl;
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::vec ym2w = ymw(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2w;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }  // end loop i 
    b1 = b0 + pinv(H) * S;
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    // if(max(abs(b1 - b0) / b1) < tol & ahatD < tol) break;
    if(max(abs(b1 - b0) / b1) < tol) break;
    b0 = b1;
  }
  arma::mat E(nx, nx, arma::fill::zeros);
  out(3) = E;
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}


// Replacing aftgee.est() when all margin = 1
// [[Rcpp::export(rng = false)]]
Rcpp::List est_No_Margin(arma::vec y,
			 arma::mat X,
			 arma::vec D,
			 arma::vec b0,
			 arma::vec nt,
			 arma::vec w,
			 std::string corstr,
			 double tol,
			 int maxit) {
  Rcpp::List out(5);
  arma::vec iter_gee(maxit, arma::fill::zeros);
  Rcpp::List histBeta(maxit);
  arma::vec b1 = b0;
  for (int j = 1; j <= maxit; j++) {
    arma::vec e = y - X * b0;
    arma::vec eres = eResC(e, D, w);
    arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
    Rcpp::List fit = gee(Ey, X, b0, nt, w, corstr, tol, maxit);
    arma::vec b1 = fit(0);
    if(max(abs(b1 - b0) / abs(b1)) < tol) break;
    b0 = b1;
    out(0) = b1;
    out(1) = j;
    iter_gee(j - 1) = fit(5);
    histBeta(j - 1) = b1; 
    out(2) =  fit(6);
  }
  out(3) = nonzeros(iter_gee);
  out(4) = histBeta;
  out.names() = Rcpp::CharacterVector::create("b", "iter", "alpha", "iter_gee", "hist_b");
  return out;
}

// for resampling
// [[Rcpp::export]]
arma::mat resampling_No_Margin(arma::vec y,
			       arma::mat X,
			       arma::vec D,
			       arma::vec b0,
			       arma::vec nt,
			       arma::vec w,
			       std::string corstr,
			       int B, 
			       double tol,
			       int maxit) {
  arma::mat out(b0.n_elem, B, arma::fill::zeros);
  for (int b = 0; b < B; b++) {
    arma::vec b1 = b0;
    arma::vec z = -log(randu(y.n_elem));
    for (int j = 1; j <= maxit; j++) {
      arma::vec e = y - X * b0;
      arma::vec eres = eResC(e, D, w);
      arma::vec Ey = D % y + (1 - D) % (eres + X * b0);
      Rcpp::List fit = gee(Ey, X, b0, nt, z % w, corstr, tol, maxit);
      arma::vec b1 = fit(0);
      if(max(abs(b1 - b0) / abs(b1)) < tol) break;
      b0 = b1;
    }
    out.col(b) = b1;
  }
  return out.t();
}
