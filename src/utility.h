#ifndef UTILITY_H
#define UTILITY_H

#include <RcppArmadillo.h>

arma::mat matvec(arma::mat x, arma::vec y); 
arma::mat matvec2(arma::mat x, arma::vec y); 
arma::vec genzi(arma::vec x);
bool iseye(const arma::mat& M);

#endif
