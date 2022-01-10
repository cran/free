// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// extern "C" {
// #include "aa.h"
// #include "aa_blas.h"
// }

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// Lasso proximal operator
// [[Rcpp::export]]
Rcpp::NumericVector prox_lasso(Rcpp::NumericVector &u,
                               int &dim_u,
                               int &unreg_p,    // the first 'unreg_p' variables are not regularized
                               double &tau_lam) {
  int j;
  
  for (j = unreg_p; j < dim_u; j++) {
    if (abs(u[j]) <= tau_lam)
      u[j] = 0.0;
    else if (u[j] > 0)
      u[j] = u[j] - tau_lam;
    else
      u[j] = u[j] + tau_lam;
  }
  return u;
}


// KM solver
// [[Rcpp::export]]
Rcpp::List REE_KM(Rcpp::NumericVector &beta,
                  int p,
                  int reg_p,
                  Rcpp::Function &U,
                  double tau,
                  double alpha,
                  std::string penalty,
                  double lambda1,
                  int maxit,
                  double tol_U,
                  double tol_beta,
                  bool verbose = false) {
  Rcpp::NumericVector beta_prev (p);
  int i;
  double s_U;
  double s_beta;
  double lam1;
  int unreg_p;
  Rcpp::NumericVector ubeta (p, 100.0); // arbitrary initial value to past the first convergence check
  Rcpp::NumericVector ubeta_prev (p);
  Rcpp::NumericVector s_beta_output (maxit);
  
  // First 'unreg_p' variables are not penalized
  unreg_p = p - reg_p;
  
  i = 0;
  while (i < maxit) {
    // store current beta in beta_prev
    beta_prev = clone(beta);
    
    // store previous ubeta in ubeta_prev
    ubeta_prev = clone(ubeta);
    
    // compute intermediates steps: U(beta)
    ubeta = U(beta);
    
    // || U^k ||_2 - || U^k-1 ||_2 convergence check
    s_U = sqrt(sum(pow(ubeta - ubeta_prev, 2)));
    s_U = s_U/p;
    if (s_U < tol_U) break;
    // end of check. 
    
    // continue to computing: beta - tau * U(beta)
    beta = beta_prev - tau * ubeta;
    
    // compute: prox_{tau lambdas}(beta - tau * U(beta) )
    if (penalty == "lasso") {
      lam1 = tau * lambda1;
      prox_lasso(beta, p, unreg_p, lam1);
    }
    
    beta = (1 - alpha) * beta_prev + alpha * beta;
    
    // convergence check
    s_beta = sqrt(sum(pow(beta - beta_prev, 2)));
    s_beta = s_beta/p;
    s_beta_output[i] = s_beta;
    
    i++;
    
    if (verbose) Rcpp::Rcout << "Mean L2 norm of coefficient update \n" << s_beta <<"\n";
    if (s_beta < tol_beta) break;
  }
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                            Rcpp::Named("iterations") = i);
                            // Rcpp::Named("Residuals") = s_beta_output);
}
