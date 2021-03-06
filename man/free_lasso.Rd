% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{free_lasso}
\alias{free_lasso}
\title{Main solver of \code{free}}
\usage{
free_lasso(
  p,
  lambda,
  est_func,
  par_init,
  alpha,
  tau,
  maxit = 1000L,
  tol_ee = 1e-06,
  tol_par = 1e-06,
  verbose = FALSE
)
}
\arguments{
\item{p}{The dimension of the dataset}

\item{lambda}{Lasso regularization coefficient}

\item{est_func}{R function, the estimating function specified by the user}

\item{par_init}{Optional, initial value for parameter update}

\item{alpha}{Tuning parameter}

\item{tau}{Tuning parameter}

\item{maxit}{Maximum iterations}

\item{tol_ee}{Convergence criterion based on the update of the estimating function}

\item{tol_par}{Convergence criterion based on the update of the parameter}

\item{verbose}{logical, print updates}
}
\value{
A list containing the regularized estimating equation estimates and the number of iterations it takes to converge.
}
\description{
Main solver of \code{free}
}
\examples{
# Standardize data
dat <- scale(mtcars)
x <- as.matrix(dat[, -1])
y <- as.vector(dat[, 1])
n <- nrow(x)
p <- ncol(x)

# Specify estimating function
ufunc <- function(b) {
  1/n * crossprod(x, (x \%*\% b - y) )
}

# Set hyperparameters
tau <- 0.6
alpha <- 0.5

# Set regularization coefficient
lambda1 <- 0
free_R1 <- free_lasso(p = p,
                      lambda = lambda1,
                      est_func = ufunc,
                      par_init = rep(0, p),
                      alpha = alpha,
                      tau = tau,
                      maxit = 10000L,
                      tol_ee = 1e-20,
                      tol_par = 1e-10,
                      verbose = FALSE)
free_R1$coefficients

# Compare with lm() - very close
lm(y~x-1)$coefficients

# Set regularization coefficient
lambda2 <- 0.7
free_R2 <- free_lasso(p = p,
                      lambda = lambda2,
                      est_func = ufunc,
                      par_init = rep(0, p),
                      alpha = alpha,
                      tau = tau,
                      maxit = 10000L,
                      tol_ee = 1e-20,
                      tol_par = 1e-10,
                      verbose = FALSE)
free_R2$coefficients

}
