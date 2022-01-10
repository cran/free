#' Regularized Estimating Equations
#'
#' Unified regularized estimating equation solver. Currently include one solver with the l1 penalty only. More solvers and penalties are under development.
#' 
#' @docType package
#' @name free-package
NULL

#' Main solver of \code{free}
#' 
#' @param p The dimension of the dataset
#' @param lambda Lasso regularization coefficient
#' @param est_func R function, the estimating function specified by the user
#' @param par_init Optional, initial value for parameter update
#' @param alpha Tuning parameter
#' @param tau Tuning parameter
#' @param maxit Maximum iterations
#' @param tol_ee Convergence criterion based on the update of the estimating function
#' @param tol_par Convergence criterion based on the update of the parameter
#' @param verbose logical, print updates
#' @returns A list containing the regularized estimating equation estimates and the number of iterations it takes to converge.
#' @examples
#' # Standardize data
#' dat <- scale(mtcars)
#' x <- as.matrix(dat[, -1])
#' y <- as.vector(dat[, 1])
#' n <- nrow(x)
#' p <- ncol(x)
#' 
#' # Specify estimating function
#' ufunc <- function(b) {
#'   1/n * crossprod(x, (x %*% b - y) )
#' }
#' 
#' # Set hyperparameters
#' tau <- 0.6
#' alpha <- 0.5
#' 
#' # Set regularization coefficient
#' lambda1 <- 0
#' free_R1 <- free_lasso(p = p,
#'                       lambda = lambda1,
#'                       est_func = ufunc,
#'                       par_init = rep(0, p),
#'                       alpha = alpha,
#'                       tau = tau,
#'                       maxit = 10000L,
#'                       tol_ee = 1e-20,
#'                       tol_par = 1e-10,
#'                       verbose = FALSE)
#' free_R1$coefficients
#' 
#' # Compare with lm() - very close
#' lm(y~x-1)$coefficients
#' 
#' # Set regularization coefficient
#' lambda2 <- 0.7
#' free_R2 <- free_lasso(p = p,
#'                       lambda = lambda2,
#'                       est_func = ufunc,
#'                       par_init = rep(0, p),
#'                       alpha = alpha,
#'                       tau = tau,
#'                       maxit = 10000L,
#'                       tol_ee = 1e-20,
#'                       tol_par = 1e-10,
#'                       verbose = FALSE)
#' free_R2$coefficients
#' 
#' @export
#' 
free_lasso <- function(p, lambda, est_func,
                       par_init, alpha, tau,
                       maxit = 1000L,
                       tol_ee = 1e-6,
                       tol_par = 1e-6,
                       verbose = FALSE) {
  
  if (missing(par_init)) par_init <- rep(0, p)
  
  free_km <- REE_KM(beta = par_init,
                    p = p,
                    reg_p = p,
                    U = est_func,
                    tau = tau,
                    alpha = alpha,
                    penalty = 'lasso',
                    lambda1 = lambda,
                    maxit = maxit,
                    tol_U = tol_ee,
                    tol_beta = tol_par,
                    verbose = verbose)
  
  return(free_km)
}