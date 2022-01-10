n <- 500
p <- 10

beta_true <- rep(0, p)

beta_nz <- c(-3, 1, 2)
p_nz <- length(beta_nz)

beta_true[1:p_nz] <- beta_nz
n_unreg <- 0

# cov_mat <- matrix(NA, nrow = p, ncol = p)
# 
# for (i in 1:p) {
#   for (j in i:p) {
#     cov_mat[i, j] <- cov_mat[j, i] <- 0.3^abs(i - j)
#   }
# }
# 
# cov_chol <- chol(cov_mat)

x <- matrix(data = rnorm(n = n * p), nrow = n, ncol = p)

y <- x %*% beta_true + rnorm(n)

ufunc <- function(b) {
  1/n * crossprod(x, (x %*% b - y) )
}

tau <- 0.5
alpha <- 0.5
lambda1 <- 0.5

free_R <- free_lasso(p = p,
                     lambda = lambda1,
                     est_func = ufunc,
                     alpha = alpha,
                     tau = tau,
                     tol_ee = 1e-20,
                     tol_par = 1e-10,
                     verbose = TRUE)

free_R
