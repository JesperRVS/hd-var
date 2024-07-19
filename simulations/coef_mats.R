## Dependencies
library("Matrix") # for triu and Diagonal

# Coefficient matrix for the VAR(1) process
# INPUTS:
#   theta:      VAR(1) coefficient, scalar
#   p:          System dimension, integer
# OUTPUT:
#   theta:      p x p VAR(1) coefficient matrix (diagonal, sparse)
coef_mat_a <- function(p) {
  theta <- 0.5
  m <- Matrix::Diagonal(n = p, x = theta)
  return(m)
}

# Coefficient matrix for the VAR(4) process
# INPUT:
#   p:          System dimension, divisible by four
# OUTPUT:
#   theta:      p x 4p coefficient matrix
coef_mat_b <- function(p) {
  if (p %% 4 != 0) {
    warning("p must be divisible by four.")
  }
  pover4 <- p / 4                           # p/4 (num blocks)
  block1 <- matrix(0.15, 4, 4)              # block in theta1
  ipover4 <- Matrix::Diagonal(pover4)       # sparse identity I_{p/4}
  theta1 <- kronecker(ipover4, block1)      # p x p block diagonal
  block4 <- matrix(-0.1, 4, 4)              # block in theta4
  theta4 <- kronecker(ipover4, block4)      # p x p block diagonal
  theta23 <- matrix(0, p, 2 * p)            # zeros in theta2 and theta3
  theta <- cbind(theta1, theta23, theta4)   # p x 4p coefficient matrix
  return(theta)
}

# Coefficient matrix for the VAR(1) process
# INPUT:
#   p:          System dimension, positive integer
# OUTPUT:
#   theta:      p x p coefficient matrix
coef_mat_c <- function(p) {
  a <- 0.4                          # parameter (hardcoded)
  theta_temp <- matrix(NA, p, p)    # p x p temporary matrix
  for (i in 1:(p - 1)) {            # fill in strict upper triangle
    for (j in 2:p) {
      theta_temp[i, j] <- ((-1)^(abs(i - j))) * (a^(abs(i - j) + 1))
    }
  }
  tri_up <- Matrix::triu(theta_temp, 1) # strict upper triangle
  di <- Matrix::Diagonal(n = p, x = a)  # diagonal (sparse)
  theta <- di + tri_up + t(tri_up)      # symmetrize
  return(theta)
}