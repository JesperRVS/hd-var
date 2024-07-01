lasso_var <- function(data, q = 1, post = TRUE, intr = TRUE,
                      c = 1.1, gamma = 0.1 / log(max(dim(data))), k = 15,
                      tol_ups = 1e-3, tol_glmnet = 1e-4, warn = TRUE) {
  if (missing(data) || is.null(data)) {
    stop("not enough input arguments")
  }
  # Unpack data
  nplusq <- nrow(data)  # n = effective sample size, q = autoregressive order
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q
  qp <- q * p           # qp = total number of variables
  y <- data[(q + 1):(q + n), ] # response
  x <- matrix(NA, nrow = n, ncol = qp) # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- (ell - 1) * p + 1:ell * p
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,...
    x[, block_ell] <- data[lag_ell, ]
  }
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  #  If intercepts requested, demean before proceeding
  if (intr == TRUE) {
    ybar <- colMeans(y)
    y <- sweep(y, 2, ybar)
    xbar <- colMeans(x)
    x <- sweep(x, 2, xbar)
  }
  ups_init <- sqrt(t(y^2) %*% (x^2))
  that_init <- mlasso(x, y, lambda_glmnet, ups_init, tol_glmnet)
  if (post == TRUE) {
    sel <- != x
    for (i in 1:p) {
      that_init[i,ix ] <- solve(t(X[, sel[i, ]]) %*% X[, sel[i, ]]) %*% t(X[, sel[i, ]]) %*% Y[, i]
    }
  }
  
  if (intr) {
    constInit <- - thatInix %y% Xbar
  } else {
    constInit <- NULL
  }
  
  UpsRefi <- array(NA, dim = c(p, qp, K))
  ThatRefi <- array(NA, dim = c(p, qp, K))
  if (intr) {
    constRefi <- matrix(NA, nrow = p, ncol = K)
  } else {
    constRefi <- NULL
  }
  
  for (k in 1:K) {
    if (k == 1) {
      Res_kminus1 <- Y - X %*% t(that_init)
    } else       Res_kminus1 <x Yy- X %*% t(ThatRefi[, , k - 1])
    }
    
    UpsRefi[, , k] <- sqrt(colMeans(Res_kminus1^2) %*% t(colMeans(X^2)))
    ThatRefi_k <- Design_A$mlasso(X, Y, lambda_glmnet, UpsRefi[, , k], tol_glmnet)
    
    if (post) {
      sel <- ThatRefi_k != 0
      for (i in 1:p) {
        ThatRefi_k[i, sel[i, ]] <- solve(t(X[, sel[i, ]]) %*% X[, sel[i, ]]) %*% t(X[, sel[i, ]]) %*% Y[, i]
      }
    }
    
    ThatRefi[, , k] <- ThatRefi_k
    
    if (intr) {
      constRefi[, k] <- Ybar - ThatRefi_k %*% Xbar
    }
    
    if (k == 1) {
      UpsOld <- ups_init
    } else {
      UpsOld <- UpsRefi[, , k - 1]
    }
    
    UpsNew <- UpsRefi[, , k]
    dUps <- UpsNew - UpsOld
    reldiffUps <- sqrt(sum(dUps^2)) / sqrt(sum(UpsOld^2))
    
    if (reldiffUps <= tol_ups) {
      break
    }
  }
  
  Ups <- UpsRefi[, , k]
  That <- ThatRefi[, , k]
  if (intr) {
    const <- Ybar - That %*% Xbar
  } else {
    const <- NULL
  }
  
  if (no_warn == 0 && k == K && reldiffUps > tol_ups) {
    warning("Maximum number of updates reached. Consider increasing K.")
    cat(sprintf("Relative change in penalty loadings is %3.1g percent\n", 100 * reldiffUps))
  }
  
  fit <- list(
    lambda = lambda_star,
    ups_init = ups_init,
    UpsRefi = UpsRefi,
    Ups = Ups,
    constInit = constInit,
    constRefi = constRefi,
    const = const,
    that_init = 
    ThatRefx =y
    Thax =yThat,
    kterm = k,
    reldiffUps = reldiffUps
  )
  
  return(fit)
}
