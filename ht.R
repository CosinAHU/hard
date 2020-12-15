# @ function `ht`
# @
# @ It computes a minimizer of the hard-thresholded least-squares problem:
# @
# @ min_beta  (2 n)^(-1) ||y - X beta||_2^2 + lambda ||p_lambda0(|beta|)||_1,
# @ params:
# @ y is an n-vector of response
# @ X is an n x p design matrix with each column vector rescaled to have L2-norm n^{1/2}
# @ lambda >= 0 is the regularization parameter
# @ lambda0 is the parameter of penalty function.
# @ varset is the potential important variables based on prior knowledge, default empty set.
# @ inival is the initial value of beta. In order to make it stable in numerical computation, if it is
# @ not specified, a random initialization using normal random numbers will be used.
# @ maxiter is the maximum number of iteration
# @ tol is the ||\hat{\beta}^{(t+1)} - \hat{\beta}^{(t)}||_2, defines the stopping criterion
# @ returns: the hard-thresholding estimate of \beta


ht <- function(X, y, lambda0 = 1e-3, lambda = 1e-2, varset = c(), inival = integer(), maxiter = 50, tol = 1e-4) {
  
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(inival)){
    inival <- rnorm(p)
  }
  XXmat <- 1/n * t(X) %*% X
  cvec <- 1/n * t(X) %*% y
  
  beta <- inival
  iter <- 1
  update <- 1
  ind <- 1:p
  varset <- union(which(inival != 0), varset)
  
  while ( (iter <= maxiter) & (update > tol) ){
    iter <- iter + 1
    betaold <- beta
    
    for (k in 1:length(ind)){
      setr <- ind
      I <- setr[k]
      setr <- setr[-k]
      
      z1 <- XXmat[I, I]
      Lam <- 1/z1
      if (length(setr) == 0) {
        z <- cvec[I]
      } else {            
        z <- (cvec[I] - XXmat[I, setr]%*%beta[setr])
      }
      
      beta[I] <- uhard(z, Lam, lambda0, lambda)
    }
    
    ind <- which(beta != 0)
    update <- sqrt(sum((beta - betaold)^2))
    
    setr <- setdiff(1:p, ind)
    
    if (length(setr) == 0){
      resc <- 0
    }else{
      resc <- abs(cvec[setr] - XXmat[setr, ind]%*%beta[ind])
    }
    
    indm <- which(resc > lambda0 + lambda)
    ind <- union(c(t(setr[indm]), ind), varset)
  }
  
  return (beta)
}

  
#@ function `uhard`
#@ univariate minimization of
#@ min_beta  (2 Lam)^(-1) (z - beta)^2 + lambda0 |beta| + p_lambda(|beta|) for hard thresholding penalty

uhard <- function(z, Lam = 1, lambda0, lambda){
  
  z <- sign(z) * max(0, abs(z) - Lam*lambda0)
  if ( Lam == 1 ){
    beta <- z*(abs(z) > lambda)
  }else{
    z0 <- sign(z) * max(0, abs(z) - Lam*lambda) / (1 - Lam)
    if (abs(z) <= lambda){
      beta <- z0
    }else if ( (abs(z) <= Lam*lambda) & ((z - z0)^2 / 2 + Lam*hard.thred(abs(z0), lambda) <= Lam*hard.thred(abs(z), lambda)) ){
      beta <- z0
    }else{
      beta <- z
    }
  }
  
  return (beta)
}


#@ function `hard.thred`
#@ hard thresholding function p_{H, \lambda}(t) = \frac{1}{2} (\lambda^2 - (\lambda - t)_{+}^2)
#@ params:
#@ t: a non-negative scalar
#@ lambda: penalty coefficient
#@ returns:
#@ the penalty function value
hard.thred <- function(t, lambda){
  
  return (1/2 * (lambda^2 - pmax(0, lambda - t)^2))
  
}

