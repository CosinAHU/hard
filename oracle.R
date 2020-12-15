#@ oracle procedure
#@ params:
#@ X: the design matrix, matrix
#@ y: the response vector, matrix
#@ betat: the true coefficient, matrix
#@ returns:
#@ the estimated coefficient using oracle procedure

oracle <- function(X, y, betat){
  idx <- which( betat != 0 )
  X.true <- X[, idx]
  beta.or <- solve(t(X.true) %*% X.true) %*% t(X.true) %*% y
  beta.oracle <- matrix(0, nrow=length(betat), ncol=1)
  beta.oracle[idx] <- beta.or
  return (beta.oracle)
}