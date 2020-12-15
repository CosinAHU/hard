################################################################
###########################simulation 1#########################
################################################################
source("./util.R")
source("./sica.R") # for SICA
source("./ht.R") # for hard-thresholding
source("./oracle.R") # for oracle procedure
library(glmnet) # for lasso


#@ (p, sigma, weak) = (1000, 0.4, 0.05)
set.seed(2020)
data1 <- generate.data(n=100, weak=0.05, p=1000, sigma=0.4, q=3)
X.train <- data1[["X.train"]]
X.test <- data1[["X.test"]]
y.train <- data1[["y.train"]]
y.test <- data1[["y.test"]]
eval.lasso <- eval.ht <- eval.sica <- eval.oracle <- matrix(0, nrow=100, ncol=8)
beta.true <- generate.beta(weak=0.05, p=1000, q=3)
beta.true <- matrix(beta.true, ncol=1)
for (iter in 1:100){
  cat("---------This is the", iter, "th iteration----------\n")
  Xtrain <- X.train[, , iter]
  ytrain <- as.matrix(y.train[, iter])
  lam <- cv.glmnet(Xtrain, ytrain, type.measure = "mse")$lambda.min
  beta.lasso <- matrix(glmnet(Xtrain, ytrain, family = "gaussian", alpha=1, lambda=lam)$beta)
  beta.ht <- ht(Xtrain ,ytrain)
  beta.sica <- t(sica(Xtrain, ytrain))
  beta.oracle <- oracle(Xtrain, ytrain, beta.true)
  eval.lasso[iter, ] <- unlist(metric(beta.lasso, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.ht[iter, ] <- unlist(metric(beta.ht, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.sica[iter, ] <- unlist(metric(beta.sica, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.oracle[iter, ] <- unlist(metric(beta.oracle, beta.true, Xtrain, ytrain, X.test, y.test))
}

out1 <- apply(eval.lasso, 2, function(x) {c(mean(x[!is.infinite(x)]), sd(x[!is.infinite(x)]))})
out2 <- apply(eval.ht, 2, function(x) {c(mean(x), sd(x))})
out3 <- apply(eval.sica, 2, function(x) {c(mean(x), sd(x))})
out4 <- apply(eval.oracle, 2, function(x) {c(mean(x), sd(x))})
out <- cbind(t(out1), t(out2), t(out3), t(out4))
write.csv(out, file="./p=1000 weak=0.05.csv", row.names = FALSE)

#@ (p, sigma, weak) = (1000, 0.4, 0.1)
data2 <- generate.data(n=100, weak=0.1, p=1000, sigma=0.4, q=3)
X.train <- data2[["X.train"]]
X.test <- data2[["X.test"]]
y.train <- data2[["y.train"]]
y.test <- data2[["y.test"]]
eval.lasso <- eval.ht <- eval.sica <- eval.oracle <- matrix(0, nrow=100, ncol=8)
beta.true <- generate.beta(weak=0.1, p=1000, q=3)
beta.true <- matrix(beta.true, ncol=1)
for (iter in 1:100){
  cat("---------This is the", iter, "th iteration----------\n")
  Xtrain <- X.train[, , iter]
  ytrain <- as.matrix(y.train[, iter])
  lam <- cv.glmnet(Xtrain, ytrain, type.measure = "mse")$lambda.min
  beta.lasso <- matrix(glmnet(Xtrain, ytrain, family = "gaussian", alpha=1, lambda=lam)$beta)
  beta.ht <- ht(Xtrain ,ytrain, lambda = 100)
  beta.sica <- t(sica(Xtrain, ytrain))
  beta.oracle <- oracle(Xtrain, ytrain, beta.true)
  eval.lasso[iter, ] <- unlist(metric(beta.lasso, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.ht[iter, ] <- unlist(metric(beta.ht, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.sica[iter, ] <- unlist(metric(beta.sica, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.oracle[iter, ] <- unlist(metric(beta.oracle, beta.true, Xtrain, ytrain, X.test, y.test))
}

out1 <- apply(eval.lasso, 2, function(x) {c(mean(x[!is.infinite(x)]), sd(x[!is.infinite(x)]))})
out2 <- apply(eval.ht, 2, function(x) {c(mean(x), sd(x))})
out3 <- apply(eval.sica, 2, function(x) {c(mean(x), sd(x))})
out4 <- apply(eval.oracle, 2, function(x) {c(mean(x), sd(x))})
out <- cbind(t(out1), t(out2), t(out3), t(out4))
write.csv(out, file="./p=1000 weak=0.1.csv", row.names = FALSE)


#@ (p, sigma, weak) = (5000, 0.3, 0.05)
data3 <- generate.data(n=100, weak=0.05, p=5000, sigma=0.3, q=3)
X.train <- data3[["X.train"]]
X.test <- data3[["X.test"]]
y.train <- data3[["y.train"]]
y.test <- data3[["y.test"]]
eval.lasso <- eval.ht <- eval.sica <- eval.oracle <- matrix(0, nrow=100, ncol=8)
beta.true <- generate.beta(weak=0.05, p=5000, q=3)
beta.true <- matrix(beta.true, ncol=1)
for (iter in 1:100){
  cat("---------This is the", iter, "th iteration----------\n")
  Xtrain <- X.train[, , iter]
  ytrain <- as.matrix(y.train[, iter])
  lam <- cv.glmnet(Xtrain, ytrain, type.measure = "mse")$lambda.min
  beta.lasso <- matrix(glmnet(Xtrain, ytrain, family = "gaussian", alpha=1, lambda=lam)$beta)
  beta.ht <- ht(Xtrain ,ytrain, lambda = 100)
  beta.sica <- t(sica(Xtrain, ytrain))
  beta.oracle <- oracle(Xtrain, ytrain, beta.true)
  eval.lasso[iter, ] <- unlist(metric(beta.lasso, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.ht[iter, ] <- unlist(metric(beta.ht, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.sica[iter, ] <- unlist(metric(beta.sica, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.oracle[iter, ] <- unlist(metric(beta.oracle, beta.true, Xtrain, ytrain, X.test, y.test))
}

out1 <- apply(eval.lasso, 2, function(x) {c(mean(x[!is.infinite(x)]), sd(x[!is.infinite(x)]))})
out2 <- apply(eval.ht, 2, function(x) {c(mean(x), sd(x))})
out3 <- apply(eval.sica, 2, function(x) {c(mean(x), sd(x))})
out4 <- apply(eval.oracle, 2, function(x) {c(mean(x), sd(x))})
out <- cbind(t(out1), t(out2), t(out3), t(out4))
write.csv(out, file="./p=5000 weak=0.05.csv", row.names = FALSE)


#@ (p, sigma, weak) = (5000, 0.3, 0.1)
data4 <- generate.data(n=100, weak=0.1, p=5000, sigma=0.3, q=3)
X.train <- data4[["X.train"]]
X.test <- data4[["X.test"]]
y.train <- data4[["y.train"]]
y.test <- data4[["y.test"]]
eval.lasso <- eval.ht <- eval.sica <- eval.oracle <- matrix(0, nrow=100, ncol=8)
beta.true <- generate.beta(weak=0.1, p=5000, q=3)
beta.true <- matrix(beta.true, ncol=1)
for (iter in 1:100){
  cat("---------This is the", iter, "th iteration----------\n")
  Xtrain <- X.train[, , iter]
  ytrain <- as.matrix(y.train[, iter])
  lam <- cv.glmnet(Xtrain, ytrain, type.measure = "mse")$lambda.min
  beta.lasso <- matrix(glmnet(Xtrain, ytrain, family = "gaussian", alpha=1, lambda=lam)$beta)
  beta.ht <- ht(Xtrain ,ytrain, lambda = 100)
  beta.sica <- t(sica(Xtrain, ytrain))
  beta.oracle <- oracle(Xtrain, ytrain, beta.true)
  eval.lasso[iter, ] <- unlist(metric(beta.lasso, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.ht[iter, ] <- unlist(metric(beta.ht, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.sica[iter, ] <- unlist(metric(beta.sica, beta.true, Xtrain, ytrain, X.test, y.test))
  eval.oracle[iter, ] <- unlist(metric(beta.oracle, beta.true, Xtrain, ytrain, X.test, y.test))
}

out1 <- apply(eval.lasso, 2, function(x) {c(mean(x[!is.infinite(x)]), sd(x[!is.infinite(x)]))})
out2 <- apply(eval.ht, 2, function(x) {c(mean(x), sd(x))})
out3 <- apply(eval.sica, 2, function(x) {c(mean(x), sd(x))})
out4 <- apply(eval.oracle, 2, function(x) {c(mean(x), sd(x))})
out <- cbind(t(out1), t(out2), t(out3), t(out4))
write.csv(out, file="./p=5000 weak=0.1.csv", row.names = FALSE)


#########################################################################
###########################Simulation 2##################################
#########################################################################

#@ compute the L2-loss of Lasso, hard, SICA and oracle procedure
#@ params:
#@ X.train: a tensor of shape (n, p, n)
#@ y.train: a matrix of shape (n, n)
#@ X.test: a matrix of 10000 samples
#@ y.test: a matrix of 10000 by 1
#@ lambda: ridge parameter
#@ returns: l2-loss given a fixed lambda
#@ note: in this simulation, we only consider (n, p, weak) = (100, 1000, 0.05)
L2.risk <- function(X.train, y.train, X.test, y.test, lambda, betat){
  n <- ncol(y.train)
  l2.loss <- matrix(0, nrow = n, ncol = 4)
  for (iter in 1:n){
    cat("---------This is the", iter, "th iteration----------\n")
    Xtrain <- X.train[, , iter]
    ytrain <- as.matrix(y.train[, iter])
    lam <- cv.glmnet(Xtrain, ytrain, type.measure = "mse")$lambda.min
    beta.lasso <- matrix(glmnet(Xtrain, ytrain, family = "gaussian", alpha=1, lambda=lam)$beta)
    beta.ht <- ht(Xtrain ,ytrain)
    beta.sica <- t(sica(Xtrain, ytrain))
    beta.oracle <- oracle(Xtrain, ytrain, betat)
    
    idx.lasso <- which(beta.lasso != 0)
    s <- length(idx.lasso)
    X.refit <- Xtrain[, idx.lasso]
    beta.refit <- solve(t(X.refit) %*% X.refit + lambda * diag(s)) %*% t(X.refit) %*% ytrain
    y.hat <- X.test[, idx.lasso] %*% beta.refit
    e <- y.test - y.hat
    loss <- sum(e^2)
    l2.loss[iter, 1] <- loss
    
    idx.ht <- which(beta.ht != 0)
    s <- length(idx.ht)
    X.refit <- Xtrain[, idx.ht]
    beta.refit <- solve(t(X.refit) %*% X.refit + lambda * diag(s)) %*% t(X.refit) %*% ytrain
    y.hat <- X.test[, idx.ht] %*% beta.refit
    e <- y.test - y.hat
    loss <- sum(e^2)
    l2.loss[iter, 2] <- loss
    
    idx.sica <- which(beta.sica != 0)
    s <- length(idx.sica)
    X.refit <- Xtrain[, idx.sica]
    beta.refit <- solve(t(X.refit) %*% X.refit + lambda * diag(s)) %*% t(X.refit) %*% ytrain
    y.hat <- X.test[, idx.sica] %*% beta.refit
    e <- y.test - y.hat
    loss <- sum(e^2)
    l2.loss[iter, 3] <- loss
    
    idx.oracle <- which(beta.oracle != 0)
    s <- length(idx.oracle)
    X.refit <- Xtrain[, idx.oracle]
    beta.refit <- solve(t(X.refit) %*% X.refit + lambda * diag(s)) %*% t(X.refit) %*% ytrain
    y.hat <- X.test[, idx.oracle] %*% beta.refit
    e <- y.test - y.hat
    loss <- sum(e^2)
    l2.loss[iter, 4] <- loss
    
  }
  return (l2.loss)
}

lam <- seq(0, 8, by=0.1)
data <- data1
X.train <- data[["X.train"]]
X.test <- data[["X.test"]]
y.train <- data[["y.train"]]
y.test <- data[["y.test"]]
beta.true <- generate.beta(weak=0.05, p=1000, q=3)

loss <- array(0, dim=c(100, 4, length(lam)))
for (k in 1:length(lam)){
  cat("Processing lambda[", k, "] = ", lam[k], "\n")
  loss[, , k] <- L2.risk(X.train, y.train, X.test, y.test, lam[k], beta.true)
}

ave.loss <- apply(loss, 2:3, mean) / 10000
par(mfrow=c(2,2))
for (i in 1:4){
  plot(lam, ave.loss[i, ], xlab=quote(lambda[1]), ylab="Prediction risk", type="l")
} 
