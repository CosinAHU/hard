#@ function `generate.beta`
#@ returns the true value of beta
#@ params:
#@ weak: the magnitute of weak signal, can be either 0.05 or 0.1
#@ p: the dimension of beta
#@ q: number of repetitions
generate.beta <- function(weak, p, q){
  
  if (weak != 0.1 & weak != 0.05){
    stop("param `weak` must be one of 0.05 and 0.1.")
  }
  
  beta.strong <- c(0.6, 0, 0, -0.6, 0, 0)
  beta.weak <- c(weak, 0, 0, -weak, 0, 0)
  p0 <- length(beta.strong)
  
  if ( p < 6*p0 ){
    stop("The dimension p is too small, please enter a larger value.")
  }
  
  # beta0 is the true coefficient
  beta0 <- c(rep(c(beta.strong, beta.weak), q), rep(0, (p - 6 * p0)))
  
  return (beta0)
  
}


#@ function `generate.data`
#@ returns the data set described in simulation 1.
#@ params:
#@ weak: the magnitute of weak signal, can be either 0.05 or 0.1
#@ p: the dimension of the true coefficients.
#@ n: the sample size.
#@ sigma: the standard deviation of \epsilon \sim N(0, \sigma^2I_n).
#@ q: repeated times of v = (\beta_{strong}^{\top}, \beta_{weak}^{\top}).
generate.data <- function(n = 100, weak, p, sigma, q){
  
  require(MASS)
  require(Rcpp)
  require(RcppArmadillo)
  
  sourceCpp(
    code = '
      #include<Rmath.h>
      #include<RcppCommon.h>
      #include<RcppArmadillo.h>
      
      // [[Rcpp::depends(RcppArmadillo)]]
      using namespace std;
      using namespace Rcpp;
      using namespace arma;
      
      // [[Rcpp::export]]
      extern "C" SEXP gen_data(int n, double sigma, arma::mat cov_mat, arma::vec betat) {
        
        int p = cov_mat.n_rows;
        arma::mat y_train(n, n, fill::zeros);
        arma::cube X_train(n, p, n, fill::zeros);
        vec mean1(p, fill::zeros);
        vec mean2(n, fill::zeros);
        vec mean3(10000, fill::zeros);
        mat S1(n, n, fill::eye);
        mat S2(10000, 10000, fill::eye);
        
        for (int i = 0; i < n; i++){
          X_train.slice(i) = trans(mvnrnd(mean1, cov_mat, n));
          y_train.col(i) = X_train.slice(i) * betat + mvnrnd(mean2, pow(sigma, 2) * S1);
        }
        
        arma::mat X_test = trans(mvnrnd(mean1, cov_mat, 10000));
        arma::vec y_test = X_test * betat + mvnrnd(mean3, pow(sigma, 2) * S2);
        
        return Rcpp::List::create(
          Rcpp::Named("X.train") = X_train, 
          Rcpp::Named("y.train") = y_train,
          Rcpp::Named("X.test") = X_test,
          Rcpp::Named("y.test") = y_test
        );
      }
    '
    )

  beta0 <- generate.beta(weak, p, q)
  # Sigma is the cov matrix that the design matrices are sampled from
  Sigma <- matrix(0, p, p)
  for (i in 1:p){
    Sigma[i, ] <- 0.5^abs((i - 1:p))
  }
  
  return (gen_data(n, sigma, Sigma, beta0))
}


#@ function `metric`
#@ returns the several metrics for evaluating the performance of different regression methods.
#@ including PE, L2-loss, L1-loss, L_{\infty}-loss, FP, FN_strong, FN_weak, \hat{\sigma}.
#@ params:
#@ betah: the estimated coefficient
#@ betat: the true coefficient
#@ X.train: the training set X in a single loop
#@ y.train: the training set y in a single loop
#@ X.test: the test set X
#@ y.test: the test set y
metric <- function(betah, betat, X.train, y.train, X.test, y.test){
  
  if ( length(betah) != length(betat) ){
    stop("The length of estimated coefficient and true coefficient must be the same.")
  }
  
  betah <- matrix(betah, ncol=1)
  betat <- matrix(betat, ncol=1)
  
  PE <- mean((y.test - X.test %*% betah)^2)  #prediction error: PE
  l2 <- sqrt(sum((betah - betat)^2))  #L2-loss
  l1 <- sum(abs(betah - betat))  #L1-loss
  linf <- max(abs(betah - betat))  #L\infty-loss
  FP <- sum((betah !=0) * (betat == 0))
  
  #FN_strong
  idx.strong <- c(1:6, 13:18, 25:30)
  FN.strong <- sum(( betah[idx.strong] ==0 ) * (betat[idx.strong] != 0))
  
  #FN_weak
  idx.weak <- c(7:12, 19:24, 31:36)
  FN.weak <- sum(( betah[idx.weak] ==0 ) * (betat[idx.weak] != 0))
  
  #sigma.hat
  n <- nrow(y.train)
  p <- sum(betah != 0)
  e <- y.train - X.train %*% betah
  sigmah <- sqrt(sum(e^2) / abs((n - p)))
  
  return (list(PE = PE, l2.loss = l2, 
               l1.loss = l1, linf.loss = linf,
               FP = FP, FN.strong = FN.strong, 
               FN.weak = FN.weak, sigma = sigmah))
  
}
