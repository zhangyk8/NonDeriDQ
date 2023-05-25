library(locpol)

GasserMullerEst = function(Y, X, xeval=NULL, h=NULL, kern=gaussK, deriv_ord=1){
  # From the paper "Estimating Regression Functions and Their Derivatives by the Kernel Method"
  # When "deriv_ord=0", "kern" should be the integration of the kernel function, e.g., "pnorm"; 
  # When "deriv_ord=1", "kern" should be the kernel function, e.g., "gaussK";
  # When "deriv_ord=2", "kern" should be the derivative of the kernel function.
  n = length(X)
  if(is.null(h)){
    require(locpol)
    h = regCVBwSelC(X, Y, deg = 0, kernel = gaussK)
    # h = thumbBw(X, Y, deg = 0, kernel = gaussK)
  }
  if(is.null(xeval)){
    xeval = sapply(1:length(X), function(i) X[i])
  }
  X_ord = sort(X)
  X_ord = c(c(-100000, X_ord), 100000)
  s_val = (X_ord[-1] + X_ord[1:(n+1)])/2
  est = sapply(xeval, function(x){
    sum(Y*(kern((x - s_val[1:n])/h) - kern((x - s_val[2:(n+1)])/h)))/(h^deriv_ord)
  })
  return(est)
}

gaussDeriv = function(x){
  -(x/sqrt(2*pi))*exp(-x^2/2)
}

gaussDeriv2 = function(x){
  ((x^1-1)/sqrt(2*pi))*exp(-x^2/2)
}

NWVarEst = function(Y, X, xeval=NULL, h_reg=NULL, h_den=NULL, kern=gaussDeriv, deriv_ord=1){
  # From the paper "Derivative Estimation in Nonparametric Regression with Random Predictor Variable".
  # When "deriv_ord=0", we implement the classical Nadaraya-Watson estimator with kernel "kern"; 
  # When "deriv_ord>=1", "kern" should be the derivatives of the kernel function with corresponding order.
  require(ks)
  n = length(X)
  if(is.null(h_reg)){
    require(locpol)
    h_reg = regCVBwSelC(X, Y, deg = 0, kernel = gaussK)
  }
  if(is.null(h_den)){
    h_den = hpi(X)
  }
  if(is.null(xeval)){
    xeval = sapply(1:length(X), function(i) X[i])
  }
  est = sapply(xeval, function(x){
    if(deriv_ord == 0){
      sum(Y*kern((x-X)/h_reg))/sum(kern((x-X)/h_reg))
    }
    else{
      den_est = kde(X, h = h_den, eval.points = X)$estimate
      sum(Y*kern((x-X)/h_reg) / den_est)/ (n*h_reg^(deriv_ord+1))
    }
  })
  return(est)
}
