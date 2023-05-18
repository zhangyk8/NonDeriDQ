library(ks)
library(locpol)

DerivEstQuotient = function(Y, X, k, h=NULL, U=NULL){
  require(ks)
  n = length(X)
  ## Estimate the CDF of input X
  if(is.null(h)){
    cdf_est = kcde(X, eval.points = X, tail.flag = "lower.tail")
    den_est = kde(X, eval.points = X)
  }else{
    cdf_est = kcde(X, h=h, eval.points = X, tail.flag = "lower.tail")
    den_est = kde(X, h = h, eval.points = X)
  }
  if(is.null(U)){
    U = cdf_est$estimate
  }
  U_ord = sort(U)
  Y_ord = Y[order(U)]
  X_ord = X[order(U)]
  den_ord = den_est$estimate[order(U)]
  
  Y_deri = numeric(n)
  for (i in 1:n){
    if (i < k+1){
      if (i == 1){
        w1 = 0
      }else{
        w1 = U_ord[(i+1):(2*i-1)] - rev(U_ord[0:(i-1)])
      }
      w2 = U_ord[(2*i):k] - U_ord[i]
      w1 = w1/(sum(w1^2) + sum(w2^2))
      w2 = w2/(sum(w1^2) + sum(w2^2))
      Y_deri[i] = sum(w1*(Y_ord[(i+1):(2*i-1)] - rev(Y_ord[0:(i-1)]))) + sum(w2*(Y_ord[(2*i):k] - Y_ord[i]))
    }else if ((i >= k+1) & (i <= n-k)){
      ws = U_ord[(i+1):(i+k)] - rev(U_ord[(i-k):(i-1)])
      ws = ws/sum(ws^2)
      Y_deri[i] = sum(ws*(Y_ord[(i+1):(i+k)] - rev(Y_ord[(i-k):(i-1)])))
    }else{
      if (i == n){
        w1 = 0
        Y1 = 0
      }else{
        w1 = U_ord[(i+1):n] - rev(U_ord[(2*i-n):(i-1)])
        Y1 = Y_ord[(i+1):n] - rev(Y_ord[(2*i-n):(i-1)])
      }
      w2 = U_ord[i] - U_ord[(i-k):(2*i-n-1)]
      w1 = w1/(sum(w1^2) + sum(w2^2))
      w2 = w2/(sum(w1^2) + sum(w2^2))
      Y_deri[i] = sum(w1*Y1) + sum(w2*(Y_ord[i] - Y_ord[(i-k):(2*i-n-1)]))
    }
  }
  lst = list("Y_deri" = Y_deri, "U_ord" = U_ord, "X_ord" = X_ord, 
             "den_ord" = den_ord)
  return(lst)
}

ErrorVar = function(Y){
  n = length(Y)
  if(n<3){
    return(NA)
  }
  return(mean((0.809*Y[1:(n-2)] - 0.5*Y[2:(n-1)] - 0.309*Y[3:n])^2))
}

OptK = function(Y, U, k_range=NULL, h=NULL, kern=gaussK){
  n = length(U)
  if(is.null(k_range)){
    k_range = seq(1, floor((n-1)/2), by = 1)
  }
  u_qry = seq(min(U), max(U), length.out = 200)
  dat = data.frame(Y = Y, U = U)
  locpol_est = locpol(Y ~ U, data = dat, bw = h, kernel = kern, deg = 3, 
                      xeval = u_qry)
  B = max(abs(locpol_est$lpFit$Y2))
  sigma_e_sq = ErrorVar(Y)
  asym_bias_sq = 9*(B^2)*(k_range^2)*(k_range + 1)^2/(16*(n+1)^2*(2*k_range+1)^2)
  asym_var = 3*sigma_e_sq*(n+1)^2/(k_range*(k_range+1)*(2*k_range+1))
  return(k_range[which.min(asym_bias_sq + asym_var)])
  # lst = list("k_exact" = k_range[which.min(asym_bias_sq + asym_var)], 
  #            "k_asymp1" = floor(2^(4/5)*sigma_e_sq^(1/5)*B^(-2/5)*n^(4/5)),
  #            "k_asymp2" = floor(2*sigma_e_sq^(1/5)*B^(-2/5)*(n+1)^(4/5)))
  # return(lst)
}

LocalPolyReg = function(Y, X, xeval=NULL, h=0.1, kern=gaussK, deg=3, deriv_ord=0){
  require(pracma)
  if(is.null(xeval)){
    xeval = sapply(1:length(X), function(i) X[i])
  }
  Y_est = numeric(length(xeval))
  for(i in 1:length(xeval)){
    weights = kern((X - xeval[i])/h)
    # Filter out the data points with zero weights to speed up regressions with kernels of local support.
    # inds = which(abs(weights) > 1e-26)
    X_mat = matrix(0, nrow = length(X), ncol = deg+1)
    for(p in 0:deg){
      X_mat[,p+1] = (X - xeval[i])^p
    }
    W = diag(weights)
    beta = inv(t(X_mat) %*% W %*% X_mat) %*% t(X_mat) %*% W %*% Y
    Y_est[i] = factorial(deriv_ord)*beta[deriv_ord+1]
  }
  return(Y_est)
}

biGaussian = function(x){
  return((2/sqrt(pi))*x^2*exp(-x^2))
}

RSSBandwidth = function(Y, U, k, h_range=NULL, deg=3){
  n = length(U)
  if(is.null(h_range)){
    h_range = seq(0.001, 1, length.out = 100)
  }
  Y_int = Y[(k+1):(n-k)]
  U_int = U[(k+1):(n-k)]
  RSS = numeric(length(h_range))
  for(i in 1:length(h_range)){
    r_est = LocalPolyReg(Y = Y_int, X = U_int, xeval = U_int, h = h_range[i], 
                         kern = biGaussian, deg = deg, deriv_ord = 0)
    RSS[i] = mean((r_est - Y_int)^2)
  }
  return(h_range[which.min(RSS)])
}


SecondDerivEstQuotient = function(Y, X, k1, k2, h=NULL, U=NULL){
  require(ks)
  n = length(X)
  ## Estimate the CDF of input X
  if(is.null(h)){
    cdf_est = kcde(X, eval.points = X, tail.flag = "lower.tail")
    den_est = kde(X, eval.points = X)
    den_deriv = kdde(X, deriv.order = 1, eval.points = X)
  }else{
    cdf_est = kcde(X, h=h, eval.points = X, tail.flag = "lower.tail")
    den_est = kde(X, h=h, eval.points = X)
    den_deriv = kdde(X, h=h, deriv.order = 1, eval.points = X)
  }
  if(is.null(U)){
    U = cdf_est$estimate
  }
  U_ord = sort(U)
  Y_ord = Y[order(U)]
  X_ord = X[order(U)]
  den_ord = den_est$estimate[order(U)]
  den_deriv_ord = den_deriv$estimate[order(U)]
  
  weights = (2*seq(1, k2, by = 1) + k1)^2
  weights = weights/sum(weights)
  
  Y_aug = c(numeric(k1+k2), Y_ord, numeric(k1+k2))
  U_aug = c((1e-26)*seq(1, k1+k2, by = 1), U_ord, (1e-26)*seq(1, k1+k2, by = 1))
  
  Y_sec = numeric(n)
  for(i in 1:n){
    ind = i + k1+k2
    Y_plus = (Y_aug[(ind+1+k1):(ind+k1+k2)] - Y_aug[(ind+1):(ind+k2)])/(U_aug[(ind+1+k1):(ind+k1+k2)] - U_aug[(ind+1):(ind+k2)])
    Y_minus = rev(Y_aug[(ind-k1-k2):(ind-1-k1)] - Y_aug[(ind-k2):(ind-1)]) / rev(U_aug[(ind-k1-k2):(ind-1-k1)] - U_aug[(ind-k2):(ind-1)])
    U_norm = U_aug[(ind+1+k1):(ind+k1+k2)] + U_aug[(ind+1):(ind+k2)] - rev(U_aug[(ind-k1-k2):(ind-1-k1)]) - rev(U_aug[(ind-k2):(ind-1)])
    Y_sec[i] = sum(2*weights*((Y_plus - Y_minus) / U_norm))
  }
  lst = list("Y_sec" = Y_sec, "U_ord" = U_ord, "X_ord" = X_ord, 
             "den_ord" = den_ord, "den_deriv_ord" = den_deriv_ord)
  return(lst)
}

OptK1K2 = function(Y, U, k1_range=NULL, k2_range=NULL, h=NULL, kern=gaussK){
  n = length(U)
  if(is.null(k1_range)){
    k1_range = seq(1, floor((n-1)/4), by = 1)
  }
  if(is.null(k2_range)){
    k2_range = seq(1, floor((n-1)/4), by = 1)
  }
  u_qry = seq(min(U), max(U), length.out = 200)
  dat = data.frame(Y = Y, U = U)
  locpol_est = locpol(Y ~ U, data = dat, bw = h, kernel = kern, deg = 4, xeval = u_qry)
  B = max(abs(locpol_est$lpFit$Y3))
  sigma_e_sq = ErrorVar(Y)
  asym_bias_sq = matrix(0, nrow = length(k1_range), ncol = length(k2_range))
  asym_var = matrix(0, nrow = length(k1_range), ncol = length(k2_range))
  for(i in 1:length(k1_range)){
    k1 = k1_range[i]
    for(j in 1:length(k2_range)){
      k2 = k2_range[j]
      k2_cub = (k2^2)*((k2+1)^2)/4
      k2_sq = k2*(k2+1)*(2*k2+1)/6
      k2_sum = k2*(k2+1)/2
      asym_bias_sq[i,j] = ((B/(n+1))*(2*k2_cub + 3*k1*k2_sq + (5*k1^2/3)*k2_sum + k1^3*k2/3)/(4*k2_sq + k1*2*k2+4*k1*k2_sum))^2
      asym_var[i,j] = 4*(n+1)^4*sigma_e_sq/(k1^2*sum((2*seq(1, k2, by = 1) + k1)^2))
    }
  }
  asym_mse = asym_bias_sq + asym_var
  k_ind = which(asym_mse == min(asym_mse), arr.ind = TRUE)
  lst = list("k1" = k_ind[1], "k2" = k_ind[2])
  return(lst)
}
