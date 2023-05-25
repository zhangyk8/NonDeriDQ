source('derivative_estimation_quotients.R')
library(pspline)

jobid = commandArgs(TRUE)
print(jobid)

n = 1000
set.seed(jobid)

## Simulation 1 (Figure 12)
X_1 = runif(n)
m_1 = cos(2*pi*X_1)^2 + log(4/3 + X_1)
eps1 = rnorm(n, mean = 0, sd = 0.1)
Y_1 = m_1 + eps1

m1Deriv = function(x_qry){
  res = -2*pi*sin(4*pi*x_qry) + 3/(3*x_qry + 4)
  return(res)
}

## Proposed estimator
U_1 = kcde(X_1, eval.points = X_1, tail.flag = "lower.tail")$estimate
k_opt1 = OptK(Y = Y_1, U = U_1, k_range = NULL, h = NULL, kern = gaussK)
# k_opt1 = 8
print(k_opt1)
deri_est1 = DerivEstQuotient(Y = Y_1, X = X_1, k = k_opt1)

h_opt1 = RSSBandwidth(Y = deri_est1$Y_deri, U = deri_est1$U_ord, k = k_opt1, 
                      h_range = seq(0.03, 0.07, by = 0.005), deg = 3)
print(h_opt1)
h_opt1 = h_opt1*1.01431

qry_range = seq(k_opt1 + 1, n-k_opt1, by = 1)
# Y_deriv_sm1 = LocalPolyReg(Y = deri_est1$Y_deri[qry_range], 
#                            X = deri_est1$U_ord[qry_range], 
#                            xeval = deri_est1$U_ord, h = h_opt1, 
#                            kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_sm1 = LocalPolyReg(Y = deri_est1$Y_deri, 
                           X = deri_est1$U_ord, 
                           xeval = deri_est1$U_ord, h = h_opt1, 
                           kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_trans1 = Y_deriv_sm1 * deri_est1$den_ord

MAE_quo = mean(abs(Y_deriv_trans1 - m1Deriv(deri_est1$X_ord))[26:675])
Med_Abs_Err = median(abs(Y_deriv_trans1 - m1Deriv(deri_est1$X_ord)))


## Proposed estimator with the oracle data-generating distribution
U_1o = X_1
k_opt1o = OptK(Y = Y_1, U = U_1o, k_range = NULL, h = NULL, kern = gaussK)
# k_opt1o = 8
print(k_opt1o)
deri_est1o = DerivEstQuotient(Y = Y_1, X = X_1, k = k_opt1o, U = U_1o)

h_opt1o = RSSBandwidth(Y = deri_est1o$Y_deri, U = deri_est1o$U_ord, k = k_opt1o, 
                       h_range = seq(0.03, 0.07, by = 0.005), deg = 3)
print(h_opt1o)
h_opt1o = h_opt1o*1.01431

qry_rangeo = seq(k_opt1o + 1, n-k_opt1o, by = 1)
# Y_deriv_sm1o = LocalPolyReg(Y = deri_est1o$Y_deri[qry_rangeo], 
#                             X = deri_est1o$U_ord[qry_rangeo], xeval = deri_est1o$U_ord,
#                             h = h_opt1o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_sm1o = LocalPolyReg(Y = deri_est1o$Y_deri, 
                            X = deri_est1o$U_ord, xeval = deri_est1o$U_ord,
                            h = h_opt1o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_trans1o = Y_deriv_sm1o

MAE_quo_ora = mean(abs(Y_deriv_trans1o - m1Deriv(deri_est1o$X_ord))[26:675])
Med_Abs_Err_ora = median(abs(Y_deriv_trans1o - m1Deriv(deri_est1o$X_ord)))


## Adaptive spline estimator
Y_pspline1 = sm.spline(x = sort(X_1), y = Y_1[order(X_1)], norder = 2)
Y_psp1 = predict(Y_pspline1, xarg = sort(X_1), nderiv = 0)
X_1_ord = sort(X_1)
gam_deriv = (Y_psp1[2:n] - Y_psp1[1:(n-1)])/(X_1_ord[2:n] - X_1_ord[1:(n-1)])
gam_deriv = (c(0, gam_deriv) + c(gam_deriv, 0))/2

MAE_sp = mean(abs(gam_deriv - m1Deriv(X_1_ord))[26:675])
Med_Abs_Err_sp = median(abs(gam_deriv - m1Deriv(X_1_ord)))

deri_res1 = data.frame(jobid = rep(jobid, 3), k_opt = c(k_opt1, k_opt1o, NA), 
                       MAE_adj = c(MAE_quo, MAE_quo_ora, MAE_sp), 
                       Med_Abs_Err = c(Med_Abs_Err, Med_Abs_Err_ora, Med_Abs_Err_sp),
                       scenario = c("Proposed", "Oracle", "GAM"))

write.csv(deri_res1, paste0("./deri_res/deri_first_order_", jobid, "_sim1.csv"), 
          row.names=FALSE)

