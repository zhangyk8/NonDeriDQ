source('derivative_estimation_quotients.R')

jobid = commandArgs(TRUE)
print(jobid)

n = 700
set.seed(jobid)

## Simulation 3 (Figure 5)
X_3 = runif(n, min = 0.25, max = 1)
m_3 = sqrt(X_3*(1 - X_3))*sin((2.1*pi)/(X_3 + 0.05))
eps3 = rnorm(n, mean = 0, sd = 0.2)
Y_3 = m_3 + eps3

m3Deriv = function(x_qry){
  res = (1-2*x_qry)*sin(2.1*pi/(x_qry + 0.05))/(2*sqrt(x_qry*(1-x_qry))) - 2.1*pi*sqrt(x_qry*(1-x_qry))*cos(2.1*pi/(x_qry + 0.05))/((x_qry + 0.05)^2)
  return(res)
}

## Proposed estimator
U_3 = kcde(X_3, eval.points = X_3, tail.flag = "lower.tail")$estimate
k_opt3 = OptK(Y = Y_3, U = U_3, k_range = NULL, h = NULL, kern = gaussK)
# k_opt3 = 8
print(k_opt3)
deri_est3 = DerivEstQuotient(Y = Y_3, X = X_3, k = k_opt3)

h_opt3 = RSSBandwidth(Y = deri_est3$Y_deri, U = deri_est3$U_ord, k = k_opt3, 
                      h_range = seq(0.03, 0.07, by = 0.005), deg = 3)
print(h_opt3)
h_opt3 = h_opt3*1.01431

qry_range = seq(k_opt3 + 1, n-k_opt3, by = 1)
Y_deriv_sm3 = LocalPolyReg(Y = deri_est3$Y_deri[qry_range],
                           X = deri_est3$U_ord[qry_range],
                           xeval = deri_est3$U_ord, h = h_opt3,
                           kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm3 = LocalPolyReg(Y = deri_est3$Y_deri, X = deri_est3$U_ord, 
#                            xeval = deri_est3$U_ord, h = h_opt3, 
#                            kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_trans3 = Y_deriv_sm3 * deri_est3$den_ord

MAE_quo = mean(abs(Y_deriv_trans3 - m3Deriv(deri_est3$X_ord))[26:675])

Med_Abs_Err = median(abs(Y_deriv_trans3 - m3Deriv(deri_est3$X_ord)))


## Proposed estimator with the oracle data-generating distribution
U_3o = (4/3)*(X_3 - 1/4)
k_opt3o = OptK(Y = Y_3, U = U_3o, k_range = NULL, h = NULL, kern = gaussK)
# k_opt3o = 8
print(k_opt3o)
deri_est3o = DerivEstQuotient(Y = Y_3, X = X_3, k = k_opt3o, U = U_3o)

h_opt3o = RSSBandwidth(Y = deri_est3o$Y_deri, U = deri_est3o$U_ord, k = k_opt3o, 
                       h_range = seq(0.03, 0.07, by = 0.005), deg = 3)
print(h_opt3o)
h_opt3o = h_opt3o*1.01431

qry_rangeo = seq(k_opt3o + 1, n-k_opt3o, by = 1)
Y_deriv_sm3o = LocalPolyReg(Y = deri_est3o$Y_deri[qry_rangeo],
                            X = deri_est3o$U_ord[qry_rangeo], xeval = deri_est3o$U_ord,
                            h = h_opt3o, kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm3o = LocalPolyReg(Y = deri_est3o$Y_deri, 
#                             X = deri_est3o$U_ord, xeval = deri_est3o$U_ord, 
#                             h = h_opt3o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_trans3o = Y_deriv_sm3o * 4/3

MAE_quo_ora = mean(abs(Y_deriv_trans3o - m3Deriv(deri_est3o$X_ord))[26:675])

Med_Abs_Err_ora = median(abs(Y_deriv_trans3o - m3Deriv(deri_est3o$X_ord)))

deri_res3 = data.frame(jobid = rep(jobid, 2), k_opt = c(k_opt3, k_opt3o), 
                      MAE_adj = c(MAE_quo, MAE_quo_ora), 
                      Med_Abs_Err = c(Med_Abs_Err, Med_Abs_Err_ora), 
                      scenario = c("Proposed", "Oracle"))

write.csv(deri_res3, paste0("./deri_res/deri_first_order_", jobid, "_sim3_new.csv"), 
          row.names=FALSE)


## Simulation 4 (Figure 6)
n = 700
set.seed(jobid)
X_4 = rnorm(n, mean = 0, sd = 0.5)
m_4 = X_4 + 2*exp(-16*X_4^2)
eps4 = rnorm(n, mean = 0, sd = 0.2)
Y_4 = m_4 + eps4

m4Deriv = function(x){
  res = 1 - 64*x*exp(-16*x^2)
  return(res)
}

## Proposed estimator
U_4 = kcde(X_4, eval.points = X_4, tail.flag = "lower.tail")$estimate
k_opt4 = OptK(Y = Y_4, U = U_4, k_range = NULL, h = NULL, kern = gaussK)
# k_opt4 = 8
print(k_opt4)
deri_est4 = DerivEstQuotient(Y = Y_4, X = X_4, k = k_opt4)

h_opt4 = RSSBandwidth(Y = deri_est4$Y_deri, U = deri_est4$U_ord, k = k_opt4, 
                      h_range = seq(0.04, 0.08, by = 0.005), deg = 3)
print(h_opt4)
h_opt4 = h_opt4*1.01431

qry_range = seq(k_opt4 + 1, n-k_opt4, by = 1)
# Y_deriv_sm4 = LocalPolyReg(Y = deri_est4$Y_deri[qry_range], X = deri_est4$U_ord[qry_range], 
#                            xeval = deri_est4$U_ord, h = h_opt4, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_sm4 = LocalPolyReg(Y = deri_est4$Y_deri, X = deri_est4$U_ord, 
                           xeval = deri_est4$U_ord, h = h_opt4, kern = gaussK, deg = 3, 
                           deriv_ord = 0)
Y_deriv_trans4 = Y_deriv_sm4 * deri_est4$den_ord

MAE_quo = mean(abs(Y_deriv_trans4 - m4Deriv(deri_est4$X_ord))[26:675])

Med_Abs_Err = median(abs(Y_deriv_trans4 - m4Deriv(deri_est4$X_ord)))

## Proposed estimator with the oracle data-generating distribution
U_4o = pnorm(X_4, mean = 0, sd = 0.5)
k_opt4o = OptK(Y = Y_4, U = U_4o, k_range = NULL, h = NULL, kern = gaussK)
# k_opt4o = 8
print(k_opt4o)
deri_est4o = DerivEstQuotient(Y = Y_4, X = X_4, k = k_opt4o, U = U_4o)

h_opt4o = RSSBandwidth(Y = deri_est4o$Y_deri, U = deri_est4o$U_ord, k = k_opt4o, 
                       h_range = seq(0.04, 0.08, by = 0.005), deg = 3)
print(h_opt4o)
h_opt4o = h_opt4o*1.01431

qry_rangeo = seq(k_opt4o + 1, n-k_opt4o, by = 1)
Y_deriv_sm4o = LocalPolyReg(Y = deri_est4o$Y_deri[qry_rangeo],
                            X = deri_est4o$U_ord[qry_rangeo], xeval = deri_est4o$U_ord,
                            h = h_opt4o, kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm4o = LocalPolyReg(Y = deri_est4o$Y_deri, 
#                             X = deri_est4o$U_ord, xeval = deri_est4o$U_ord, 
#                             h = h_opt4o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_trans4o = Y_deriv_sm4o * dnorm(deri_est4o$X_ord, mean = 0, sd = 0.5)


MAE_quo_ora = mean(abs(Y_deriv_trans4o - m4Deriv(deri_est4o$X_ord))[26:675])

Med_Abs_Err_ora = median(abs(Y_deriv_trans4o - m4Deriv(deri_est4o$X_ord)))

deri_res4 = data.frame(jobid = rep(jobid, 2), k_opt = c(k_opt4, k_opt4o), 
                      MAE_adj = c(MAE_quo, MAE_quo_ora), 
                      Med_Abs_Err = c(Med_Abs_Err, Med_Abs_Err_ora),
                      scenario = c("Proposed", "Oracle"))

write.csv(deri_res4, paste0("./deri_res/deri_first_order_", jobid, "_sim4_new.csv"), 
          row.names=FALSE)

