source('derivative_estimation_quotients.R')

jobid = commandArgs(TRUE)
print(jobid)

n = 700
set.seed(jobid)
X_6 = runif(n, min = 0, max = 1)
m_6 = 8*exp(-(1-5*X_6)^3*(1-7*X_6))
eps6 = rnorm(n, mean = 0, sd = 0.1)
Y_6 = m_6 + eps6

m6SecDeriv = function(x){
  res = 32*((70*x - 11)^2*(5*x-1)^3 - 525*x + 90)*(5*x - 1)*exp(-(1-5*x)^3*(1-7*x))
  return(res)
}

## Proposed estimator
U_6 = kcde(X_6, eval.points = X_6, tail.flag = "lower.tail")$estimate
k_opt6 = OptK(Y = Y_6, U = U_6, k_range = seq(1, 300, by = 1), h = NULL, kern = gaussK)
print(k_opt6)
deri_est6 = DerivEstQuotient(Y = Y_6, X = X_6, k = k_opt6)

k_opt_2d6 = OptK1K2(Y = Y_6, U = U_6, k1_range = seq(1, 100, by = 1), 
                    k2_range = seq(1, 100, by = 1), h = NULL, kern = gaussK)
# k_opt_2d6 = list("k1" = 15, "k2" = 8)
print(k_opt_2d6)
deri_sec6 = SecondDerivEstQuotient(Y = Y_6, X = X_6, k1 = k_opt_2d6$k1, 
                                   k2 = k_opt_2d6$k2)

h_opt6 = RSSBandwidth(Y = deri_est6$Y_deri, U = deri_est6$U_ord, k = k_opt6, 
                      h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt6)
h_opt6 = h_opt6*1.01431
h_opt_2d6 = RSSBandwidth(Y = deri_sec6$Y_sec, U = deri_sec6$U_ord, k = k_opt_2d6$k1 + k_opt_2d6$k2, 
                         h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt_2d6)
h_opt_2d6 = h_opt_2d6*1.01431

qry_range = seq(k_opt_2d6$k1 + k_opt_2d6$k2 + 1, 700 - k_opt_2d6$k1 - k_opt_2d6$k2, by = 1)
Y_deriv_sm6 = LocalPolyReg(Y = deri_est6$Y_deri[qry_range],
                           X = deri_est6$U_ord[qry_range], xeval = deri_est6$U_ord,
                           h = h_opt6, kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_sm6 = LocalPolyReg(Y = deri_sec6$Y_sec[qry_range], X = deri_sec6$U_ord[qry_range],
                         xeval = deri_sec6$U_ord, h = h_opt_2d6,
                         kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm6 = LocalPolyReg(Y = deri_est6$Y_deri, 
#                            X = deri_est6$U_ord, xeval = deri_est6$U_ord, 
#                            h = h_opt6, kern = gaussK, deg = 3, deriv_ord = 0)
# Y_sec_sm6 = LocalPolyReg(Y = deri_sec6$Y_sec, X = deri_sec6$U_ord, 
#                          xeval = deri_sec6$U_ord, h = h_opt_2d6, 
#                          kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_trans6 = Y_sec_sm6 * deri_sec6$den_ord^2 + Y_deriv_sm6 * deri_sec6$den_deriv_ord

MAE_quo = mean(abs(Y_sec_trans6 - m6SecDeriv(deri_sec6$X_ord))[31:670])
Med_Abs_Err = median(abs(Y_sec_trans6 - m6SecDeriv(deri_sec6$X_ord)))


## Proposed estimator with the oracle data-generating distribution
U_6o = X_6
k_opt6o = OptK(Y = Y_6, U = U_6o, k_range = NULL, h = NULL, kern = gaussK)
print(k_opt6o)
deri_est6o = DerivEstQuotient(Y = Y_6, X = X_6, k = k_opt6o, U = U_6o)
k_opt_2d6o = OptK1K2(Y = Y_6, U = U_6o, k1_range = seq(1, 100, by = 1), 
                     k2_range = seq(1, 100, by = 1), h = NULL, kern = gaussK)
# k_opt_2d6 = list("k1" = 15, "k2" = 8)
print(k_opt_2d6o)
deri_sec6o = SecondDerivEstQuotient(Y = Y_6, X = X_6, k1 = k_opt_2d6o$k1, 
                                    k2 = k_opt_2d6o$k2, U = U_6o)

h_opt6o = RSSBandwidth(Y = deri_est6o$Y_deri, U = deri_est6o$U_ord, k = k_opt6o, 
                       h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt6o)
h_opt6o = h_opt6o*1.01431
h_opt_2d6o = RSSBandwidth(Y = deri_sec6o$Y_sec, U = deri_sec6o$U_ord, 
                          k = k_opt_2d6o$k1 + k_opt_2d6o$k2, 
                          h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt_2d6o)
h_opt_2d6o = h_opt_2d6o*1.01431

qry_range_o = seq(k_opt_2d6o$k1 + k_opt_2d6o$k2 + 1, 700 - k_opt_2d6o$k1 - k_opt_2d6o$k2, by = 1)
Y_deriv_sm6o = LocalPolyReg(Y = deri_est6o$Y_deri[qry_range_o],
                            X = deri_est6o$U_ord[qry_range_o],
                            xeval = deri_est6o$U_ord, h = h_opt6o,
                            kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_sm6o = LocalPolyReg(Y = deri_sec6o$Y_sec[qry_range_o],
                          X = deri_sec6o$U_ord[qry_range_o], xeval = deri_est6o$U_ord,
                          h = h_opt_2d6o, kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm6o = LocalPolyReg(Y = deri_est6o$Y_deri, 
#                             X = deri_est6o$U_ord, 
#                             xeval = deri_est6o$U_ord, h = h_opt6o, 
#                             kern = gaussK, deg = 3, deriv_ord = 0)
# Y_sec_sm6o = LocalPolyReg(Y = deri_sec6o$Y_sec, 
#                           X = deri_sec6o$U_ord, xeval = deri_est6o$U_ord, 
#                           h = h_opt_2d6o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_trans6o = Y_sec_sm6o + Y_deriv_sm6o * 0

MAE_quo_ora = mean(abs(Y_sec_trans6o - m6SecDeriv(deri_sec6o$X_ord))[31:670])
Med_Abs_Err_ora = median(abs(Y_sec_trans6o - m6SecDeriv(deri_est6o$X_ord)))

deri_sec_res6 = data.frame(jobid = rep(jobid, 2), k_opt = c(k_opt6, k_opt6o), 
                       k1_opt = c(k_opt_2d6$k1, k_opt_2d6o$k1),
                       k2_opt = c(k_opt_2d6$k2, k_opt_2d6o$k2),
                       MAE_adj = c(MAE_quo, MAE_quo_ora), 
                       Med_Abs_Err = c(Med_Abs_Err, Med_Abs_Err_ora), 
                       scenario = c("Proposed", "Oracle"))

write.csv(deri_sec_res6, paste0("./deri_res/deri_second_order_", jobid, "_sim6_new.csv"), 
          row.names=FALSE)


## Non-uniform distribution of X
n = 700
set.seed(jobid)
X_4 = rnorm(n, mean = 0, sd = 0.5)
m_4 = X_4 + 2*exp(-16*X_4^2)
eps4 = rnorm(n, mean = 0, sd = 0.2)
Y_4 = m_4 + eps4

m4SecDeriv = function(x){
  res = 64*(32*x^2-1)*exp(-16*x^2)
  return(res)
}

## Proposed estimator
U_4 = kcde(X_4, eval.points = X_4, tail.flag = "lower.tail")$estimate
k_opt4 = OptK(Y = Y_4, U = U_4, k_range = seq(1, 300, by = 1), h = NULL, kern = gaussK)
print(k_opt4)
deri_est4 = DerivEstQuotient(Y = Y_4, X = X_4, k = k_opt4)

k_opt_2d4 = OptK1K2(Y = Y_4, U = U_4, k1_range = seq(1, 100, by = 1), 
                    k2_range = seq(1, 100, by = 1), h = NULL, kern = gaussK)
# k_opt_2d4 = list("k1" = 15, "k2" = 8)
print(k_opt_2d4)
deri_sec4 = SecondDerivEstQuotient(Y = Y_4, X = X_4, k1 = k_opt_2d4$k1, 
                                   k2 = k_opt_2d4$k2)

h_opt4 = RSSBandwidth(Y = deri_est4$Y_deri, U = deri_est4$U_ord, k = k_opt4, 
                      h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt4)
h_opt4 = h_opt4*1.01431
h_opt_2d4 = RSSBandwidth(Y = deri_sec4$Y_sec, U = deri_sec4$U_ord, 
                         k = k_opt_2d4$k1 + k_opt_2d4$k2, 
                         h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt_2d4)
h_opt_2d4 = h_opt_2d4*1.01431

qry_range = seq(k_opt_2d4$k1 + k_opt_2d4$k2 + 1, 700 - k_opt_2d4$k1 - k_opt_2d4$k2, by = 1)
Y_deriv_sm4 = LocalPolyReg(Y = deri_est4$Y_deri[qry_range],
                           X = deri_est4$U_ord[qry_range], 
                           xeval = deri_est4$U_ord,
                           h = h_opt4, kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_sm4 = LocalPolyReg(Y = deri_sec4$Y_sec[qry_range], X = deri_sec4$U_ord[qry_range],
                         xeval = deri_sec4$U_ord, h = h_opt_2d4,
                         kern = gaussK, deg = 3, deriv_ord = 0)
# Y_deriv_sm4 = LocalPolyReg(Y = deri_est4$Y_deri,
#                            X = deri_est4$U_ord, 
#                            xeval = deri_est4$U_ord,
#                            h = h_opt4, kern = gaussK, deg = 3, deriv_ord = 0)
# Y_sec_sm4 = LocalPolyReg(Y = deri_sec4$Y_sec, X = deri_sec4$U_ord,
#                          xeval = deri_sec4$U_ord, h = h_opt_2d4,
#                          kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_trans4 = Y_sec_sm4 * deri_sec4$den_ord^2 + Y_deriv_sm4 * deri_sec4$den_deriv_ord

MAE_quo = mean(abs(Y_sec_trans4 - m4SecDeriv(deri_sec4$X_ord))[31:670])
Med_Abs_Err = median(abs(Y_sec_trans4 - m4SecDeriv(deri_sec4$X_ord)))


## Proposed estimator with the oracle data-generating distribution
U_4o = pnorm(X_4, mean = 0, sd = 0.5)
k_opt4o = OptK(Y = Y_4, U = U_4o, k_range = NULL, h = NULL, kern = gaussK)
print(k_opt4o)
deri_est4o = DerivEstQuotient(Y = Y_4, X = X_4, k = k_opt4o, U = U_4o)
k_opt_2d4o = OptK1K2(Y = Y_4, U = U_4o, k1_range = seq(1, 100, by = 1), 
                     k2_range = seq(1, 100, by = 1), h = NULL, kern = gaussK)
# k_opt_2d4o = list("k1" = 15, "k2" = 8)
print(k_opt_2d4o)
deri_sec4o = SecondDerivEstQuotient(Y = Y_4, X = X_4, k1 = k_opt_2d4o$k1, 
                                    k2 = k_opt_2d4o$k2, U = U_4o)

h_opt4o = RSSBandwidth(Y = deri_est4o$Y_deri, U = deri_est4o$U_ord, k = k_opt4o, 
                       h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt4o)
h_opt4o = h_opt4o*1.01431
h_opt_2d4o = RSSBandwidth(Y = deri_sec4o$Y_sec, U = deri_sec4o$U_ord, 
                          k = k_opt_2d4o$k1 + k_opt_2d4o$k2, 
                          h_range = seq(0.03, 0.1, by = 0.005), deg = 3)
print(h_opt_2d4o)
h_opt_2d4o = h_opt_2d4o*1.01431

qry_range_o = seq(k_opt_2d4o$k1 + k_opt_2d4o$k2 + 1, 700 - k_opt_2d4o$k1 - k_opt_2d4o$k2, by = 1)
# Y_deriv_sm4o = LocalPolyReg(Y = deri_est4o$Y_deri,
#                             X = deri_est4o$U_ord,
#                             xeval = deri_est4o$U_ord, h = h_opt4o,
#                             kern = gaussK, deg = 3, deriv_ord = 0)
# Y_sec_sm4o = LocalPolyReg(Y = deri_sec4o$Y_sec,
#                           X = deri_sec4o$U_ord, 
#                           xeval = deri_est4o$U_ord,
#                           h = h_opt_2d4o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_deriv_sm4o = LocalPolyReg(Y = deri_est4o$Y_deri[qry_range_o],
                            X = deri_est4o$U_ord[qry_range_o],
                            xeval = deri_est4o$U_ord, h = h_opt4o,
                            kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_sm4o = LocalPolyReg(Y = deri_sec4o$Y_sec[qry_range_o],
                          X = deri_sec4o$U_ord[qry_range_o], 
                          xeval = deri_est4o$U_ord,
                          h = h_opt_2d4o, kern = gaussK, deg = 3, deriv_ord = 0)
Y_sec_trans4o = Y_sec_sm4o * dnorm(deri_est4o$X_ord, mean = 0, sd = 0.5)^2 - Y_deriv_sm4o * (deri_est4o$X_ord/(sqrt(2*pi)*0.5^3))*exp(-deri_est4o$X_ord^2/(2*0.5^2))

MAE_quo_ora = mean(abs(Y_sec_trans4o - m4SecDeriv(deri_sec4o$X_ord))[31:670])
Med_Abs_Err_ora = median(abs(Y_sec_trans4o - m4SecDeriv(deri_est4o$X_ord)))

deri_sec_res4 = data.frame(jobid = rep(jobid, 2), k_opt = c(k_opt4, k_opt4o), 
                           k1_opt = c(k_opt_2d4$k1, k_opt_2d4o$k1),
                           k2_opt = c(k_opt_2d4$k2, k_opt_2d4o$k2),
                           MAE_adj = c(MAE_quo, MAE_quo_ora), 
                           Med_Abs_Err = c(Med_Abs_Err, Med_Abs_Err_ora), 
                           scenario = c("Proposed", "Oracle"))

write.csv(deri_sec_res4, paste0("./deri_res/deri_second_order_", jobid, "_sim4_new.csv"), 
          row.names=FALSE)
