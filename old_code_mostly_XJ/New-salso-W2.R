# New salso with W2
# example of mixture model simulation data


rm(list = ls())

library(AntMAN)
library(salso)
library(invgamma)
library(fossil)
library(Rcpp)
# source("W2-functions.R")
sourceCpp("W2-functions.cpp")
source("Salso-functions-generalformat.R")
set.seed(143)

# Simulate data

# data generation
N = 1000
c_alloc_true = sample(1:5, N, replace = TRUE)
hist(c_alloc_true)
K=5
y = rep(NA, N)
m = c(7,10,0,5,-2)
# m = c(-2,0,3)
s = c(0.4, 1, 0.3,1.1, 1)
for (i in 1:K) {
  I = c_alloc_true==i
  y[I] = rnorm(sum(I), m[i], s[i])
}
hist(y)


# prior
alpha = 1
kappa0 = 0.01
theta = 0
a = 2
b = 5


# s2 = rinvgamma(3, shape=a, scale = b)
# s = sqrt(s2)
# m = rep(0,3)
# m[1] = rnorm(n=1, mean=theta, sd=s[1]/sqrt(kappa0))
# m[2] = rnorm(n=1, mean=theta, sd=s[2]/sqrt(kappa0))
# m[3] = rnorm(n=1, mean=theta, sd=s[3]/sqrt(kappa0))
# for (i in 1:K) {
#   I = c_alloc_true==i
#   y[I] = rnorm(sum(I), m[i], s[i])
# }
# hist(y)


# run ANTMAN to get some MCMC samples
mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=theta, k0=kappa0, nu0=a, sig02=b)
mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AntMAN::AM_mcmc_fit(
  y = y,
  mix_kernel_hyperparams = mixture_uvn_params,
  mix_components_prior =components_prior,
  # mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
# plot (fit)

eam = AM_clustering(fit)  # M*n matrix
AMcluster = AM_salso(eam, "binder")
AMcluster2 = AM_salso(eam, "VI")
hist(AMcluster)

# test 
prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)
res_ini = Initialization_phase(eam, y, prior_list,W2)
res_swt = Sweetening_phase(res_ini, y, prior_list,W2)
res_zea = Zealous_updates_phase(res_swt, y, prior_list,W2, n_max=30, k=10)

# hist(c_alloc_true)
# hist(res_zea)
# hist(AMcluster)
# W2(y, c_alloc_true, prior_list)  #[1] 0.05082871


print("RI")
print(rand.index(c_alloc_true,AMcluster))
print(rand.index(c_alloc_true,AMcluster2))
print(rand.index(c_alloc_true,res_ini))
print(rand.index(c_alloc_true,res_swt))
print(rand.index(c_alloc_true,res_zea))
print("ARI")
print(adj.rand.index(c_alloc_true,AMcluster))
print(adj.rand.index(c_alloc_true,AMcluster2))
print(adj.rand.index(c_alloc_true,res_ini))
print(adj.rand.index(c_alloc_true,res_swt))
print(adj.rand.index(c_alloc_true,res_zea))

# [1] "RI"
# [1] 0.9214716
# [1] 0.9186845
# [1] 0.9271795
# [1] 0.9321516
# [1] 0.9321516
# [1] "ARI"
# [1] 0.7664984
# [1] 0.7566688
# [1] 0.7797114
# [1] 0.8500205
# [1] 0.8500205

# 
# start_time <- Sys.time()
# res_salso = Salso_whole_procedure(eam, y,prior_list, n_runs=3, n_maxzea=0)
# end_time <- Sys.time()
# print(end_time - start_time)

# # test W2
# sourceCpp("W2-functions.cpp")
# source("W2-functions.R")
# W2(y, eam[5,], prior_list)  # -0.7717862 -- R   [1] -0.7717906 -- Rcpp
# res_ini = Initialization_phase(eam, y, prior_list)
# print(sum(AMcluster==res_ini))  # n_maxzea=0 all equal
# print(rand.index(c_alloc_true,res_ini))  # 0.8
# 
# 
# F_inv( 300/300,c(0.4, 0.3, 0.3),m, s) #-11.91708-R-Rcpp
# x = seq(from = -20, to = 100,length.out=500)
# y = numeric(length(x))
# for (i in 1:length(x)) {
#   y[i] = F(x[i], c(0.4, 0.3, 0.3), m, s)
# }
# 
# plot(x,y)


# # print(sum(AMcluster==res_swt))
# print(sum(AMcluster==res_salso))  # n_maxzea=0 all equal
# 
# # print(rand.index(c_alloc_true,res_swt))  # 0.8
# print(rand.index(c_alloc_true,AMcluster))  # 0.8
