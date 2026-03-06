rm(list = ls())

library(AntMAN)
library(salso)
library(invgamma)
library(fossil)
source("salso_functions.R")
source("W2-functions.R")

# set.seed(123)
# n=50
# p_true=0.3
# mean1_true=1
# sd1_true=1
# mean2_true=5
# sd2_true=1/2
# observed_data = rep(NA,n)
# for (i in 1:n) {
#   U = runif(n=1,0,1)
#   if(U<p_true){
#     observed_data[i] = rnorm(n=1, mean=mean1_true, sd = sd1_true)
#   }else{
#     observed_data[i] = rnorm(n=1, mean=mean2_true, sd = sd2_true)
#   }
# }

set.seed(123)

# data generation
N = 300
c_alloc_true = sample(1:3, N, replace = TRUE)
hist(c_alloc_true)
K=3
y = rep(NA, N)
m = c(-2,0,3)
s = c(0.4, 1, 1.1)
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


mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=theta, k0=kappa0, nu0=a, sig02=b)
mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=5, verbose=1)
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
hist(AMcluster)


# start_time <- Sys.time()
# mysalso <- Salso_whole_procedure(eam, n_runs=5)
# end_time <- Sys.time()
# print(end_time - start_time)

true_salso <- salso(eam,loss="binder")

print(true_salso)
print(AMcluster)
print(sum(true_salso==AMcluster))
# print(mysalso)
# print(sum(AMcluster==mysalso))
# print(sum(true_salso==mysalso))
# library(doParallel)
# registerDoParallel(cores=4)
# library(foreach)
# 
# system.time(Salso_parallel(eam, n_runs=5,Salso_ith_run))


prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)

# test W2
c_alloc = AMcluster

# W2(data=y, c_alloc, prior_list)

res = min_W2(eam, y, prior_list) 
print(res)
print("How many cluster labels are the same?")
print(sum(true_salso==res))

# which((true_salso - res)!=1)

print(rand.index(c_alloc_true,res)) # 0.8033445
print(rand.index(c_alloc_true,true_salso)) # 0.8033445
