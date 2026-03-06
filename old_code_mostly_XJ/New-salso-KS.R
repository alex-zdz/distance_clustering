# New salso with W2
# example of mixture model simulation data


rm(list = ls())

library(AntMAN)
library(salso)
library(invgamma)
library(fossil)
source("KS-functions.R")
source("Salso-functions-generalformat.R")
set.seed(143)

# Simulate data

# data generation
N = 100
c_alloc_true = sample(1:5, N, replace = TRUE)
hist(c_alloc_true)
K=5
y = rep(NA, N)
# m = c(20,10,0,5,-2)
# m = c(-2,0,3)
m = c(7,10,0,5,-2)
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
prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)


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


D = KS_distance_mixture(y, eam[1,], prior_list)

res_ini = Initialization_phase(eam, y, prior_list,KS_distance_mixture)
res_swt = Sweetening_phase(res_ini, y, prior_list,KS_distance_mixture)
res_zea = Zealous_updates_phase(res_swt, y, prior_list,KS_distance_mixture, n_max=10, k=10)
res_salso = Salso_whole_procedure(eam, y,prior_list,KS_distance_mixture, n_runs=10, n_max=10 )


KS_distance_mixture(y, res_ini, prior_list)
KS_distance_mixture(y, res_swt, prior_list)
KS_distance_mixture(y, res_zea, prior_list)
KS_distance_mixture(y, res_salso, prior_list)



print("RI")
print(rand.index(c_alloc_true,AMcluster))
print(rand.index(c_alloc_true,AMcluster2))
print(rand.index(c_alloc_true,res_ini))
print(rand.index(c_alloc_true,res_swt))
print(rand.index(c_alloc_true,res_zea))
print(rand.index(c_alloc_true,res_salso))

print("ARI")
print(adj.rand.index(c_alloc_true,AMcluster))
print(adj.rand.index(c_alloc_true,AMcluster2))
print(adj.rand.index(c_alloc_true,res_ini))
print(adj.rand.index(c_alloc_true,res_swt))
print(adj.rand.index(c_alloc_true,res_zea))
print(adj.rand.index(c_alloc_true,res_salso))


# [1] "RI"
# [1] 0.7024242
# [1] 0.7024242
# [1] 0.8444444
# [1] 0.8915152
# [1] 0.8915152
# [1] 0.8874747
# [1] "ARI"
# [1] 0.4035498
# [1] 0.4035498
# [1] 0.5721194
# [1] 0.7707583
# [1] 0.7707583
# [1] 0.7733318