rm(list = ls())

library(AntMAN)
library(salso)
library(invgamma)
library(fossil)
library(sn)
# old functions by XJ
#source("C:/Users/alexm/NUS Dropbox/Alexander Mozdzen/distance_clustering/code/KS-functions-1d.R")
source("C:/Users/alexm/NUS Dropbox/Alexander Mozdzen/distance_clustering/code/Salso-functions-generalformat.R")
# new function
source("main_raw.R")
source("KS-functions-1d_raw.R")

#source("MCMC-fit-plot.R")

# Load the new functions into a separate environment to avoid naming issues
#ax_env <- new.env()
#sys.source("KS-functions-1d_raw.R", envir = ax_env)

N = 100

# prior
alpha = 1
kappa0 = 0.01
theta = 0
a = 2
b = 5
prior_list = list("alpha" = alpha, "kappa0" = kappa0, "theta" = theta, "a" = a, "b" = b)

# Simulate data  
sample_i = data_generation(N, type="Gaussian mixture")
data_i = sample_i$data
c_true_i = sample_i$c_true

# run ANTMAN to get some MCMC samples
mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm(m0=theta, k0=kappa0, nu0=a, sig02=b)
mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)
fit <- AntMAN::AM_mcmc_fit(y = data_i,
                           mix_kernel_hyperparams = mixture_uvn_params,
                           mix_components_prior =components_prior,
                           # mix_weight_prior = weights_prior,
                           mcmc_parameters = mcmc_params)
eam = AM_clustering(fit)  # M*n matrix
# Compute Binder and VI estimates:
AMcluster_binder = AM_salso(eam, "binder")
AMcluster_VI = AM_salso(eam, "VI")
# ARI_binder[i_replication] = adj.rand.index(c_true_i,AMcluster_binder)
# ARI_VI[i_replication] = adj.rand.index(c_true_i,AMcluster_VI)
# RI_binder[i_replication] = rand.index(c_true_i,AMcluster_binder)
# RI_VI[i_replication] = rand.index(c_true_i,AMcluster_VI)

res_KS_xj = Salso_whole_procedure(eam, data_i, prior_list, KS_distance_mixture, n_runs=3, n_max=3)
# since AM output starts at 0
c_sample <- apply(eam, 2, function(x) x + 1)

#environment(mde_main) <- ax_env

res_KS_ax = mde_main(c_sample = c_sample , data = data_i, prior_list,
                     distance_function = KS_distance_mixture,
                     n_runs = 3, n_max = 4, n_split = 4)

# Speed comparison
#library(microbenchmark)
# benchmark_results <- microbenchmark(
#   res_KS_xj = Salso_whole_procedure(eam, data_i, prior_list, KS_distance_mixture, n_runs = 3, n_max = 3),
#   res_KS_ax = mde_main(c_sample = c_sample, data = data_i, prior_list,
#                        distance_function = KS_distance_mixture, n_runs = 3, n_max = 3, n_split = 3),
#   times = 2L # adjust the number of repetitions if needed
# )
# print(benchmark_results)


adj.rand.index(c_true_i,AMcluster_binder)

# results seed dependent!
adj.rand.index(c_true_i,res_KS_xj)

adj.rand.index(c_true_i,res_KS_ax)





