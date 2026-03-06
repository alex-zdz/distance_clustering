library(AntMAN)
library(salso)
source("salso_functions.R")

set.seed(123)
n=50
p_true=0.3
mean1_true=1
sd1_true=1
mean2_true=5
sd2_true=1/2
observed_data = rep(NA,n)
for (i in 1:n) {
  U = runif(n=1,0,1)
  if(U<p_true){
    observed_data[i] = rnorm(n=1, mean=mean1_true, sd = sd1_true)
  }else{
    observed_data[i] = rnorm(n=1, mean=mean2_true, sd = sd2_true)
  }
}

mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=1, k0=0.3333333, nu0=4.222222, sig02=3.661027)
mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AntMAN::AM_mcmc_fit(
  y = observed_data,
  mix_kernel_hyperparams = mixture_uvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
# plot (fit)

eam = AM_clustering(fit)  # M*n matrix
AMcluster = AM_salso(eam, "binder")

start_time <- Sys.time()
mysalso <- Salso_whole_procedure(eam, n_runs=5)
end_time <- Sys.time()
print(end_time - start_time)

true_salso <- salso(eam,loss="binder")

print(true_salso)
print(AMcluster)
print(mysalso)
print(sum(AMcluster==mysalso))
print(sum(true_salso==mysalso))
# library(doParallel)
# registerDoParallel(cores=4)
# library(foreach)
# 
# system.time(Salso_parallel(eam, n_runs=5,Salso_ith_run))
