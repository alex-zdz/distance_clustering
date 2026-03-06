# 2d example

rm(list = ls())

library(MASS)
library(matlib)
library(AntMAN)
# library(salso)

source("Pearson-functions-2d.R")
source("KS-functions-2d.R")
source("Salso-functions-generalformat.R")

set.seed(123)


# data generation
N = 30
d = 2
c_alloc_true = sample(1:3, N, replace = TRUE)
hist(c_alloc_true)
K=3
y = matrix(NA, N, d)
mu = c(-2,0,3)
s = c(0.4, 1, 1.1)
for (i in 1:K) {
  I = c_alloc_true==i
  y[I,] = mvrnorm(sum(I), rep(mu[i],d), s[i]*diag(d))
}
hist(y[,1])
hist(y[,2])


# prior
alpha = 1
kappa0 = 0.01
m0 = rep(0,d)
nu0 = d+2
Lambda0 = diag(d)*10
prior_list = list("alpha"=alpha, "kappa0"=kappa0, "m0"=m0, "nu0"=nu0, "Lambda0"=Lambda0)



c_alloc = c_alloc_true
# posterior
alpha_post=numeric(K)
cluster=sort(unique(c_alloc))
# alpha_post
 ##mean of Dir
###parameters for kernel
# kappa_post = rep(NA,3)
mu_post = matrix(NA, K,d)
# sigma2_post = rep(NA,3)
Lambda_post = vector("list",K)
Sigma_post = vector("list",K)

Sk_calculation <- function(yk_mat){
  d = dim(yk_mat)[2]
  Sk = matrix(0, d,d)
  col_mean = matrix(colMeans(yk_mat), nrow = d, byrow = TRUE)
  for (i in 1:dim(yk_mat)[1]) {
    yk_mat_rowi = matrix(yk_mat[i,], nrow = d, byrow = TRUE)
    Sk = Sk + (yk_mat_rowi-col_mean) %*% t(yk_mat_rowi-col_mean)
  }
  return(Sk)
}



#     
for(i in 1:K){
  I = c_alloc==cluster[i]
  nk=sum(I)
  yk=y[I,]
  Sk = Sk_calculation(yk)  # [d,d]
  alpha_post[i]=alpha+sum(c_alloc==cluster[i])
  
  kappa_post_k = kappa0 +nk
  nu_post_k = nu0 + nk
  mu_post[i,] = (kappa0 * m0 + colSums(yk))/kappa_post_k
  Lambda_post_k = Lambda0 + Sk + t(t(colMeans(yk)-m0)) %*% t((colMeans(yk)-m0)) * (kappa0*nk/kappa_post_k)
  Lambda_post[[i]] = Lambda_post_k
  Sigma_post[[i]] = (Lambda_post_k)/(nu_post_k - d -1)
}

weights=alpha_post/sum(alpha_post)

# 
# library(mltools)
# library(data.table)
# 
# 
# dt <- data.table(x = y[,1], y = y[,2])
# ecdf_2d = empirical_cdf(dt, ubounds = data.table(x = -50, y = -50))


# Pearson_distance_mixture_2d(observed_data=y, c_alloc_true, prior_list)

# run ANTMAN to get some MCMC samples
mixture_mvn_params = AntMAN::AM_mix_hyperparams_multinorm(mu0=m0, ka0=kappa0, nu0=nu0, Lam0=Lambda0)
mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AntMAN::AM_mcmc_fit(
  y = y,
  mix_kernel_hyperparams = mixture_mvn_params,
  mix_components_prior =components_prior,
  # mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
# plot (fit)

eam = AM_clustering(fit)  # M*n matrix 
AMcluster = AM_salso(eam, "binder")
AMcluster2 = AM_salso(eam, "VI")
hist(AMcluster)

KS_distance_mixture_2d(y, c_alloc, prior_list)
KS_distance_mixture_2d(y, eam[2,], prior_list)
KS_distance_mixture_2d(y, sample(1:3, N, replace = TRUE), prior_list)

Pearson_distance_mixture_2d(y, eam[2,], prior_list, bins=seq(-10,10,by=0.5))
Pearson_distance_mixture_2d(y, c_alloc, prior_list, bins=seq(-10,10,by=0.5))

start_time <- Sys.time()
c_first = Initialization_phase(eam[1:10,], y, prior_list,KS_distance_mixture_2d)
# result = Salso_whole_procedure(eam, y,prior_list, Pearson_distance_mixture_2d, n_runs=1, n_max=3)
end_time <- Sys.time()
print(end_time - start_time)
#   # test posterior 1d and 2d
# alpha = 1
# kappa0 = 0.01
# theta = 0
# a = 2
# b = 5
# prior_list_1d = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)
# 
# alpha = 1
# kappa0 = 0.01
# m0 = rep(0,d)
# nu0 = d+2
# Lambda0 = diag(d) * 10
# prior_list_2d = list("alpha"=alpha, "kappa0"=kappa0, "m0"=m0, "nu0"=nu0, "Lambda0"=Lambda0)
# 
# posterior_1d = mixture_posterior(c_alloc, y[,1], prior_list_1d)
# posterior_list_test = mixture_posterior_2d(c_alloc, observed_data, prior_list_2d)
# 
# Pearson_distance_mixture_2d(y, c_alloc, prior_list_2d, bins=seq(-10,10,by=0.5))
# 
# Pearson_distance_mixture(y[,1], eam[1,], prior_list_1d)
