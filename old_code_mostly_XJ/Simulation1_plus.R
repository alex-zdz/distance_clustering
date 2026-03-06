# this file contains the original simulations by xj as well as an additional run with the functions contained in "salso_functions"
# which seem more accurate w.r.t. the salso paper
# Next step is to use the more accurate functions in our algorithm


rm(list = ls())

library(AntMAN)
library(salso)
library(invgamma)
library(fossil)
library(sn)
library(Rcpp)
library(ggplot2)
library(patchwork)
library(tidyr)
sourceCpp("W2-functions.cpp")
source("Pearson-functions.R")
source("KS-functions-1d.R")
source("Salso-functions-generalformat.R")
source("MCMC-fit-plot.R")
set.seed(143)

# Load the new functions into a separate environment to avoid naming issues
env_new <- new.env()
sys.source("test files/salso_functions.R", envir = env_new)

data_generation <- function(N, type="Gaussian mixture"){
  
  if(type=="Gaussian mixture"){
    # data generation
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    # hist(c_alloc_true)
    
    y = rep(NA, N)
    # m = c(20,10,0,5,-2)
    m = c(-2,0,3)
    # m = c(7,10,0,5,-2)
    s = c(0.4, 1, 0.3)
    for (i in 1:K) {
      I = c_alloc_true==i
      y[I] = rnorm(sum(I), m[i], s[i])
    }
    # hist(y)
  }else if(type=="Skew Gaussian mixture"){
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    y = rep(NA, N)
    m = c(-2,0,3)
    s = c(0.4, 1, 0.3)
    a = c(1,10,4)
    
    for (i in 1:K) {
      I = c_alloc_true==i
      y[I] = rsn(sum(I),  m[i], s[i], a[i])
    }
    # plot(density(y))
  }else if(type=="Skew-symmetric mixture"){
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    y = rep(NA, N)
    
    I2 = c_alloc_true==2
    y[I2] = rsn(sum(I2),  0, 5, 4)
    I3 = c_alloc_true==3
    y[I3] = rnorm(sum(I3),  -4, sqrt(0.5))
    
    I1 = c_alloc_true==1
    set1 = which(I1)
    w11 = c(0.364,0.212,0.424)
    c11 = sample(1:3, sum(I1), replace = TRUE, prob=w11)
    I11 = c11==1
    y[set1[I11]] = rsn(sum(I11),  2.5, 1, -10)
    I12 = c11==2
    y[set1[I12]] = rnorm(sum(I12),  2.325, sqrt(0.2))
    I13 = c11==3
    y[set1[I13]] = rnorm(sum(I13),  1.085, sqrt(0.7))
    
    # plot(density(y))
  }

  
  return(list(data=y, c_true = c_alloc_true))
}



N = 100

# prior
alpha = 1
kappa0 = 0.01
theta = 0
a = 2
b = 5
prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)

# 22.01.2026 Alex: testing XJs methods with own data generating function (no hardcoded means)
generate_mixture_data <- function(
    N,
    K = 3,
    dim = 1,
    alpha = 1,
    mu_true = NULL,      # vector for dim = 1, matrix (K x dim) for dim > 1
    Sigma_true = NULL         # vector for dim = 1, list of covariance matrices for dim > 1
) {
  #--- Component weights from Dirichlet(alpha,...,alpha) ---
  gamma_draws <- rgamma(K, shape = alpha, rate = 1)
  weights_true <- gamma_draws / sum(gamma_draws)
  
  #--- Assign clusters ---
  cluster_true <- sample(1:K, N, replace = TRUE, prob = weights_true)
  
  #--- Default mu_true ---
  if (dim == 1) {
    if (is.null(mu_true)) mu_true <- seq(-K, K, length.out = K)
  } else {
    if (is.null(mu_true)) {
      mu_true <- matrix(0, nrow = K, ncol = dim)
      for (i in 1:K) mu_true[i, ] <- rep(i, dim)  # simple separated defaults, consider making random
    }
  }
  
  #--- Default covariances ---
  if (dim == 1) {
    if (is.null(Sigma_true)) Sigma_true <- rep(1, K)
  } else {
    if (is.null(Sigma_true)) {
      Sigma_true <- vector("list", K)
      for (i in 1:K) Sigma_true[[i]] <- diag(dim) * 0.1     # identity covariance
    }
  }
  
  #--- Initialize data ---
  y <- matrix(NA, nrow = N, ncol = dim)
  
  #--- Generate Gaussian data ---
  for (i in 1:K) {
    idx <- cluster_true == i
    ni  <- sum(idx)
    
    if (ni > 0) {
      if (dim == 1) {
        y[idx, ] <- rnorm(ni, mean = mu_true[i], sd = Sigma_true[i])
      } else {
        y[idx, ] <- MASS::mvrnorm(n = ni, mu = mu_true[i, ], Sigma = Sigma_true[[i]])
      }
    }
  }
  
  #--- Return vector for 1d case ---
  if (dim == 1) y <- as.numeric(y)
  
  return(list(
    data = y,
    cluster_true = cluster_true,
    weights_true = weights_true,
    mu_true = mu_true,
    Sigma_true = Sigma_true
  ))
}


n_replication=10
#new test
n_replication=2

ARI_binder = numeric(n_replication)
ARI_VI = numeric(n_replication)
RI_binder = numeric(n_replication)
RI_VI = numeric(n_replication)
RI_distance_Pearson = numeric(n_replication)
ARI_distance_Pearson = numeric(n_replication)
RI_distance_W2 = numeric(n_replication)
ARI_distance_W2 = numeric(n_replication)
RI_distance_KS = numeric(n_replication)
ARI_distance_KS = numeric(n_replication)

#pdf("Boxplot-GaussianMixture-N100.pdf")
start_time <- Sys.time()

for (i_replication in 1:n_replication){
  print("i_replication=")
  print(i_replication)
  
  # sample_i = data_generation(N, type="Skew-symmetric mixture")
  sample_i = data_generation(N, type="Gaussian mixture")
  # sample_i = data_generation(N, type="Skew Gaussian mixture")
  
  data_i = sample_i$data
  c_true_i = sample_i$c_true
  
  # 22.01.2026 different data generation:
  N <- 100
  K <- 3
  D = 2
  alpha_0 <- 1
  
  simdata <- generate_mixture_data(
    N,
    K,
    dim = D,
    alpha = alpha_0,
    mu_true = NULL,      
    Sigma_true = NULL 
  )
  
  data_i <- simdata$data
  c_true_i <- simdata$cluster_true
  
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
  AMcluster_binder = AM_salso(eam, "binder")
  
  # binder_salso = salso(eam, loss = binder(a = 0.1))
  # case_match(c(binder_salso), 1~1, 2 ~ 3, 3 ~2)
  
  AMcluster_VI = AM_salso(eam, "VI")
  ARI_binder[i_replication] = adj.rand.index(c_true_i,AMcluster_binder)
  ARI_VI[i_replication] = adj.rand.index(c_true_i,AMcluster_VI)
  RI_binder[i_replication] = rand.index(c_true_i,AMcluster_binder)
  RI_VI[i_replication] = rand.index(c_true_i,AMcluster_VI)
  
  print("Start Pearson")
  res_Pearson = Salso_whole_procedure(eam, data_i,prior_list, Pearson_distance_mixture, n_runs=3, n_max=3)
  print("Start W2")
  res_W2 = Salso_whole_procedure(eam, data_i, prior_list, W2, n_runs=3, n_max=3)
  print("Start KS")
  res_KS = Salso_whole_procedure(eam, data_i, prior_list, KS_distance_mixture, n_runs=3, n_max=3)
  
  # Binder estimate using the salso_functions
  binder_2 <- env_new$Salso_whole_procedure(eam, 2)
  
  RI_distance_Pearson[i_replication] = rand.index(c_true_i,res_Pearson)
  ARI_distance_Pearson[i_replication] = adj.rand.index(c_true_i,res_Pearson)
  RI_distance_W2[i_replication] = rand.index(c_true_i,res_W2)
  ARI_distance_W2[i_replication] = adj.rand.index(c_true_i,res_W2)
  RI_distance_KS[i_replication] = rand.index(c_true_i,res_KS)
  ARI_distance_KS[i_replication] = adj.rand.index(c_true_i,res_KS)
  
  # fit_binder = MCMC_fit(n_iter=5000, AMcluster_binder,data_i, prior_list, len_range=5000) 
  # fit_VI = MCMC_fit(n_iter=5000, AMcluster_VI,data_i, prior_list, len_range=5000) 
  # fit_W2 = MCMC_fit(n_iter=5000, res_W2,data_i, prior_list, len_range=5000) 
  # fit_Pearson = MCMC_fit(n_iter=5000, res_Pearson,data_i, prior_list, len_range=5000) 
  # fit_KS = MCMC_fit(n_iter=5000, res_KS,data_i, prior_list, len_range=5000) 
  # 
  # x_range = x_range = seq(-10,10,length.out=5000)
  # hist(data_i, prob=TRUE, ylim=c(0,0.35))
  # lines(x_range, fit_binder, col="#F8766D",lwd=3)
  # lines(x_range, fit_VI, col="#A3A500",lwd=3)
  # lines(x_range, fit_Pearson, col="#00BF7D",lwd=3)
  # lines(x_range, fit_W2, col="#00B0F6",lwd=3)
  # lines(x_range, fit_KS, col="#E76BF3",lwd=3)
  # legend("topright", c("Binder", "VI", "Pearson","W2", "KS"), col=c("#F8766D","#A3A500", "#00BF7D","#00B0F6","#E76BF3"),lty=1, lwd=3, cex=1)
  
  # ggplot(home_data, aes(x = price, y = after_stat(density))) +
  #   geom_histogram() +
  #   geom_vline(aes(xintercept = mean_price), price_stats, color = "red", linewidth = 2) +
  #   geom_density(color = "green", linewidth = 2)
  
}
end_time <- Sys.time()
print(end_time - start_time)

ARI100 = c(ARI_binder,ARI_VI,ARI_distance_Pearson, ARI_distance_W2,ARI_distance_KS)
RI100 = c(RI_binder,RI_VI,RI_distance_Pearson, RI_distance_W2,RI_distance_KS)
name = c(rep("Binder",n_replication), rep("VI",n_replication), rep("Pearson",n_replication),rep("W2",n_replication),rep("KS",n_replication))

# name=c("binder","VI","Pearson",“W2”)
# value=ARI

library(dplyr)
library(ggplot2)

# pdf("Boxplot-SkewSymmetricMixture-N100.pdf")
# pdf("Boxplot-SkewGaussianMixture-N100.pdf")


data = data.frame("name"=name, "value"=ARI100)
data$name <- factor(data$name, levels=c("Binder", "VI", "Pearson","W2","KS"))
p1 = data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("ARI boxplot") +
  xlab("")

data = data.frame("name"=name, "value"=RI100)
data$name <- factor(data$name, levels=c("Binder", "VI", "Pearson","W2","KS"))
p2 = data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("RI boxplot") +
  xlab("")


p1|p2

dev.off()



library(scales)
show_col(hue_pal()(5))
