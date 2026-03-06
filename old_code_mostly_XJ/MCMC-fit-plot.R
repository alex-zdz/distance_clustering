# MCMC with fixed c

library(LaplacesDemon)

calculate_n0 <- function(c_vec_0,M){
  n0_vec = rep(0,M)
  for(k in 1:M){
    n0_vec[k] = length(which(c_vec_0==k))
  }
  return(n0_vec)
}


omega_update <- function(c, alpha){
  a_vec = calculate_n0(c,max(c)) + rep(alpha, max(c))
  omega_new = rdirichlet(n=1, a_vec) 
  return(omega_new)
}


mu_update <- function(mu_unique, sigmasq_unique, Y, c, kappa0){
  K0 = max(c)
  mu_unique_new = rep(NA,K0)
  
  for(k in 1:K0){
    set_k = which(c==k)
    n_k = length(set_k)
    A = sum(Y[set_k])/sigmasq_unique[k] + kappa0/sigmasq_unique[k]
    C = n_k/sigmasq_unique[k] + kappa0/sigmasq_unique[k]
    
    var_star = 1/C
    mean_star = A/C
    mu_unique_new[k] = rnorm(1, mean_star, sqrt(var_star))
  }
  return(mu_unique_new)
}

sigmasq_update <- function(sigmasq_unique, mu_unique, Y, c, a, b, theta, kappa0){
  K0 = max(c)
  sigmasq_unique_new = rep(NA, K0)
  
  for(k in 1:K0){
    set_k = which(c==k)
    n_k = length(set_k)
    alpha_star = n_k/2 + 1/2 + a
    beta_star = 0.5*sum((Y[set_k]-mu_unique[k])^2)+b + 0.5*kappa0*(mu_unique[k]-theta)
    sigmasq_unique_new[k] = invgamma::rinvgamma(n=1, shape=alpha_star, rate=beta_star)
  }
  return(sigmasq_unique_new)
}


MCMC_fit <- function(n_iter, c_fixed, Y, prior_list, len_range){
  
  K0 = max(c_fixed)
  kappa0 = prior_list$kappa0
  theta = prior_list$theta
  alpha = prior_list$alpha
  a = prior_list$a
  b = prior_list$b
  
  mu_mcmc = matrix(NA, nrow=n_iter, ncol=K0)
  sigmasq_mcmc = matrix(NA, nrow=n_iter, ncol=K0)
  omega_mcmc = matrix(NA, nrow=n_iter, ncol=K0)
  
  mu_unique = rep(1, K0)
  sigmasq_unique = rep(1, K0)
  
  f_Gm = function(x,w,u,s) sum( w*dnorm(x,mean=u,sd=s) )
  x_range = seq(-10,10,length.out=len_range)
  fit_all = rep(0, len_range)
    
  for (i in 1:n_iter) {
    mu_unique = mu_update(mu_unique, sigmasq_unique, Y, c_fixed, kappa0)
    sigmasq_unique = sigmasq_update(sigmasq_unique, mu_unique, Y, c_fixed, a, b, theta, kappa0)
    omega = omega_update(c_fixed, alpha)
    
    mu_mcmc[i,] = mu_unique
    sigmasq_mcmc[i,] = sigmasq_unique
    omega_mcmc[i,] = omega
    
    fit_i = lapply(x_range, f_Gm, omega, mu_unique, sqrt(sigmasq_unique))
    fit_all = fit_all + as.vector(unlist(fit_i))
  }
  
  fit_avar = fit_all/n_iter
  
  # plot(mu_mcmc[,1])
  # plot(sigmasq_mcmc[,1])
  # plot(omega_mcmc[,1])
  # hist(Y, prob=TRUE, ylim=c(0,0.5))
  # lines(x_range, fit_avar, col="red")
  
  return(fit_avar)
}

# 
# data_generation <- function(N, type="Gaussian mixture"){
#   
#   if(type=="Gaussian mixture"){
#     # data generation
#     weights_true = c(0.45,0.25,0.3)
#     K=3
#     c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
#     # hist(c_alloc_true)
#     
#     y = rep(NA, N)
#     # m = c(20,10,0,5,-2)
#     m = c(-2,0,3)
#     # m = c(7,10,0,5,-2)
#     s = c(0.4, 1, 0.3)
#     for (i in 1:K) {
#       I = c_alloc_true==i
#       y[I] = rnorm(sum(I), m[i], s[i])
#     }
#     hist(y)
#   }else if(type=="Skew Gaussian mixture"){
#     weights_true = c(0.45,0.25,0.3)
#     K=3
#     c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
#     y = rep(NA, N)
#     m = c(-2,0,3)
#     s = c(0.4, 1, 0.3)
#     a = c(1,10,4)
#     
#     for (i in 1:K) {
#       I = c_alloc_true==i
#       y[I] = rsn(sum(I),  m[i], s[i], a[i])
#     }
#     plot(density(y))
#   }else if(type=="Skew-symmetric mixture"){
#     weights_true = c(0.45,0.25,0.3)
#     K=3
#     c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
#     y = rep(NA, N)
#     
#     I2 = c_alloc_true==2
#     y[I2] = rsn(sum(I2),  0, 5, 4)
#     I3 = c_alloc_true==3
#     y[I3] = rnorm(sum(I3),  -4, sqrt(0.5))
#     
#     I1 = c_alloc_true==1
#     set1 = which(I1)
#     w11 = c(0.364,0.212,0.424)
#     c11 = sample(1:3, sum(I1), replace = TRUE, prob=w11)
#     I11 = c11==1
#     y[set1[I11]] = rsn(sum(I11),  2.5, 1, -10)
#     I12 = c11==2
#     y[set1[I12]] = rnorm(sum(I12),  2.325, sqrt(0.2))
#     I13 = c11==3
#     y[set1[I13]] = rnorm(sum(I13),  1.085, sqrt(0.7))
#     
#     plot(density(y))
#   }
#   
#   
#   return(list(data=y, c_true = c_alloc_true))
# }

# 
# # prior
# alpha = 1
# kappa0 = 0.01
# theta = 0
# a = 2
# b = 5
# prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)
# 
# sample_i = data_generation(N, type="Gaussian mixture")
# data_i = sample_i$data
# c_true_i = sample_i$c_true
# 
# fit_avar = MCMC_fit(n_iter=1000, c_true_i,data_i, prior_list, len_range=5000) 
