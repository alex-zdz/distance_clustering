rm(list = ls())

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

# # posterior
# alpha_post=numeric(K)
# cluster=sort(unique(c_alloc))
# for(i in 1:K){
#   alpha_post[i]=alpha+sum(c_alloc==cluster[i])
# }
# alpha_post
# weights=alpha_post/sum(alpha_post) ##mean of Dir
# ###parameters for kernel
# kappa_post = rep(NA,3)
# mu_post = rep(NA,3)
# sigma2_post = rep(NA,3)
# a_post = rep(NA,3)
# b_post = rep(NA,3)
#     
# for(i in 1:K){
#   I = c_alloc==cluster[i]
#   nk=sum(I)
#   yk=y[I]
#   kappa_post[i] = kappa0 +nk
#   mu_post[i] = (kappa0* theta + sum(yk))/kappa_post[i]
#   a_post[i] = a+ nk*0.5
#   b_post[i] = b + 0.5*sum((yk-mean(yk))^2)+0.5*kappa0*nk*(mean(yk)-theta)^2/kappa_post[i]
#   sigma2_post[i] = b_post[i]/(a_post[i] -1)
# }



T_function <- function(x){
  Tx = - (2*pi)^(-1/2) * exp(-x^2/2)
  return(Tx)
}

delta_T <- function(xi_n, xi_n_minus, mu_k, sigma2_k){
  result = T_function((xi_n-mu_k)/sigma2_k) - T_function((xi_n_minus-mu_k)/sigma2_k)
  return(result)
}

delta_F <- function(xi_n, xi_n_minus, mu_k, sigma2_k){
  f1 = pnorm(xi_n, mean = mu_k, sd = sqrt(sigma2_k), lower.tail = TRUE, log.p = FALSE)
  f2 = pnorm(xi_n_minus, mean = mu_k, sd = sqrt(sigma2_k), lower.tail = TRUE, log.p = FALSE)
  result = f1 - f2
  return(result)
}

F = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )
# provide an initial bracket for the quantile. default is c(-1000,1000). 
F_inv = function(p,w,u,s,br=c(-1000,1000))
{
  G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}

W2 <- function(data, c_alloc, prior_list){
  
  # prior
  alpha = prior_list$alpha
  kappa0 = prior_list$kappa0
  theta = prior_list$theta
  a = prior_list$a
  b = prior_list$b
  
  
  cluster=sort(unique(c_alloc))
  K = length(cluster)
  N = length(data)
  data_ordered = sort(data)
  
  alpha_post = numeric(K)
  mu_post = numeric(K)
  sigma2_post = numeric(K)

  
  for (k in 1:K) {
    I = c_alloc==cluster[k]
    nk=sum(I)
    yk=data[I]
    alpha_post[k]=alpha+sum(c_alloc==cluster[k])
    kappa_post = kappa0 +nk
    mu_post[k] = (kappa0* theta + sum(yk))/kappa_post
    a_post = a + nk*0.5
    b_post = b + 0.5*sum((yk-mean(yk))^2)+0.5*kappa0*nk*(mean(yk)-theta)^2/kappa_post
    sigma2_post[k] = b_post/(a_post -1)
  }
  weights=alpha_post/sum(alpha_post)
  # mu_post = m
  # sigma2_post = s
  
  sum_F = 0
  sum_T = 0
  part2 = numeric(K)
  part3 = numeric(K)
  for (i_n in 1:N)  {
    xi_n = F_inv(i_n/N,weights,mu_post, sqrt(sigma2_post))
    xi_n_minus = F_inv((i_n-1)/N,weights,mu_post, sqrt(sigma2_post))
    for (k in 1:K){
      sum_F = sum_F + data_ordered[i_n] * delta_F(xi_n, xi_n_minus, mu_post[k], sigma2_post[k])
      sum_T = sum_T + data_ordered[i_n] * delta_T(xi_n, xi_n_minus, mu_post[k], sigma2_post[k])
    }
    part2[k] = mu_post[k]^2 + sigma2_post[k]
    part3[k] = mu_post[k] * sum_F + sigma2_post[k] * sum_T
  }
  
  result = sum(data ^ 2) + sum(weights * part2) - 2 * sum(weights * part3)
  return(result)
}



prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)

# test W2
c_alloc = sample(1:3, N, replace = TRUE)
W2(data=y, c_alloc, prior_list)
W2(data=y, c_alloc_true, prior_list)


# X = c(rnorm(5000), rnorm(2500,mean=2,sd=1),rnorm(2000,mean=5,sd=1),rnorm(500,mean=10,sd=1))
# quantile(X,.95)
# F_inv(.95,c(.5,.25,.2,.05),c(0,2,5,10),c(1,1,1,1))
