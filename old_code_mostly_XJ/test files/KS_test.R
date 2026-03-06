# Kolmogorov Smirnov distance


set.seed(1)

n=100
y = rnorm(n, 2, 1)
Fn = ecdf(y)
knots_vec = knots(Fn)

ecdf_func <- function(data) { 
  Length <- length(data) 
  sorted <- sort(data) 
  
  ecdf <- rep(0, Length) 
  for (i in 1:n) { 
    ecdf[i] <- sum(sorted <= data[i]) / Length 
  } 
  return(ecdf) 
} 

F_N_ecdf <- ecdf_func(y) 
F_N_ecdf_sort = sort(F_N_ecdf)

plot(y, F_N_ecdf)

plot(Fn)

F_parametric <- function(x, mean, sd){
  f = pnorm(x, mean = mean, sd = sd)
  return(f)
}


mu = 2
sigma = 1

distance_vec = numeric(n-1)
distance_vec[1] = abs(F_parametric(knots_vec[1], mu, sigma) - F_N_ecdf_sort[1])

for(i in 2:n){
  x_i = knots_vec[i-1]
  x_i_next = knots_vec[i]
  
  d1 = abs(F_parametric(x_i, mu, sigma) - F_N_ecdf_sort[i])
  d2 = abs(F_parametric(x_i_next, mu, sigma) - F_N_ecdf_sort[i])
  distance_vec[i] = max(d1, d2)
}

KS = max(distance_vec)
print(KS)

KS_distance_uninorm <- function(observed_data, mu0, sd0){
  
  ecdf_func <- function(data) { 
    Length <- length(data) 
    sorted <- sort(data) 
    
    ecdf <- rep(0, Length) 
    for (i in 1:n) { 
      ecdf[i] <- sum(sorted <= data[i]) / Length 
    } 
    return(ecdf) 
  } 
  
  n = length(observed_data)
  FN = ecdf_func(observed_data)
  FN_sort = sort(FN)
  
  knots_vec = sort(observed_data)
  distance_vec = numeric(n-1)
  distance_vec[1] = abs(pnorm(knots_vec[1], mu0, sd0) - FN_sort[1])
  
  for(i in 2:n){
    x_i = knots_vec[i-1]
    x_i_next = knots_vec[i]
    
    d1 = abs(pnorm(x_i, mu0, sd0) - FN_sort[i])
    d2 = abs(pnorm(x_i_next, mu0, sd0) - FN_sort[i])
    distance_vec[i] = max(d1, d2)
  }
  
  KS = max(distance_vec)
  return(KS)
}

KS_test1 = KS_distance_uninorm(y, 2,1)
print(KS_test1)


Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

KS_distance_mixture <- function(observed_data, w, u, s, cdf_func){
  
  ecdf_func <- function(data) { 
    Length <- length(data) 
    sorted <- sort(data) 
    ecdf <- rep(0, Length) 
    for (i in 1:n) { 
      ecdf[i] <- sum(sorted <= data[i]) / Length 
    } 
    return(ecdf) 
  } 
  
  n = length(observed_data)
  FN = ecdf_func(observed_data)
  FN_sort = sort(FN)
  knots_vec = sort(observed_data)
  distance_vec = numeric(n-1)
  distance_vec[1] = abs(cdf_func(knots_vec[1], w, u, s) - FN_sort[1])
  for(i in 2:n){
    x_i = knots_vec[i-1]
    x_i_next = knots_vec[i]
    d1 = abs(cdf_func(x_i,w, u, s) - FN_sort[i])
    d2 = abs(cdf_func(x_i_next,w, u, s) - FN_sort[i])
    distance_vec[i] = max(d1, d2)
  }
  KS = max(distance_vec)
  return(KS)
}



mixture_posterior <- function(c_alloc, data, prior_list){
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
  
  # posterior mean
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
  
  posterior_list = list("weights"=weights, "mu_post"=mu_post, "sigma2_post"=sigma2_post)
  return(posterior_list)
}



# data generation
N = 300
c_alloc_true = sample(1:5, N, replace = TRUE)
hist(c_alloc_true)
K=5
y = rep(NA, N)
m = c(20,10,0,5,-2)
# m = c(-2,0,3)
s = c(0.4, 1, 0.3,1.1, 1)
for (i in 1:K) {
  I = c_alloc_true==i
  y[I] = rnorm(sum(I), m[i], s[i])
}
hist(y)
prior_list = list("alpha"=alpha, "kappa0"=kappa0, "theta"=theta, "a"=a, "b"=b)

c_alloc_test = sample(1:5, N, replace = TRUE)
posterior_list = mixture_posterior(c_alloc_test, y, prior_list)
posterior_list2 = mixture_posterior(c_alloc_true, y, prior_list)

KS_test2 = KS_distance_mixture(y, posterior_list$weights, posterior_list$mu_post, posterior_list$sigma2_post, Mixture_normal_cdf)
print(KS_test2)
KS_test3 = KS_distance_mixture(y, posterior_list2$weights, posterior_list2$mu_post, posterior_list2$sigma2_post, Mixture_normal_cdf)
print(KS_test3)
