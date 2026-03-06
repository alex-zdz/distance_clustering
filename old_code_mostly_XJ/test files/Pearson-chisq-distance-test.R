# Pearson chi-sq distance


n=100
y = rnorm(n, 2, 1)
Fn = ecdf(y)
knots(Fn)

bins = seq(-5,5,by=1)
n_bins = length(bins)+1
ff = numeric(n_bins)
gg = numeric(n_bins)
gg[1] = sum(knots(Fn)<=bins[1])/n
gg[n_bins] = sum(knots(Fn)>bins[n_bins-1])/n
ff[1] = pnorm(bins[1])
ff[n_bins] = 1-pnorm(bins[n_bins-1])

for (i in 2:(n_bins-1)) {
  gg[i] = sum((knots(Fn)<=bins[i])*  (knots(Fn)>bins[i-1]))/n
  ff[i] = pnorm(bins[i]) - pnorm(bins[i-1])
}


# gg[2] = sum((knots(Fn)<=bins[2]) * (knots(Fn)>bins[1]))/n

Pearson_distance_result1 = sum((gg-ff)^2/ff)


Pearson_distance_uninorm <- function(observed_data, mu0, sd0){
  # Pearson distance between observed data and an univariate normal distribution N(mu0, sd0)
  
  n = length(observed_data)
  Fn = ecdf(observed_data)
  
  bins = seq(-50,50,by=0.5)
  n_bins = length(bins)+1
  ff = numeric(n_bins)
  gg = numeric(n_bins)
  gg[1] = sum(knots(Fn)<=bins[1])/n
  gg[n_bins] = sum(knots(Fn)>bins[n_bins-1])/n
  ff[1] = pnorm(bins[1], mean=mu0, sd=sd0)
  ff[n_bins] = 1-pnorm(bins[n_bins-1], mean=mu0, sd=sd0)
  
  for (i in 2:(n_bins-1)) {
    gg[i] = sum((knots(Fn)<=bins[i])*  (knots(Fn)>bins[i-1]))/n
    ff[i] = pnorm(bins[i], mean=mu0, sd=sd0) - pnorm(bins[i-1], mean=mu0, sd=sd0)
  }
  

  
  result = sum((gg-ff)^2/ff)
  return(result) 
}

Pearson_distance_result2 = Pearson_distance_uninorm(y, 0, 1)



Pearson_distance_mixture <- function(observed_data, w, u, s, cdf_func){
  # Pearson distance between observed data and an univariate mixture normal distribution with weights w, mean vec u, sd vec s
  
  n = length(observed_data)
  Fn = ecdf(observed_data)
  
  bins = seq(-5,5,by=1)
  n_bins = length(bins)+1
  ff = numeric(n_bins)
  gg = numeric(n_bins)
  gg[1] = sum(knots(Fn)<=bins[1])/n
  gg[n_bins] = sum(knots(Fn)>bins[n_bins-1])/n
  ff[1] = cdf_func(bins[1], w, u, s)
  ff[n_bins] = 1-cdf_func(bins[n_bins-1],w, u, s)
  
  for (i in 2:(n_bins-1)) {
    gg[i] = sum((knots(Fn)<=bins[i])*  (knots(Fn)>bins[i-1]))/n
    ff[i] = cdf_func(bins[i], w, u, s) - cdf_func(bins[i-1], w, u, s)
  }
  
  result = sum((gg-ff)^2/ff)
  return(result) 
  
}


Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )



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


c_alloc_test = sample(1:5, N, replace = TRUE)
posterior_list = mixture_posterior(c_alloc_test, y, prior_list)

Pearson_distance_result3 = Pearson_distance_mixture(y, posterior_list$weights, posterior_list$mu_post, posterior_list$sigma2_post, Mixture_normal_cdf)
