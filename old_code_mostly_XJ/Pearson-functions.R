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
    a_post = a + (nk + 1) * 0.5 # if mean depends on sigma^2, need to add 1!
    b_post = b + 0.5*sum((yk-mean(yk))^2)+0.5*kappa0*nk*(mean(yk)-theta)^2/kappa_post
    sigma2_post[k] = b_post/(a_post -1) 
  }
  weights= alpha_post/sum(alpha_post)
  
  posterior_list = list("weights"=weights, "mu_post"=mu_post, "sigma2_post"=sigma2_post)
}

Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )


Pearson_distance_mixture <- function(observed_data, c_alloc, prior_list){
  # Pearson distance between observed data and an univariate mixture normal distribution with weights w, mean vec u, sd vec s
  
  posterior_list = mixture_posterior(c_alloc, observed_data, prior_list)
  w = posterior_list$weights
  u = posterior_list$mu_post
  s = posterior_list$sigma2_post
  
  
  n = length(observed_data)
  Fn = ecdf(observed_data)
  
  # bins1 = seq(-50,-15,by=1)
  # bins2 = seq(-15,15,by=0.2)
  # bins3 = seq(15,50,by=1)
  # bins = c(bins1, bins2, bins3)
  bins = seq(-50,50,by=0.5)
  n_bins = length(bins)+1
  ff = numeric(n_bins)
  gg = numeric(n_bins)
  gg[1] = sum(knots(Fn)<=bins[1])/n
  gg[n_bins] = sum(knots(Fn)>bins[n_bins-1])/n
  ff[1] = Mixture_normal_cdf(bins[1], w, u, s)
  ff[n_bins] = 1-Mixture_normal_cdf(bins[n_bins-1],w, u, s)

  
  for (i in 2:(n_bins-1)) {
    gg[i] = sum((knots(Fn)<=bins[i])*  (knots(Fn)>bins[i-1]))/n
    ff[i] = Mixture_normal_cdf(bins[i], w, u, s) - Mixture_normal_cdf(bins[i-1], w, u, s)
  }
  
  for (i in 1:n_bins) {
    if(ff[i]==0){
      ff[i] = 1e-50
    }
  }
  
  # print(gg)
  # print(ff)
  
  result = sum((gg-ff)^2/ff)
  
  if(is.na(result)){
    print("check pearson")
    print(((gg-ff)^2/ff))
    print(ff[n_bins-1])
  }

  # print(ff[n_bins])
  # print(Mixture_normal_cdf(bins[n_bins-1],w, u, s))
  
  return(result) 
  
}


