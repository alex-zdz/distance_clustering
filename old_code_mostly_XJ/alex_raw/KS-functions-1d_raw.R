Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x, mean = u, sd = sqrt(s)) )

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
    mu_post[k] = (kappa0 * theta + sum(yk))/kappa_post
    a_post = a + (nk + 1) * 0.5 # if mean depends on sigma^2, need to add 1!
    b_post = b + 0.5 * sum((yk-mean(yk))^2) + 0.5*kappa0*nk*(mean(yk)-theta)^2/kappa_post
    sigma2_post[k] = b_post/(a_post -1) 
  }
  weights= alpha_post/sum(alpha_post)   
  
  posterior_list = list("weights" = weights, "mu_post" = mu_post, "sigma2_post" = sigma2_post)
  return(posterior_list)
}


KS_distance_mixture <- function(observed_data, c_alloc, prior_list){
  
  posterior_list = mixture_posterior(c_alloc, observed_data, prior_list)
  w = posterior_list$weights
  u = posterior_list$mu_post
  s = posterior_list$sigma2_post
  
  n = length(observed_data)
  sorted <- sort(observed_data) 
  FN <- rep(0, n) 
  for (i in 1:n) { 
    FN[i] <- sum(sorted <= observed_data[i]) / n 
  } 
    
  FN_sort = sort(FN)
  knots_vec = sort(observed_data)
  distance_vec = numeric(n-1)
  distance_vec[1] = abs(Mixture_normal_cdf(knots_vec[1], w, u, s) - FN_sort[1])
  for(i in 2:n){
    x_i = knots_vec[i-1]
    x_i_next = knots_vec[i]
    d1 = abs(Mixture_normal_cdf(x_i,w, u, s) - FN_sort[i])
    d2 = abs(Mixture_normal_cdf(x_i_next,w, u, s) - FN_sort[i])
    distance_vec[i] = max(d1, d2)
  }
  KS = max(distance_vec)
  return(KS)
}



