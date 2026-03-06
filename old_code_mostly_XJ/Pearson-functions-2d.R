library(mvtnorm)
library('plot.matrix')
library(mltools)
library(data.table)

mixture_posterior_2d <- function(c_alloc, data, prior_list){
  # prior
  alpha = prior_list$alpha
  kappa0 = prior_list$kappa0
  m0 = prior_list$m0
  nu0 = prior_list$nu0
  Lambda0 = prior_list$Lambda0
  
  
  cluster=sort(unique(c_alloc))
  K = length(cluster)
  N = dim(data)[1]
  d = dim(data)[2]
  

  alpha_post=numeric(K)
  mu_post = matrix(NA, K,d)
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
  
  # posterior mean
  for(i in 1:K){
    I = c_alloc==cluster[i]
    nk=sum(I)
    yk=y[I,]
    Sk = Sk_calculation(yk)  # [d,d]
    alpha_post[i]=alpha+nk
    
    kappa_post_k = kappa0 +nk
    nu_post_k = nu0 + nk
    mu_post[i,] = (kappa0 * m0 + colSums(yk))/kappa_post_k
    Lambda_post_k = Lambda0 + Sk + t(t(colMeans(yk)-m0)) %*% t((colMeans(yk)-m0)) * (kappa0*nk/kappa_post_k)
    Lambda_post[[i]] = Lambda_post_k
    Sigma_post[[i]] = (Lambda_post_k)/(nu_post_k - d -1)
  }
  
  weights=alpha_post/sum(alpha_post)
  
  posterior_list = list("weights"=weights, "mu_post"=mu_post, "sigma2_post"=Sigma_post)
  return(posterior_list)
}

Mixture_normal_cdf_1d = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=sqrt(s)) )


Mixture_normal_cdf_2d <- function(x,w,u,s){
  K = length(s)
  cdf = 0
  for (k in 1:K) {
    mu_k = u[k,]
    Sigma_k = s[[k]]
    cdf = cdf + w[k] * max(pmvnorm(lower=-Inf,upper=x, mean=mu_k,sigma=Sigma_k),0)
  }
  return(cdf)
} 


Pearson_distance_mixture_2d <- function(observed_data, c_alloc, prior_list, bins=seq(-5,5,by=0.5)){
  # Pearson distance between observed data and an 2d multivariate mixture normal distribution with weights w, mean u, Sigma list s
  
  posterior_list = mixture_posterior_2d(c_alloc, observed_data, prior_list)
  w = posterior_list$weights     # vec: len=K
  u = posterior_list$mu_post     # mat: [K, d], each row is the mean of a d dimensional normal component
  s = posterior_list$sigma2_post # list: len=K
  
  u_marginal_1 = u[,1]  # the marginal mean
  u_marginal_2 = u[,2]  # the marginal mean
  s_marginal_1 = c()    # the marginal sigmasq
  s_marginal_2 = c()    # the marginal sigmasq
  for (k in 1:length(w)){
    s_marginal_1 = c(s_marginal_1, s[[k]][1,1])
    s_marginal_2 = c(s_marginal_2, s[[k]][2,2])
  }
  
  
  n = dim(observed_data)[1]
  
  # bins1 = seq(-50,-15,by=1)
  # bins2 = seq(-15,15,by=0.2)
  # bins3 = seq(15,50,by=1)
  # bins = c(bins1, bins2, bins3)
  # bins = seq(-50,50,by=0.5)
  # bins = seq(-5,5,by=0.5)
  n_bins = length(bins)+1
  
  ff = matrix(NA,n_bins,n_bins)  # ff[i,j] represents the cdf of space A_ij
  gg = matrix(NA,n_bins,n_bins)  # gg[i,j] represents the Empirical cdf of space A_ij
  
  dt = data.table(x = observed_data[,1], y = observed_data[,2])
  gg[1,1] = empirical_cdf(dt, ubounds = data.table(x = bins[1], y = bins[1]))$CDF
  ff[1,1] = Mixture_normal_cdf_2d(c(bins[1],bins[1]), w, u, s)
  
  for (i in 2:(n_bins-1)) {
    gg[i,1] = empirical_cdf(dt, ubounds = data.table(x = bins[i], y = bins[1]))$CDF - 
              empirical_cdf(dt, ubounds = data.table(x = bins[i-1], y = bins[1]))$CDF
    ff[i,1] = Mixture_normal_cdf_2d(c(bins[i],bins[1]), w, u, s) - 
              Mixture_normal_cdf_2d(c(bins[i-1],bins[1]), w, u, s)
  }
  for (j in 2:(n_bins-1)) {
    gg[1,j] = empirical_cdf(dt, ubounds = data.table(x = bins[1], y = bins[j]))$CDF - 
              empirical_cdf(dt, ubounds = data.table(x = bins[1], y = bins[j-1]))$CDF
    ff[1,j] = Mixture_normal_cdf_2d(c(bins[1],bins[j]), w, u, s) - 
              Mixture_normal_cdf_2d(c(bins[1],bins[j-1]), w, u, s)
  }
  gg[n_bins,n_bins] = 1 - empirical_cdf(dt, ubounds = data.table(x = bins[n_bins-1], y = bins[n_bins-1]))$CDF
  gg[1,n_bins] = empirical_cdf(observed_data[,1], ubounds = bins[1])$CDF - 
                  empirical_cdf(dt, ubounds = data.table(x = bins[1], y = bins[n_bins-1]))$CDF
  gg[n_bins,1] = empirical_cdf(observed_data[,2], ubounds = bins[1])$CDF - 
                empirical_cdf(dt, ubounds = data.table(x = bins[n_bins-1], y = bins[1]))$CDF
  
  
  ff[n_bins,n_bins] = 1-Mixture_normal_cdf_2d(c(bins[n_bins-1],bins[n_bins-1]),w, u, s)
  ff[1,n_bins] = Mixture_normal_cdf_1d(bins[1],w,u_marginal_1,s_marginal_1) - 
                  Mixture_normal_cdf_2d(c(bins[1],bins[n_bins-1]),w, u, s)
  ff[n_bins,1] = Mixture_normal_cdf_1d(bins[1],w,u_marginal_2,s_marginal_2) - 
                  Mixture_normal_cdf_2d(c(bins[n_bins-1],bins[1]),w, u, s)
   
  
  for (i in 2:(n_bins-1)) {
    for (j in 2:(n_bins-1)){
      gg[i,j] = empirical_cdf(dt, ubounds = data.table(x = bins[i], y = bins[j]))$CDF - 
                empirical_cdf(dt, ubounds = data.table(x = bins[i], y = bins[j-1]))$CDF - 
                empirical_cdf(dt, ubounds = data.table(x = bins[i-1], y = bins[j]))$CDF +
                empirical_cdf(dt, ubounds = data.table(x = bins[i-1], y = bins[j-1]))$CDF
      ff[i,j] = Mixture_normal_cdf_2d(c(bins[i],bins[j]), w, u, s) - 
                Mixture_normal_cdf_2d(c(bins[i],bins[j-1]), w, u, s) - 
                Mixture_normal_cdf_2d(c(bins[i-1],bins[j]), w, u, s) + 
                Mixture_normal_cdf_2d(c(bins[i-1],bins[j-1]), w, u, s)
    }
  }
  
  for (i in 2:(n_bins-1)) {
    gg[i,n_bins] = empirical_cdf(observed_data[,1], ubounds = bins[i])$CDF - 
                    empirical_cdf(observed_data[,1], ubounds = bins[i-1])$CDF - 
                    empirical_cdf(dt, ubounds = data.table(x = bins[i], y = bins[n_bins-1]))$CDF + 
                    empirical_cdf(dt, ubounds = data.table(x = bins[i-1], y = bins[n_bins-1]))$CDF
    ff[i,n_bins] = Mixture_normal_cdf_1d(bins[i], w, u_marginal_1, s_marginal_1) - 
                    Mixture_normal_cdf_1d(bins[i-1], w, u_marginal_1, s_marginal_1) - 
                    Mixture_normal_cdf_2d(c(bins[i],bins[n_bins-1]), w, u, s) + 
                    Mixture_normal_cdf_2d(c(bins[i-1],bins[n_bins-1]), w, u, s)
      
  }
  for (j in 2:(n_bins-1)) {
    gg[n_bins,j] = empirical_cdf(observed_data[,2], ubounds = bins[j])$CDF - 
                    empirical_cdf(observed_data[,2], ubounds = bins[j-1])$CDF - 
                    empirical_cdf(dt, ubounds = data.table(x = bins[n_bins-1], y = bins[j]))$CDF + 
                    empirical_cdf(dt, ubounds = data.table(x = bins[n_bins-1], y = bins[j-1]))$CDF
    ff[n_bins,j] = Mixture_normal_cdf_1d(bins[j], w, u_marginal_2, s_marginal_2) - 
                    Mixture_normal_cdf_1d(bins[j-1], w, u_marginal_2, s_marginal_2) - 
                    Mixture_normal_cdf_2d(c(bins[n_bins-1],bins[j]), w, u, s) + 
                    Mixture_normal_cdf_2d(c(bins[n_bins-1],bins[j-1]), w, u, s)
  }
  
  for (i in 1:n_bins) {
    for (j in 1:n_bins){
      if(ff[i,j]<0){
        ff[i,j] = 0
      }
      if(gg[i,j]<0){
        gg[i,j] = 0
      }
    }
  }
  
  
  ff_de = ff
  for (i in 1:n_bins) {
    for (j in 1:n_bins){
      if(ff_de[i,j]==0){
        ff_de[i,j] = 1e-50
      }
    }
  }
  
  
  # # # print("gg")
  # plot(gg)
  # # # print("ff")
  # plot.new()
  # plot.new()
  # plot(ff)
  
  result = sum((gg-ff)^2/ff_de)
  # print("result")
  # plot(((gg-ff)^2/ff_de))
  # 
  
  
  if(is.na(result)){
    print("check pearson")
    print(((gg-ff)^2/ff))
    print(ff[n_bins-1])
  }
  
  # print(ff[n_bins])
  # print(Mixture_normal_cdf(bins[n_bins-1],w, u, s))
  
  return(result) 
  
}


