library(mvtnorm)
# library('plot.matrix')
library(mltools)
library(data.table)

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






KS_distance_mixture_2d <- function(observed_data, c_alloc, prior_list){
  
  posterior_list = mixture_posterior_2d(c_alloc, observed_data, prior_list)
  w = posterior_list$weights     # vec: len=K
  u = posterior_list$mu_post     # mat: [K, d], each row is the mean of a d dimensional normal component
  s = posterior_list$sigma2_post # list: len=K
  
  
  n = dim(observed_data)[1]
  F_P = matrix(NA, n, n)
  F_N = matrix(NA, n, n)
  
  knots_d1 = sort(observed_data[,1])
  knots_d2 = sort(observed_data[,2])
  dt = data.table(x = observed_data[,1], y = observed_data[,2])
  # FN = ecdf_func(observed_data)
  # FN_sort = sort(FN)
  # knots_vec = sort(observed_data)
  # distance_vec = numeric(n-1)
  # distance_vec[1] = abs(Mixture_normal_cdf(knots_vec[1], w, u, s) - FN_sort[1])
  
  
  # Build a grid based on the 1st and 2nd dimension of all the data points
  # Calculate parametric CDF F_P, and empirical CDF F_N on the grid
  for(i in 1:n){
    x_i = knots_d1[i]
    for (j in 1:n) {
      y_j = knots_d2[j]
      F_P[i,j] = Mixture_normal_cdf_2d(x=c(x_i,y_j),w,u,s)
      F_N[i,j] = empirical_cdf(dt, ubounds = data.table(x = x_i, y = y_j))$CDF
    }
    # x_i = knots_vec[i-1]
    # x_i_next = knots_vec[i]
    # d1 = abs(Mixture_normal_cdf(x_i,w, u, s) - FN_sort[i])
    # d2 = abs(Mixture_normal_cdf(x_i_next,w, u, s) - FN_sort[i])
    # distance_vec[i] = max(d1, d2)
  }
  
  distance_mat = matrix(NA, n, n)
  for(i in 2:n){
    for (j in 2:n) {
      d1 = abs(F_N[i,j] - F_P[i,j])
      d2 = abs(F_N[i-1,j] - F_P[i,j])
      d3 = abs(F_N[i,j-1] - F_P[i,j])
      d4 = abs(F_N[i-1,j-1] - F_P[i,j])
      distance_mat[i,j] = max(d1,d2,d3,d4)
    }
  }
  distance_mat[1,1] = abs(F_N[1,1] - F_P[1,1])
  for (i in 2:n) {
    d1 = abs(F_N[i,1] - F_P[i,1])
    d2 = abs(F_N[i-1,1] - F_P[i,1])
    distance_mat[i,1] = max(d1,d2)
  }
  for (j in 2:n) {
    d1 = abs(F_N[1,j] - F_P[1,j])
    d2 = abs(F_N[1,j-1] - F_P[1,j])
    distance_mat[1,j] = max(d1,d2)
  }
  
  KS = max(distance_mat)
  return(KS)
}



