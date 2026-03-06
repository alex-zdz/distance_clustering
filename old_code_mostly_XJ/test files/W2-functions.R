T_function <- function(x){
  Tx = - (2*pi)^(-1/2) * exp(-x^2/2)
  return(Tx)
}

delta_T <- function(xi_n, xi_n_minus, mu_k, sigma2_k){
  result = T_function((xi_n-mu_k)/sqrt(sigma2_k)) - T_function((xi_n_minus-mu_k)/sqrt(sigma2_k))
  return(result)
}

delta_F <- function(xi_n, xi_n_minus, mu_k, sigma2_k){
  f1 = pnorm(xi_n, mean = mu_k, sd = sqrt(sigma2_k), lower.tail = TRUE, log.p = FALSE)
  f2 = pnorm(xi_n_minus, mean = mu_k, sd = sqrt(sigma2_k), lower.tail = TRUE, log.p = FALSE)
  result = f1 - f2
  return(result)
}

F = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

G_func <- function(x, func_param) {
  w = func_param$w
  u = func_param$u
  s = func_param$s
  p = func_param$p
  res = F(x,w,u,s) - p
  return(res)
}

myUniroot <- function(func_param, lower, upper, tol = 1.0e-9) {
  f_lower <- G_func(lower, func_param)
  f_upper <- G_func(upper, func_param)
  
  while (abs(upper - lower) > tol) {
    mid <- (lower + upper) / 2.0
    f_mid <- G_func(mid, func_param)
    
    if (f_mid == 0.0) {
      return(mid)
    } else if (f_lower * f_mid <= 0) {
      upper <- mid
    } else {
      lower <- mid
      f_lower <- f_mid
    }
  }
  return((lower + upper) / 2.0)
}

F_inv = function(p,w,u,s,br=c(-1000,1000)){
  F_param = list("p"=p, "w"=w, "u"=u, "s"=s)
  result = myUniroot(F_param, br[1], br[2], tol = 1.0e-9)
  return(result)
}


# provide an initial bracket for the quantile. default is c(-1000,1000).
# F_inv = function(p,w,u,s,br=c(-1000,1000)){
#   G = function(x) F(x,w,u,s) - p
#   result = tryCatch({
#     uniroot(G,br)$root
#   }, warning = function(w) {
#     print("warning")
#   }, error = function(e) {
#     print("error")
#     x = seq(from = -100, to = 100,length.out=500)
#     y = numeric(length(x))
#     for (i in 1:length(x)) {
#       y[i] = G(x[i])
#     }
#     plot(x,y)
#   })
#   return( result )
# }

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
  # mu_post = m
  # sigma2_post = s
  # print("weights")
  # print(weights)

  F_mat = matrix(NA, nrow=N, ncol=K)
  T_mat = matrix(NA, nrow=N, ncol=K)
  part2 = numeric(K)
  part3 = numeric(K)
  for (i_n in 1:N)  {
    xi_n = F_inv(i_n/N,weights,mu_post, sqrt(sigma2_post))
    xi_n_minus = F_inv((i_n-1)/N,weights,mu_post, sqrt(sigma2_post))
    
    # print("input4-R")
    # print(sqrt(sigma2_post))
    # print("xi_n")
    # print(xi_n)
    # print("xi_n_minus")
    # print(xi_n_minus)
    
    for (k in 1:K){
      F_mat[i_n, k] = delta_F(xi_n, xi_n_minus, mu_post[k], sigma2_post[k])
      T_mat[i_n, k] = delta_T(xi_n, xi_n_minus, mu_post[k], sigma2_post[k])
    }
  }
  
  sum_F_k = rep(0, K)
  sum_T_k = rep(0, K)
  for (k in 1:K) {
    for (i_n in 1:N) {
      sum_F_k[k] = sum_F_k[k] + data_ordered[i_n] * F_mat[i_n, k]
      sum_T_k[k] = sum_T_k[k] + data_ordered[i_n] * T_mat[i_n, k]
    }
    part2[k] = mu_post[k]^2 + sigma2_post[k]
    part3[k] = mu_post[k] * sum_F_k[k] + sqrt(sigma2_post[k]) * sum_T_k[k]
    # print(mu_post[k] * sum_F_k[k])
  }

  
  # print("F_mat--R")
  # print(F_mat[1:3,1:K])
  # print("T_mat--R")
  # print(T_mat[1:3,1:K])
  
  # print("sum_F_k")
  # print(sum_F_k)
  # print("part3")
  # print(part3)
  
  # result = - 2 * sum(weights * part3)
  result = (1/N) * sum(data ^ 2) + sum(weights * part2) - 2 * sum(weights * part3)
  return(result)
}


min_W2 <- function(c_samples, data, prior_list){
  # c_samples: M*n matrix
  M = dim(c_samples)[1]
  W2_vec = numeric(M)
  for (i in 1:M) {
    if(i%%100==0){
      print(i)
    }
    W2_vec[i] = W2(data=data, c_samples[i,], prior_list)
  }
  min_m = which(W2_vec==min(W2_vec))[1]
  # print(min_m)
  return(c_samples[min_m,])
}
