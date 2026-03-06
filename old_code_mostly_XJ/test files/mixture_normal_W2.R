# Generate some data from the true model
set.seed(123)
n=500
p_true=0.3
mean1_true=1
sd1_true=1
mean2_true=5
sd2_true=1/2
observed_data = rep(NA,n)
for (i in 1:n) {
  U = runif(n=1,0,1)
  if(U<p_true){
    observed_data[i] = rnorm(n=1, mean=mean1_true, sd = sd1_true)
  }else{
    observed_data[i] = rnorm(n=1, mean=mean2_true, sd = sd2_true)
  }
}

mixture_model_cdf <- function(x, weight1, mean1, sd1, mean2, sd2) {
  f = weight1 * pnorm(x, mean1, sd1) + (1-weight1) * pnorm(x, mean2, sd2)
  return(f)
}

Wasserstein2_sq <- function(observed_data, mixture_model_cdf, param_list){
  weight1 = param_list$weight1
  mean1 = param_list$mean1
  sd1 = param_list$sd1
  mean2 = param_list$mean2
  sd2 = param_list$sd2
  
  Y_ordered = sort(observed_data)
  vec_FminusG = rep(0, n)
  for (i in 1:n) {
    F_i = mixture_model_cdf(Y_ordered[i], weight1, mean1, sd1, mean2, sd2)
    G_i = i/n
    vec_FminusG[i] = F_i - G_i
  }
  W2_sq = sum((vec_FminusG)^2)
  return(W2_sq)
}

param_list = list(weight1=0.7,mean1=1,sd1=1,mean2=5, sd2=0.5)

Wasserstein2_sq(observed_data,mixture_model_cdf,param_list)
