library(minpack.lm)

# Define a mixture model
mixture_model_pdf <- function(x, weight1, mean1, sd1, mean2, sd2) {
  f = weight1 * dnorm(x, mean1, sd1) + (1-weight1) * dnorm(x, mean2, sd2)
  return(f)
}
mixture_model_cdf <- function(x, weight1, mean1, sd1, mean2, sd2) {
  f = weight1 * pnorm(x, mean1, sd1) + (1-weight1) * pnorm(x, mean2, sd2)
  return(f)
}


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

G_n = ecdf(observed_data)
x <- seq(-5, 10, length.out = 200)
G_n_x = G_n(x)


# Use nlsLM to optimize the parameters
initial_params <- list(weight1=0.5, mean1=1, sd1=1.2, mean2=5, sd2=1)  # Initial parameter values

fit <- nlsLM(G_n_x ~ mixture_model_cdf(x, weight1, mean1, sd1, mean2, sd2), start = initial_params, lower = c(0,-Inf,0, -Inf,0), upper = c(1, Inf, Inf, Inf, Inf))

# Print the result
print(summary(fit))


# Formula: G_n_x ~ mixture_model_cdf(x, weight1, mean1, sd1, mean2, sd2)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# weight1 0.303694   0.000952   319.0   <2e-16 ***
#   mean1   0.978473   0.008411   116.3   <2e-16 ***
#   sd1     1.102149   0.010592   104.1   <2e-16 ***
#   mean2   5.042021   0.002094  2407.6   <2e-16 ***
#   sd2     0.506103   0.002761   183.3   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.003249 on 195 degrees of freedom
# 
# Number of iterations to convergence: 7 
# Achieved convergence tolerance: 1.49e-08