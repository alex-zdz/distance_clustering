

#Create grid:
# Define parameter grids
method_vals   <- c("KS", "Pearson", "W2")
n_runs_vals   <- c(1, 3)         # small to moderate c(1, 10, 100) 
n_sweet_vals  <- c(1, 5, 10)    # short vs. longer sweeps
n_ms_vals <- c(1, 5, 10)         # max local moves
n_merge_vals  <- c(1, 10)         
n_split_vals  <- c(1, 10)
data_lengths  <- c(100, 5e2)

# Create full parameter grid
param_grid <- expand.grid(
  method   = method_vals,
  n_runs   = n_runs_vals,
  n_sweet  = n_sweet_vals,
  n_ms = n_ms_vals,
  n_merge  = n_merge_vals,
  n_split  = n_split_vals,
  data_lengths = data_lengths,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

current_date <- format(Sys.Date(), "%d%m%Y")
# Save the grid
dir.create(paste0("results/simstudy/", current_date), showWarnings = FALSE, recursive = TRUE)
saveRDS(param_grid, paste0("results/simstudy/", current_date,"/grid.rds"))

# if(FALSE){
#   run = 1
#   dataset = "simdata_1"
#   n_save = 5e2;  n_burn = 10e2; n_thin = 1
# }

grid_eval <- function(run, param_grid){
  
  # Load libraries
  library(AntMAN)
  library(Rcpp)
  library(RcppArmadillo)
  library(fossil)
  # Source functions
  source("src/mixture_utils.R")
  Rcpp::sourceCpp("src/W2-functions.cpp")
  
  #pars <- as.list(param_grid[run, , drop = FALSE ])
  method   <- param_grid[run,]$method
  n_runs   <- param_grid[run,]$n_runs
  n_sweet  <- param_grid[run,]$n_sweet
  tol      <- param_grid[run,]$tol
  n_ms <- param_grid[run,]$n_ms
  n_merge  <- param_grid[run,]$n_merge
  n_split  <- param_grid[run,]$n_split
  N <- param_grid[run,]$data_lengths
  
  K <- 3
  alpha_0 <- 1
  # fix the seed
  set.seed(1)
  
  simdata <- generate_mixture_data(
    N = N,
    K = K,
    type = "Gaussian",
    alpha = alpha_0,
    means = NULL,
    sds = NULL,
    skews = NULL
  )
  
  y <- simdata$data
  cluster_true <- simdata$cluster_true
  
  # -----------------------------
  # Prior hyperparameters
  # -----------------------------
  alpha <- 1
  kappa_0 <- 1
  sig2_0 <- 1
  mu_0 <- 0
  nu_0 <- 1
  
  # -----------------------------
  # MCMC settings
  # -----------------------------
  n_iter <- 2e4
  burn_in <- n_iter / 4
  thin <- 10
  
  # -----------------------------
  # Define ANTMAN priors & parameters
  # -----------------------------
  mixture_uvn_params <- AntMAN::AM_mix_hyperparams_uninorm(
    m0 = mu_0, k0 = kappa_0, nu0 = nu_0, sig02 = sig2_0
  )
  
  mcmc_params <- AntMAN::AM_mcmc_parameters(
    niter = n_iter, burnin = burn_in, thin = thin, verbose = 1
  )
  
  components_prior <- AntMAN::AM_mix_components_prior_pois(Lambda = K)
  
  # -----------------------------
  # Run MCMC and cluster
  # -----------------------------
  mix_post_draws <- AntMAN::AM_mcmc_fit(
    y = y,
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior = components_prior,
    mcmc_parameters = mcmc_params
  )
  
  clustering_matrix <- AM_clustering(mix_post_draws) + 1 # M*n matrix
  
  posterior_params <- prepare_AM_posterior_params(mix_post_draws)
  
  prior_list <- list(
    alpha_0 = 1,      
    kappa_0 = 1,    
    mu_0    = 0,       
    a_0     = 2,      
    b_0     = 1   
  )
  
  
  ################################################################################
  
  print(y)
  
  startt <- Sys.time()
  
  err_clustering <- "no_error"
  err_ari        <- "no_error"
  startt <- Sys.time()
  
  # Run clustering safely
  results <- tryCatch({
    run_clustering(
      data              = y,
      clustering_matrix = clustering_matrix,
      posterior_params  = posterior_params,
      prior_list        = prior_list,
      method            = method,
      n_runs            = n_runs,
      n_sweet           = n_sweet,
      tol               = 1e-10,
      n_ms              = n_ms,
      n_merge           = n_merge,
      n_split           = n_split
    )
  }, error = function(e) {
    err_clustering <- conditionMessage(e)
    NULL
  })
  
  elapsed_sec <- as.numeric(Sys.time() - startt, units = "secs")
  
  # Compute Adjusted Rand Index safely
  adj_rand <- tryCatch({
    if (!is.null(results) && !is.null(results$clustering)) {
      adj.rand.index(cluster_true, results$clustering)
    } else {
      NA_real_
    }
  }, error = function(e) {
    err_adj_rand <- conditionMessage(e)
    NA_real_
  })
  
  # -----------------------------
  # Store results in a data frame
  # -----------------------------
  df_out <- data.frame(
    run              = run,
    time_sec         = elapsed_sec,
    adj_rand         = adj_rand,
    err_clustering   = err_clustering,
    err_adj_rand     = err_ari,
    stringsAsFactors = FALSE
  )
  
  current_date <- format(Sys.Date(), "%d%m%Y")
  #filename <- paste0( paste0(colnames(param_grid[run,]), "_"), paste0(param_grid[run,]), collapse = "_" )
  filename <- paste0(paste0(param_grid[run,]), collapse = "_" )
  
  saveRDS(df_out, paste0("results/simstudy/", current_date,"/",filename,".rds"))
  
}



# Load and source, simulate data

library(pbapply)
library(parallel)
parallel::detectCores()
cl <- parallel::makeCluster(parallel::detectCores())

op <- pboptions(type="timer")

# Quick run for testing

#system.time(pblapply(1:10, grid_eval, param_grid = param_grid, cl = cl))

system.time(pblapply(1:nrow(param_grid), grid_eval, param_grid = param_grid, cl = cl))

parallel::stopCluster(cl)




