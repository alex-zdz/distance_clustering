# =============================================================================
# Simulation study – 2D Gaussian mixtures – evaluating search-based clustering
# =============================================================================
# Purpose: Compare performance of Wasserstein discrepancy + refinement strategies
# on two synthetic 2D scenarios using posterior-informed initializations
# Date: March 2026
# Author: Alex
#
# Dependencies:
# - AntMAN
# - mclust (for adjustedRandIndex)
# - salso (for binder, VI)
# - tidyverse, pbapply, parallel
# - Custom source files (updated versions without 'version' parameter)

# ---- 1. Clear environment & load libraries ----
rm(list = ls())
library(AntMAN)
library(mclust)       # adjustedRandIndex
library(salso)        # binder, VI
library(tidyverse)
library(pbapply)
library(parallel)

# ---- 2. Source the three main function files ----
source("src/search_algorithm.R")
source("src/mixture_utils.R")
source("src/distances.R")

# ---- 3. Data generation functions ----
generate_scenario1 <- function() {
  set.seed(1)
  generate_mixture_data(
    N = 100,
    K = 3,
    dim = 2,
    mu_true = matrix(c(-2, 0, 0, 0, 2, 0), nrow = 3, byrow = TRUE),
    Sigma_true = lapply(1:3, function(i) diag(c(0.2, 0.2)))
  )
}

generate_scenario2 <- function() {
  set.seed(2)
  generate_mixture_data(
    N = 150,
    K = 4,
    dim = 2,
    mu_true = matrix(c(-4, 0, -1, 0, 1, 0, 4, 0), nrow = 4, byrow = TRUE),
    Sigma_true = lapply(1:4, function(i) diag(c(0.5, 0.5)))
  )
}

# ---- 4. AntMAN posterior draws ----
get_antman_posterior <- function(y, K_true) {
  D <- ncol(y)
  alpha_0 <- 1
  kappa_0 <- 0.01
  nu_0    <- D + 2
  S_0     <- diag(D) * 2
  
  n_iter  <- 20000
  burn_in <- n_iter / 4
  thin    <- 10
  
  mixture_mvn_params <- AM_mix_hyperparams_multinorm(
    mu0 = rep(0, D), ka0 = kappa_0, nu0 = nu_0, Lam0 = S_0
  )
  
  mcmc_params <- AM_mcmc_parameters(
    niter = n_iter, burnin = burn_in, thin = thin, verbose = 0
  )
  
  components_prior <- AM_mix_components_prior_pois(Lambda = K_true)
  
  mix_post_draws <- AM_mcmc_fit(
    y = y,
    mix_kernel_hyperparams = mixture_mvn_params,
    mix_components_prior   = components_prior,
    mcmc_parameters        = mcmc_params
  )
  
  eam <- AM_clustering(mix_post_draws)
  clustering_matrix <- eam + 1L
  
  list(
    mix_post_draws    = mix_post_draws,
    clustering_matrix = clustering_matrix
  )
}

# ---- 5. Fixed prior list ----
prior_list <- list(
  alpha_0 = 1,
  kappa_0 = 0.01,
  mu_0    = c(0, 0),
  nu_0    = 4,
  S_0     = diag(2) * 5
)

# ---- 6. Parameter grid (no version anymore) ----
methods     <- "Wasserstein"
clusterwise <- c(FALSE, TRUE)
run_grid <- expand.grid(
  method      = methods,
  clusterwise = clusterwise,
  n_sweet     = c(1, 3),
  n_ms        = c(0, 3),
  n_merge     = 5,
  n_split     = 5,
  nu_0        = c(4, 7),
  scenario    = c("scen1", "scen2"),
  seed        = 1,
  stringsAsFactors = FALSE
)

cat("Total configurations:", nrow(run_grid), "\n")

# ---- 7. Precompute data, posteriors and fixed theta matrices ----
cat("Generating data, AntMAN posteriors and fixed projection directions ...\n")

out_dir <- "results/simstudy/2D"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

scen1_data <- generate_scenario1()
scen1_post <- get_antman_posterior(scen1_data$data, length(unique(scen1_data$cluster_true)))

scen2_data <- generate_scenario2()
scen2_post <- get_antman_posterior(scen2_data$data, length(unique(scen2_data$cluster_true)))

# Save full posterior objects once (needed for later Binder/VI loss computation in dashboard)
saveRDS(scen1_data, file.path(out_dir, "scen1_data.rds"))
saveRDS(scen2_data, file.path(out_dir, "scen2_data.rds"))
saveRDS(scen1_post, file.path(out_dir, "scen1_post.rds"))
saveRDS(scen2_post, file.path(out_dir, "scen2_post.rds"))
cat("Saved posterior objects: scen1_post.rds and scen2_post.rds\n")

# Fixed theta — same for all configurations within a scenario (L = 20)
L_fixed <- 20
set.seed(123)  # reproducible across runs
theta_scen1 <- generate_theta_normal(L = L_fixed, d = 2)
theta_scen2 <- generate_theta_normal(L = L_fixed, d = 2)

# ---- 8. Parallel worker function ----
run_one_config <- function(i, run_grid, prior_list, theta_scen1, theta_scen2) {
  row <- run_grid[i, ]
  set.seed(row$seed)
  
  # Select scenario data, posterior and theta
  if (row$scenario == "scen1") {
    scen_data <- scen1_data
    post      <- scen1_post
    theta     <- theta_scen1
  } else {
    scen_data <- scen2_data
    post      <- scen2_post
    theta     <- theta_scen2
  }
  
  y <- scen_data$data
  cluster_true <- scen_data$cluster_true
  
  # Override nu_0 for this configuration
  prior_list$nu_0 <- row$nu_0
  
  # Run clustering
  fit <- tryCatch(
    run_clustering(
      data              = y,
      clustering_matrix = post$clustering_matrix,
      posterior_draws   = post$mix_post_draws,
      prior_list        = prior_list,
      method            = row$method,
      clusterwise       = row$clusterwise,
      estimator         = "mle",
      n_runs            = 1,
      n_sweet           = row$n_sweet,
      tol               = 1e-10,
      n_ms              = row$n_ms,
      n_merge           = row$n_merge,
      n_split           = row$n_split,
      M                 = 2000,
      alpha_0           = prior_list$alpha_0,
      nu_0              = row$nu_0,
      L                 = nrow(theta),   # overridden by theta anyway
      theta             = theta
    ),
    error = function(e) {
      message(sprintf("Config %d (%s) failed: %s", i, row$scenario, e$message))
      return(NULL)
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  # Best run
  best_idx <- which.min(fit$distance_record)
  best_D   <- fit$distance_record[best_idx]
  best_c   <- fit$final_clustering
  
  # ARI vs truth
  ari <- adjustedRandIndex(best_c, cluster_true)
  
  # Posterior consensus clusterings (for ARI against consensus)
  binder_clust <- salso(post$clustering_matrix, loss = binder())
  vi_clust     <- salso(post$clustering_matrix, loss = VI())
  ari_binder   <- adjustedRandIndex(best_c, binder_clust)
  ari_vi       <- adjustedRandIndex(best_c, vi_clust)
  
  # Counts
  n_sweet_acc <- fit$total_accepted_sweets[best_idx]
  n_merge_acc <- if ("total_accepted_merges" %in% names(fit)) fit$total_accepted_merges[best_idx] else 0L
  n_split_acc <- if ("total_accepted_splits" %in% names(fit)) fit$total_accepted_splits[best_idx] else 0L
  
  # Result row
  res <- tibble(
    scenario           = row$scenario,
    method             = row$method,
    clusterwise        = row$clusterwise,
    n_sweet            = row$n_sweet,
    n_ms               = row$n_ms,
    n_merge            = row$n_merge,
    n_split            = row$n_split,
    nu_0               = row$nu_0,
    seed               = row$seed,
    final_D            = best_D,
    ARI                = ari,
    ARI_Binder         = ari_binder,
    ARI_VI             = ari_vi,
    n_accepted_sweets  = n_sweet_acc,
    n_accepted_merges  = n_merge_acc,
    n_accepted_splits  = n_split_acc,
    final_clust        = list(best_c),
    .rows = 1L
  )
  
  # Save immediately
  fname <- sprintf(
    "sim_%s_wass_clustw%s_nsweet%d_nms%d_merge%d_split%d_nu0%d.rds",
    row$scenario,
    ifelse(row$clusterwise, "T", "F"),
    row$n_sweet,
    row$n_ms,
    row$n_merge,
    row$n_split,
    row$nu_0
  )
  
  saveRDS(res, file.path(out_dir, fname))
  
  invisible(res)
}

# ---- 9. Run in parallel ----
cl <- makeCluster(detectCores() - 1)   # leave one core free

clusterExport(cl, varlist = c(
  "run_grid", "prior_list",
  "scen1_data", "scen1_post", "theta_scen1",
  "scen2_data", "scen2_post", "theta_scen2",
  "run_one_config", "out_dir"
))

clusterEvalQ(cl, {
  library(AntMAN)
  library(mclust)
  library(salso)
  library(tidyverse)
  source("src/search_algorithm.R")
  source("src/mixture_utils.R")
  source("src/distances.R")
})

pboptions(type = "timer")
cat("Starting parallel run (", nrow(run_grid), " configurations) ...\n")

system.time({
  res_list <- pblapply(
    X = seq_len(nrow(run_grid)),
    FUN = run_one_config,
    run_grid    = run_grid,
    prior_list  = prior_list,
    theta_scen1 = theta_scen1,
    theta_scen2 = theta_scen2,
    cl = cl
  )
})

stopCluster(cl)

cat("All configurations processed and saved individually.\n")
cat("Results are in: ", out_dir, "\n")
cat("Posterior objects saved as: scen1_post.rds and scen2_post.rds\n")