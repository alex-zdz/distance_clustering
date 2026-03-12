# =============================================================================
# Simulation study – Gaussian mixtures in variable dimension D
# evaluating search-based clustering with sliced Wasserstein
# =============================================================================
# Purpose: Compare performance of Wasserstein discrepancy + refinement strategies
#          across different data dimensions using posterior-informed initializations
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

# ---- 3. Data generation functions (parameterized by D) ----
generate_scenario1 <- function(D) {
  set.seed(1)
  generate_mixture_data(
    N = 100,
    K = 3,
    dim = D,
    mu_true = matrix(seq(-2, 2, length.out = 3*D), nrow = 3, byrow = TRUE),
    Sigma_true = lapply(1:3, function(i) diag(D) * 0.2)
  )
}

generate_scenario2 <- function(D) {
  set.seed(2)
  generate_mixture_data(
    N = 150,
    K = 4,
    dim = D,
    mu_true = matrix(seq(-4, 4, length.out = 4*D), nrow = 4, byrow = TRUE),
    Sigma_true = lapply(1:4, function(i) diag(D) * 0.5)
  )
}

# ---- 4. AntMAN posterior draws (dimension-aware) ----
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

# ---- 5. Base prior list (will be adjusted per D and per config) ----
prior_list_base <- list(
  alpha_0 = 1,
  kappa_0 = 0.01,
  nu_0    = 4,
  S_0     = NULL   # will be set per D
)

# ---- 6. Parameter grid ----
methods     <- "Wasserstein"
clusterwise <- c(FALSE, TRUE)
dimensions  <- c(5, 10)   # ← customize here

run_grid <- expand.grid(
  method      = methods,
  clusterwise = clusterwise,
  n_sweet     = c(1, 3),
  n_ms        = c(3),
  n_merge     = 5,
  n_split     = 5,
  nu_0        = c(2, 5),
  D           = dimensions,
  scenario    = c("scen1", "scen2"),
  seed        = 1,
  stringsAsFactors = FALSE
)

cat("Total configurations:", nrow(run_grid), "\n")

# ---- 7. Precompute per dimension: data, posteriors, theta ----
cat("Precomputing data, posteriors and projection directions per dimension ...\n")

out_dir_base <- "results/simstudy/variableD"
if (!dir.exists(out_dir_base)) dir.create(out_dir_base, recursive = TRUE)

# Storage for parallel access
all_data  <- list()
all_post  <- list()
all_theta <- list()

for (d in unique(run_grid$D)) {
  cat("  Dimension D =", d, "...\n")
  
  out_dir_d <- file.path(out_dir_base, sprintf("D%02d", d))
  if (!dir.exists(out_dir_d)) dir.create(out_dir_d, recursive = TRUE)
  
  # Generate data
  scen1_data_d <- generate_scenario1(d)
  scen2_data_d <- generate_scenario2(d)
  
  saveRDS(scen1_data_d, file.path(out_dir_d, "scen1_data.rds"))
  saveRDS(scen2_data_d, file.path(out_dir_d, "scen2_data.rds"))
  
  # AntMAN posteriors
  scen1_post_d <- get_antman_posterior(scen1_data_d$data, length(unique(scen1_data_d$cluster_true)))
  scen2_post_d <- get_antman_posterior(scen2_data_d$data, length(unique(scen2_data_d$cluster_true)))
  
  saveRDS(scen1_post_d, file.path(out_dir_d, "scen1_post.rds"))
  saveRDS(scen2_post_d, file.path(out_dir_d, "scen2_post.rds"))
  
  # Fixed theta (same for all configs in this D)
  L_fixed <- 20
  set.seed(123 + d)  # reproducible but different per D
  theta_d <- generate_theta_normal(L = L_fixed, d = d)
  
  # Store for parallel workers
  all_data[[as.character(d)]]  <- list(scen1 = scen1_data_d, scen2 = scen2_data_d)
  all_post[[as.character(d)]]  <- list(scen1 = scen1_post_d, scen2 = scen2_post_d)
  all_theta[[as.character(d)]] <- theta_d
}

# ---- 8. Parallel worker function ----
run_one_config <- function(i, run_grid, prior_list_base) {
  row <- run_grid[i, ]
  set.seed(row$seed)
  
  d_str <- as.character(row$D)
  out_dir <- file.path("results/simstudy/variableD", sprintf("D%02d", row$D))
  
  # Select correct data/post/theta for this D
  if (row$scenario == "scen1") {
    scen_data <- all_data[[d_str]]$scen1
    post      <- all_post[[d_str]]$scen1
  } else {
    scen_data <- all_data[[d_str]]$scen2
    post      <- all_post[[d_str]]$scen2
  }
  
  theta <- all_theta[[d_str]]
  
  y <- scen_data$data
  cluster_true <- scen_data$cluster_true
  
  # Prepare prior list for this D and nu_0
  prior_list <- prior_list_base
  prior_list$mu_0 <- rep(0, row$D)
  prior_list$S_0  <- diag(row$D) * 5
  prior_list$nu_0 <- row$D + row$nu_0
  
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
      nu_0              = row$D + row$nu_0,
      L                 = nrow(theta),
      theta             = theta
    ),
    error = function(e) {
      message(sprintf("Config %d (D=%d, %s) failed: %s", i, row$D, row$scenario, e$message))
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
  
  # Posterior consensus ARI (for dashboard reference)
  binder_clust <- salso(post$clustering_matrix, loss = binder())
  vi_clust     <- salso(post$clustering_matrix, loss = VI())
  ari_binder   <- adjustedRandIndex(best_c, binder_clust)
  ari_vi       <- adjustedRandIndex(best_c, vi_clust)
  
  # Counts
  n_sweet_acc <- fit$total_accepted_sweets[best_idx]
  n_merge_acc <- if ("total_accepted_merges" %in% names(fit)) fit$total_accepted_merges[best_idx] else 0L
  n_split_acc <- if ("total_accepted_splits" %in% names(fit)) fit$total_accepted_splits[best_idx] else 0L
  
  # Result tibble
  res <- tibble(
    D                  = row$D,
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
  
  # Save
  fname <- sprintf(
    "sim_D%02d_%s_wass_clustw%s_nsweet%d_nms%d_merge%d_split%d_nu0%d.rds",
    row$D, row$scenario,
    ifelse(row$clusterwise, "T", "F"),
    row$n_sweet, row$n_ms, row$n_merge, row$n_split, row$nu_0
  )
  
  saveRDS(res, file.path(out_dir, fname))
  
  invisible(res)
}

# ---- 9. Run in parallel ----
cl <- makeCluster(detectCores() - 1)

clusterExport(cl, varlist = c(
  "run_grid", "prior_list_base",
  "all_data", "all_post", "all_theta",
  "run_one_config"
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
    run_grid       = run_grid,
    prior_list_base = prior_list_base,
    cl = cl
  )
})

stopCluster(cl)

cat("All configurations processed and saved.\n")
cat("Results are in subdirectories: results/simstudy/variableD/Dxx/\n")
cat("Posterior and data objects saved per D (scen1_post.rds, scen2_post.rds, etc.)\n")