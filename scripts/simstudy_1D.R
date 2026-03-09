# =============================================================================
# Simulation study – 1D Gaussian mixtures – evaluating search-based clustering
# =============================================================================
# Purpose:   Compare performance of different discrepancy measures + refinement strategies
#            on two synthetic 1D scenarios using posterior-informed initializations
# Date:      March 2026
# Author:    Alex
#
# Dependencies:
# - AntMAN
# - mclust (for adjustedRandIndex)
# - salso (for Binder & VI)
# - tidyverse, pbapply, parallel
# - Custom source files

# ---- 1. Clear environment & load libraries ----
rm(list = ls())
library(AntMAN)
library(mclust)         # adjustedRandIndex
library(salso)          # binder, VI
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
  N <- 100
  K <- 3
  mu_true    <- c(-2, 0, 1)
  sd_true    <- c(1, 2, 1) * 0.2
  weights_true <- c(0.3, 0.5, 0.2)
  cluster_true <- rep(1:K, round(weights_true * N))
  y <- numeric(N)
  for (k in 1:K) {
    idx <- cluster_true == k
    ni  <- sum(idx)
    y[idx] <- rnorm(ni, mean = mu_true[k], sd = sd_true[k])
  }
  list(y = y, cluster_true = cluster_true, K = K)
}

generate_scenario2 <- function() {
  set.seed(2)
  N <- 150
  K <- 4
  mu_true    <- c(-4, -1, 1, 4)
  sd_true    <- rep(0.5, 4)
  weights_true <- rep(0.25, 4)
  cluster_true <- rep(1:K, round(weights_true * N))
  y <- numeric(N)
  for (k in 1:K) {
    idx <- cluster_true == k
    ni  <- sum(idx)
    y[idx] <- rnorm(ni, mean = mu_true[k], sd = sd_true[k])
  }
  list(y = y, cluster_true = cluster_true, K = K)
}

# ---- 4. AntMAN posterior draws ----
get_antman_posterior <- function(y, K_true) {
  alpha_0 <- 1
  kappa_0 <- 0.01
  nu_0    <- 3
  sig2_0  <- 2
  mu_0    <- 0
  
  n_iter   <- 20000
  burn_in  <- n_iter / 4
  thin     <- 10
  
  mixture_uvn_params <- AM_mix_hyperparams_uninorm(
    m0 = mu_0, k0 = kappa_0, nu0 = nu_0, sig02 = sig2_0
  )
  
  mcmc_params <- AM_mcmc_parameters(
    niter = n_iter, burnin = burn_in, thin = thin, verbose = 0
  )
  
  components_prior <- AM_mix_components_prior_pois(Lambda = K_true)
  
  mix_post_draws <- AM_mcmc_fit(
    y = y,
    mix_kernel_hyperparams = mixture_uvn_params,
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
  mu_0    = 0,
  a_0     = 2,
  b_0     = 5
)

# ---- 6. Parameter grid ----
methods     <- c("KS", "Pearson", "Wasserstein")
clusterwise <- c(FALSE, TRUE)
versions    <- c("grid", "quantile")

run_grid <- expand.grid(
  method      = methods,
  clusterwise = clusterwise,
  version     = versions,
  n_sweet     = c(1, 5),
  n_ms        = c(0, 5),
  n_merge     = 5,
  n_split     = 5,
  alpha_0     = c(2, 5),
  scenario    = c("scen1", "scen2"),
  seed        = 1,
  stringsAsFactors = FALSE
) |>
  filter(
    !(method %in% c("KS", "Pearson") & version == "quantile")
  ) |>
  filter(
    !(method == "Wasserstein" & clusterwise == TRUE & version == "grid")
  ) %>% filter(
    !(method == "Wasserstein" & clusterwise == FALSE & version == "quantile")
  )

dim(run_grid)

# ---- 7. Parallel worker function ----
run_one_config <- function(i, run_grid, prior_list) {
  row <- run_grid[i, ]
  set.seed(row$seed)
  
  # Create output directory (safe to call multiple times)
  out_dir <- "results/simstudy/1D"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Select scenario
  if (row$scenario == "scen1") {
    dat <- generate_scenario1()
  } else {
    dat <- generate_scenario2()
  }
  y            <- dat$y
  cluster_true <- dat$cluster_true
  K_true       <- dat$K
  
  # AntMAN initialization
  post <- get_antman_posterior(y, K_true)
  
  # Run clustering
  fit <- tryCatch(
    run_clustering(
      data               = y,
      clustering_matrix  = post$clustering_matrix,
      posterior_params   = post$mix_post_draws,
      prior_list         = prior_list,
      method             = row$method,
      clusterwise        = row$clusterwise,
      estimator          = "mle",
      n_runs             = 2,
      n_sweet            = row$n_sweet,
      tol                = 1e-10,
      n_ms               = row$n_ms,
      n_merge            = row$n_merge,
      n_split            = row$n_split,
      version            = row$version,
      M                  = 2000,
      alpha_0            = row$alpha_0
    ),
    error = function(e) {
      # Minimal error reporting to console (will appear in pblapply output)
      message(sprintf("Config %d failed: %s", i, e$message))
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
  
  # Posterior-based summaries
  binder_clust <- salso(post$clustering_matrix, loss = binder())
  vi_clust     <- salso(post$clustering_matrix, loss = VI())
  ari_binder   <- adjustedRandIndex(binder_clust, cluster_true)
  ari_vi       <- adjustedRandIndex(vi_clust, cluster_true)
  
  # Counts
  n_sweet_acc <- fit$total_accepted_sweets[best_idx]
  n_merge_acc <- if ("total_accepted_merges" %in% names(fit)) fit$total_accepted_merges[best_idx] else 0L
  n_split_acc <- if ("total_accepted_splits" %in% names(fit)) fit$total_accepted_splits[best_idx] else 0L
  
  # Create result tibble
  res <- tibble(
    scenario          = row$scenario,
    method            = row$method,
    clusterwise       = row$clusterwise,
    version           = row$version,
    n_sweet           = row$n_sweet,
    n_ms              = row$n_ms,
    n_merge           = row$n_merge,
    n_split           = row$n_split,
    alpha_0           = row$alpha_0,
    seed              = row$seed,
    final_D           = best_D,
    ARI               = ari,
    ARI_Binder        = ari_binder,
    ARI_VI            = ari_vi,
    n_accepted_sweets = n_sweet_acc,
    n_accepted_merges = n_merge_acc,
    n_accepted_splits = n_split_acc,
    final_clust       = list(best_c),
    .rows             = 1L
  )
  
  # Save immediately – descriptive filename
  fname <- sprintf(
    "sim_%s_%s_cw%s_v%s_ns%d_a0%d.rds",
    row$scenario,
    row$method,
    ifelse(row$clusterwise, "T", "F"),
    row$version,
    row$n_sweet,
    row$alpha_0
  )
  
  saveRDS(res, file.path(out_dir, fname))
  
}

# ---- 8. Run in parallel ----
cl <- makeCluster(detectCores())
clusterExport(cl, c("run_grid", "prior_list", "generate_scenario1", "generate_scenario2",
                    "get_antman_posterior", "run_clustering", "run_one_config"))
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
    X   = seq_len(nrow(run_grid)),
    FUN = run_one_config,
    run_grid   = run_grid,
    prior_list = prior_list,
    cl = cl
  )
})

stopCluster(cl)

cat("All configurations processed and saved individually.\n")
cat("Files are in: results/simstudy/1D/\n")
cat("Number of successful runs:", sum(!sapply(res_list, is.null)), "\n")