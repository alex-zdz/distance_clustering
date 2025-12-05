#' Generate Mixture Model Data
#'
#' Generates synthetic data from a mixture model (Gaussian or skew-Gaussian).
#' Component weights are sampled from a Dirichlet distribution via normalized gamma draws.
#'
#' @param N Integer. Number of observations to generate.
#' @param K Integer. Number of mixture components. Default = 3.
#' @param type Character. Type of mixture distribution:
#'   - `"Gaussian"` (default)
#'   - `"Skew Gaussian"`
#'   - `"Skew-symmetric"` (not yet implemented).
#' @param alpha Numeric. Concentration parameter for the Dirichlet prior on weights.
#' @param means Numeric vector of length \code{K}. Component means. 
#'   Defaults to \code{seq(-2, 2, length.out = K)} if not provided.
#' @param sds Numeric vector of length \code{K}. Component standard deviations. 
#'   Defaults to \code{rep(1, K)} if not provided.
#' @param skews Numeric vector of length \code{K}. Component skewness parameters (only used for skew-Gaussian).
#'   Defaults to \code{rep(0, K)} if not provided.
#'
#' @return A list with three elements:
#'   \item{data}{Numeric vector of simulated values.}
#'   \item{cluster_true}{Integer vector of true cluster allocations.}
#'   \item{weights_true}{Numeric vector of true component weights.}
#'
#' @examples
#' set.seed(123)
#' sim <- generate_mixture_data(N = 200, K = 3, type = "Gaussian", alpha = 5)
#' hist(sim$data, breaks = 30)
#' table(sim$cluster_true)
#'
#' @export
generate_mixture_data <- function(
    N,
    K = 3,
    type = "Gaussian",
    alpha = 1,
    means = NULL,
    sds = NULL,
    skews = NULL
) {
  #--- Generate component weights from Dirichlet ---
  gamma_draws <- rgamma(K, shape = alpha, rate = 1)
  weights_true <- gamma_draws / sum(gamma_draws)
  
  #--- Assign clusters ---
  cluster_true <- sample(1:K, N, replace = TRUE, prob = weights_true)
  
  #--- Set default parameters if not provided ---
  if (is.null(means)) means <- seq(-K, K, length.out = K)
  if (is.null(sds)) sds <- rep(1, K)
  if (type %in% c("Skew Gaussian", "Skew-symmetric") & is.null(skews)) skews <- rep(0, K)
  
  #--- Initialize data vector ---
  y <- numeric(N)
  
  #--- Generate data based on type ---
  if (type == "Gaussian") {
    for (i in 1:K) {
      idx <- cluster_true == i
      y[idx] <- rnorm(sum(idx), mean = means[i], sd = sds[i])
    }
    
  } else if (type == "Skew Gaussian") {
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Package 'sn' is required for skew-normal generation. Install it via install.packages('sn').")
    }
    for (i in 1:K) {
      idx <- cluster_true == i
      y[idx] <- sn::rsn(sum(idx), xi = means[i], omega = sds[i], alpha = skews[i])
    }
    
  } else if (type == "Skew-symmetric") {
    stop("Skew-symmetric mixture not implemented in general form yet.")
    
  } else {
    stop("Unknown mixture type. Choose 'Gaussian', 'Skew Gaussian', or 'Skew-symmetric'.")
  }
  
  return(list(data = y, cluster_true = cluster_true, weights_true = weights_true))
}



#Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=sqrt(s)) )
mixture_cdf <- function(x, weight, mean, sig2) sapply(x, function(xi) sum(weight * pnorm(xi, mean = mean, sd = sqrt(sig2))))

KS_mixture_distance <- function(data, weight, mean, sig2) {
  as.numeric(ks.test(data, mixture_cdf, weight, mean, sig2)$statistic)
}

# Pearson distance new:
Pearson_mixture_distance <- function(data, weight, mean, sig2) {
  n <- length(data)
  Fn <- ecdf(data)
  bins <- seq(-50, 50, by = 0.5)
  n_bins <- length(bins) + 1
  
  # empirical bin probabilities
  bin_counts <- tabulate(cut(knots(Fn), breaks = c(-Inf, bins, Inf)), nbins = n_bins)
  gg <- bin_counts / n
  
  # mixture bin probabilities
  cdf_vals <- mixture_cdf(bins, weight, mean, sig2)  # vectorized
  ff <- numeric(n_bins)
  ff[1] <- cdf_vals[1]
  ff[2:(n_bins-1)] <- diff(cdf_vals)
  ff[n_bins] <- 1 - cdf_vals[n_bins-1]
  
  # avoid zero division
  ff[ff == 0] <- 1e-50
  
  # Pearson χ² statistic
  sum((gg - ff)^2 / ff)
}



#' Compute Distance Between Data and Mixture
#'
#' Wrapper function for computing distances between observed data and a Gaussian mixture.
#' Currently supports the Kolmogorov–Smirnov (KS) distance.
#'
#' @param data Numeric vector. Observed data.
#' @param w Numeric vector of component weights.
#' @param mean Numeric vector of component means.
#' @param sig2 Numeric vector of component variances.
#' @param method Character. Distance method to use. Currently only `"KS"` is implemented.
#'
#' @return Numeric scalar. Distance value.
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100)
#' compute_distance(data, w = c(0.5, 0.5), mean = c(0, 2), sig2 = c(1, 1), method = "KS")
#'
#' @export
compute_distance <- function(data, weight, mean, sig2, method = "KS") {
  if (method == "KS") {
    KS_mixture_distance(data, weight, mean, sig2)
  }else if (method == "W2") {
    W2(data, weight, mean, sig2)
  } else if (method == "Pearson") {
    Pearson_mixture_distance(data, weight, mean, sig2)
  } else {
    stop("Unknown distance method")
  }
}

# Prepare AntMAN output:
prepare_AM_posterior_params <- function(mix_post_draws) {
  list(
    weight = lapply(mix_post_draws$W, unlist),
    mu_post = lapply(mix_post_draws$mu, unlist),
    sigma2_post = lapply(mix_post_draws$sig2, unlist)
  )
}

# Calculate the posterior expectations of a mixture model conditional on the chosen clustering
mixture_posterior <- function(c_alloc, data, prior_list) {
  # Extract priors
  alpha_0 <- prior_list$alpha_0
  kappa_0 <- prior_list$kappa_0
  mu_0    <- prior_list$mu_0
  a_0     <- prior_list$a_0
  b_0     <- prior_list$b_0
  
  # Number of clusters (assumed labeled 1...K)
  K <- max(c_alloc)
  
  # Preallocate
  alpha_post  <- numeric(K)
  mu_post     <- numeric(K)
  sigma2_post <- numeric(K)
  
  for (k in seq_len(K)) {
    y_k <- data[c_alloc == k]
    n_k <- length(y_k)
    
    # Posterior for weights
    alpha_post[k] <- alpha_0 + n_k
    
    # Posterior for mean
    kappa_post <- kappa_0 + n_k
    mu_post[k] <- (kappa_0 * mu_0 + sum(y_k)) / kappa_post
    
    # Posterior for variance (# Old code was using mean(y_k) instead of mu_post[k])
    a_post <- a_0 + 0.5 * (n_k + 1)  # "+1" accounts for mean uncertainty
    b_post <- b_0 + 0.5 * sum((y_k - mu_post[k])^2) +
      0.5 * kappa_0 * (mu_post[k] - mu_0)^2
    
    sigma2_post[k] <- b_post / (a_post - 1)  # posterior mean of Inv-Gamma
  }
  
  weights_post <- alpha_post / sum(alpha_post)  # expected Dirichlet
  
  list(
    weights_post   = weights_post,
    mu_post   = mu_post,
    sigma2_post = sigma2_post
  )
}

fulfill_gap_label <- function(c_vec) {
  unique_labels <- sort(unique(c_vec))
  relabeled_vec <- match(c_vec, unique_labels)
  return(relabeled_vec)
}

# Initialization Phase of the Algorithm
initialize_clustering <- function(data,
                                  clustering_matrix,
                                  posterior_params = NULL,
                                  prior_list = NULL,
                                  method = "KS") {
  
  n_candidates <- nrow(clustering_matrix)
  distance <- numeric(n_candidates)
  
  for (i in seq_len(n_candidates)) {
    if (!is.null(posterior_params)) {
      # use given posterior parameters
      w_i    <- posterior_params$weight[[i]]
      mu_i   <- posterior_params$mu_post[[i]]
      sig2_i <- posterior_params$sigma2_post[[i]]
    } else {
      # compute posterior parameters from clustering
      c_i <- clustering_matrix[i, ]
      param_post <- mixture_posterior(c_i, data, prior_list)
      w_i    <- param_post$weights
      mu_i   <- param_post$mu_post
      sig2_i <- param_post$sigma2_post
    }
    
    distance[i] <- compute_distance(data, w_i, mu_i, sig2_i, method)
  }
  
  best_index <- which.min(distance)
  D_current <- distance[best_index]
  c_current <- clustering_matrix[best_index,]
  
  return(list(D_current = D_current, c_current = c_current))
}




#Sweetening Phase
sweetening <- function(c_current, data, prior_list, D_current,
                       n_sweet = 100, tol = 1e-10, method = "KS") {
  
  run_sweet <- 0
  
  while (run_sweet < n_sweet) {
    run_sweet <- run_sweet + 1
    
    # Loop over subjects in random order
    random_perm <- sample(seq_along(c_current))
    
    for (i in random_perm) {
      c_minus <- c_current[-i]
      K_sweet <- max(c_minus) + 1
      losses <- numeric(K_sweet)
      
      for (k in seq_len(K_sweet)) {
        c_candidate <- append(c_minus, k, after = i - 1)
        param_post <- mixture_posterior(c_candidate, data, prior_list)
        losses[k] <- compute_distance(data,
                                      weight = param_post$weights,
                                      mean = param_post$mu_post,
                                      sig2 = param_post$sigma2_post,
                                      method = method)
      }
      
      # Reassign i
      c_current[i] <- which.min(losses)
    }
    
    D_new <- min(losses)
    
    if (abs(D_new - D_current) < tol) {
      D_current <- D_new
      break
    }
    
    D_current <- D_new
  }
  
  return(list(c_current = c_current, D_current = D_current, n_sweet = n_sweet))
}

# Merge and Split Phase
merge_split_phase <- function(c_current,
                              D_current,
                              data,
                              prior_list,
                              n_ms   = 1,
                              n_merge = 1,
                              n_split = 1,
                              method  = "KS") {
  
  n_cluster <- max(c_current)
  n_merge_accept <- 0
  n_split_accept <- 0
  
  if (n_ms > 0) {
    for (iter in seq_len(n_ms)) {
      
      # --- Merging step ---
      if (n_cluster > 1) {
        pairs <- t(combn(1:n_cluster, 2))
        merge_pairs <- pairs[sample(nrow(pairs), min(n_merge, nrow(pairs))), , drop = FALSE]
        
        for (j in seq_len(nrow(merge_pairs))) {
          clusters_to_merge <- merge_pairs[j, ]
          c_merge <- c_current
          c_merge[c_merge == clusters_to_merge[1]] <- clusters_to_merge[2]
          c_merge <- fulfill_gap_label(c_merge)
          
          param_post <- mixture_posterior(c_merge, data, prior_list)
          D_new <- compute_distance(data,
                                    weight    = param_post$weights,
                                    mean = param_post$mu_post,
                                    sig2 = param_post$sigma2_post,
                                    method = method)
          if (D_new < D_current) {
            c_current <- c_merge
            n_cluster <- max(c_current)
            D_current <- D_new
            n_merge_accept <- n_merge_accept + 1
          }
        }
      }
      
      # --- Splitting step ---
      c_split <- c_current
      
      for (j in seq_len(min(n_split, n_cluster))) {
        cl_to_split <- sample(1:n_cluster, 1)
        idx_split <- which(c_current == cl_to_split)
        
        # random binary split
        for (i in idx_split) {
          if (runif(1) < 0.5) {
            c_split[i] <- n_cluster + 1
          }
        }
        c_split <- fulfill_gap_label(c_split)
        
        param_post <- mixture_posterior(c_split, data, prior_list)
        D_new <- compute_distance(data,
                                  weight    = param_post$weights,
                                  mean = param_post$mu_post,
                                  sig2 = param_post$sigma2_post,
                                  method = method)
        if (D_new < D_current) {
          c_current <- c_split
          n_cluster <- max(c_current)
          D_current <- D_new
          n_split_accept <- n_split_accept + 1
        }
      }
    }
  }
  
  list(
    c_current      = c_current,
    D_current      = D_current,
    n_cluster      = n_cluster,
    n_merge_accept = n_merge_accept,
    n_split_accept = n_split_accept
  )
}

# Main clustering
run_clustering <- function(data,
                           clustering_matrix,
                           posterior_params,   
                           prior_list,
                           method         = "KS",
                           n_runs         = 10,
                           n_sweet = 100,
                           tol            = 1e-10,
                           n_ms          = 1,
                           n_merge        = 1,
                           n_split        = 1) {
  
  # storage (per-run)
  distance_record <- numeric(n_runs)           # stores current distance after each run
  total_accepted_sweets <- integer(n_runs)
  total_accepted_merges <- integer(n_runs)           # accepted merges per run
  total_accepted_splits <- integer(n_runs)           # accepted splits per run
  
  
  for (run in seq_len(n_runs)) {
    print(run)
    ## ---------------------------
    ## 1) Initialization phase
    ## ---------------------------
    init_res <- initialize_clustering(
      data              = data,
      clustering_matrix = clustering_matrix,   # list of candidate clusterings
      posterior_params  = posterior_params,
      method            = method
    )
    
    c_current <- init_res$c_current
    D_current  <- init_res$D_current
    
    ## ---------------------------
    ## 2) Sweetening phase
    ## ---------------------------
    sweet_res <- sweetening(
      c_current   = c_current,
      data        = data,
      prior_list  = prior_list,
      D_current   = D_current,
      n_sweet = n_sweet,
      tol         = tol,
      method      = method
    )
    c_current <- sweet_res$c_current
    D_current <- sweet_res$D_current
    total_accepted_sweets[run] <- sweet_res$n_sweet
    
    ## ---------------------------
    ## 3) Merge–Split phase
    ## ---------------------------
    ms_res <- merge_split_phase(
      c_current  = c_current,
      D_current  = D_current,
      data       = data,
      prior_list = prior_list,
      n_ms       = n_ms,
      n_merge    = n_merge,
      n_split    = n_split,
      method     = method
    )
    
    # update clustering and distance
    c_current <- ms_res$c_current
    D_current <- ms_res$D_current
    # record accepted merges/splits
    total_accepted_merges[run] <-  ms_res$n_merge_accept
    total_accepted_splits[run] <- ms_res$n_split_accept 
    
    # record distance
    distance_record[run] <- D_current
  }
  
  # return final clustering and diagnostics
  list(
    clustering = c_current,
    distance_record = distance_record,         
    total_accepted_merges = total_accepted_merges,
    total_accepted_splits = total_accepted_splits,
    total_accepted_sweets = total_accepted_sweets
  )
}
