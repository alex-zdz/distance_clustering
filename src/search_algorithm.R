########################################################################
# Functions directly associated with the Search Algorithm
########################################################################

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
                                  method = "KS",
                                  clusterwise = FALSE,
                                  estimator = "posterior",
                                  version = "grid",          # new – only relevant for Wasserstein
                                  M = 2000,
                                  alpha_0 = 2) {                # new – only relevant for Wasserstein
  
  estimator <- match.arg(estimator, choices = c("posterior", "mle"))
  n_candidates <- nrow(clustering_matrix)
  distance <- numeric(n_candidates)
  
  for (i in seq_len(n_candidates)) {
    c_i <- clustering_matrix[i, ]
    
    if (!is.null(posterior_params) && estimator == "posterior") {
      w_i   <- posterior_params$weight_post[[i]]
      mu_i  <- posterior_params$mu_post[[i]]
      var_i <- posterior_params$var_post[[i]]
    } else {
      param <- get_mixture_params(
        clustering    = c_i,
        data       = data,
        estimator  = estimator,
        prior_list = prior_list,
        alpha_0    = alpha_0
      )
      w_i   <- param$weight
      mu_i  <- param$mean
      var_i <- param$var
    }
    
    distance[i] <- compute_distance(
      data       = data,
      weight     = w_i,
      mean       = mu_i,
      var        = var_i,
      method     = method,
      clusterwise = clusterwise,
      c_alloc    = if (clusterwise) c_i else NULL,
      version    = version,
      M          = M
    )
  }
  
  best_index <- which.min(distance)
  
  list(
    D_current = distance[best_index],
    c_current = clustering_matrix[best_index, ]
  )
}


# Sweetening Phase
sweetening <- function(c_current, data, prior_list, D_current,
                       n_sweet = 100, tol = 1e-10,
                       method = "KS", clusterwise = FALSE,
                       estimator = c("posterior", "mle"),
                       version = "grid", M = 2000, alpha_0 = 2) {
  
  estimator <- match.arg(estimator)
  run_sweet <- 0
  
  while (run_sweet < n_sweet) {
    run_sweet <- run_sweet + 1
    random_perm <- sample(seq_along(c_current))
    
    for (i in random_perm) {
      c_minus <- fulfill_gap_label(c_current[-i])
      K_sweet <- max(c_minus) + 1
      losses <- numeric(K_sweet)
      
      for (k in seq_len(K_sweet)) {
        c_candidate <- append(c_minus, k, after = i - 1)
        
        param <- get_mixture_params(
          c_candidate, data,
          estimator = estimator,
          prior_list = prior_list,
          alpha_0 = alpha_0
        )
        
        losses[k] <- compute_distance(
          data       = data,
          weight     = param$weight,
          mean       = param$mean,
          var        = param$var,
          method     = method,
          clusterwise = clusterwise,
          c_alloc    = if (clusterwise) c_candidate else NULL,
          version    = version,
          M          = M
        )
      }
      
      c_current[i] <- which.min(losses)
      c_current <- fulfill_gap_label(c_current)
    }
    
    D_new <- min(losses)
    
    if (abs(D_new - D_current) < tol) {
      D_current <- D_new
      break
    }
    
    D_current <- D_new
  }
  
  list(c_current = c_current, D_current = D_current, n_sweet = run_sweet)
}


# Merge and Split Phase
merge_split_phase <- function(c_current,
                              D_current,
                              data,
                              prior_list,
                              n_ms = 1,
                              n_merge = 1,
                              n_split = 1,
                              method = "KS",
                              clusterwise = FALSE,
                              estimator = c("posterior", "mle"),
                              version = "grid",
                              M = 2000) {
  
  estimator <- match.arg(estimator)
  n_cluster <- max(c_current)
  n_merge_accept <- 0
  n_split_accept <- 0
  
  for (iter in seq_len(n_ms)) {
    
    ## --- Merge ---
    if (n_cluster > 1) {
      pairs <- t(combn(1:n_cluster, 2))
      merge_pairs <- pairs[sample(nrow(pairs), min(n_merge, nrow(pairs))), , drop = FALSE]
      
      for (j in seq_len(nrow(merge_pairs))) {
        c_merge <- c_current
        c_merge[c_merge == merge_pairs[j, 1]] <- merge_pairs[j, 2]
        c_merge <- fulfill_gap_label(c_merge)
        
        param <- get_mixture_params(c_merge, data, estimator, prior_list)
        
        D_new <- compute_distance(
          data       = data,
          weight     = param$weight,
          mean       = param$mean,
          var        = param$var,
          method     = method,
          clusterwise = clusterwise,
          c_alloc    = if (clusterwise) c_merge else NULL,
          version    = version,
          M          = M
        )
        
        if (D_new < D_current) {
          c_current <- c_merge
          n_cluster <- max(c_current)
          D_current <- D_new
          n_merge_accept <- n_merge_accept + 1
        }
      }
    }
    
    ## --- Split ---
    for (j in seq_len(min(n_split, n_cluster))) {
      cl_to_split <- sample(seq_len(n_cluster), 1)
      idx_split <- which(c_current == cl_to_split)
      c_split <- c_current
      
      for (i in idx_split)
        if (runif(1) < 0.5) c_split[i] <- n_cluster + 1
      
      c_split <- fulfill_gap_label(c_split)
      
      param <- get_mixture_params(c_split, data, estimator, prior_list)
      D_new <- compute_distance(
        data       = data,
        weight     = param$weight,
        mean       = param$mean,
        var        = param$var,
        method     = method,
        clusterwise = clusterwise,
        c_alloc    = if (clusterwise) c_split else NULL,
        version    = version,
        M          = M
      )
      
      if (D_new < D_current) {
        c_current <- c_split
        n_cluster <- max(c_current)
        D_current <- D_new
        n_split_accept <- n_split_accept + 1
      }
    }
  }
  
  list(
    c_current = c_current,
    D_current = D_current,
    n_cluster = n_cluster,
    n_merge_accept = n_merge_accept,
    n_split_accept = n_split_accept
  )
}


# Entire clustering algorithm
run_clustering <- function(data,
                           clustering_matrix,
                           posterior_draws,
                           prior_list,
                           method      = "KS",
                           clusterwise = FALSE,
                           estimator   = c("posterior", "mle"),
                           n_runs      = 10,
                           n_sweet     = 100,
                           tol         = 1e-10,
                           n_ms        = 1,
                           n_merge     = 1,
                           n_split     = 1,
                           version     = "grid",    # new
                           M           = 2000,
                           alpha_0     = 2
                           ) {    # new
  
  estimator <- match.arg(estimator)
  
  posterior_params <- if (estimator == "posterior")
    prepare_AM_posterior_params(posterior_draws) else NULL
  
  c_record               <- matrix(NA, n_runs, ncol(clustering_matrix))
  distance_record        <- numeric(n_runs)
  total_accepted_merges  <- integer(n_runs)
  total_accepted_splits  <- integer(n_runs)
  total_accepted_sweets  <- integer(n_runs)
  
  init <- initialize_clustering(
    data, clustering_matrix,
    posterior_params, prior_list,
    method, clusterwise, estimator,
    version = version, M = M, alpha_0 = alpha_0
  )
  
  c_current <- init$c_current
  D_current <- init$D_current
  
  for (run in seq_len(n_runs)) {
    
    sweet <- sweetening(
      c_current, data, prior_list, D_current,
      n_sweet, tol, method, clusterwise, estimator,
      version = version, M = M, alpha_0 = alpha_0
    )
    
    ms <- merge_split_phase(
      sweet$c_current, sweet$D_current,
      data, prior_list,
      n_ms, n_merge, n_split,
      method, clusterwise, estimator,
      version = version, M = M
    )
    
    c_current <- ms$c_current
    D_current <- ms$D_current
    
    c_record[run, ]             <- c_current 
    distance_record[run]        <- D_current
    total_accepted_merges[run]  <- ms$n_merge_accept
    total_accepted_splits[run]  <- ms$n_split_accept
    total_accepted_sweets[run]  <- sweet$n_sweet
    
  }
  
  final_clustering <- c_record[which.min(distance_record), ] 
  
  list(
    final_clustering           = final_clustering,
    c_record                   = c_record,
    distance_record            = distance_record,
    total_accepted_merges      = total_accepted_merges,
    total_accepted_splits      = total_accepted_splits,
    total_accepted_sweets      = total_accepted_sweets
  )
}

