
generate_mixture_data <- function(
    N,
    K = 3,
    dim = 1,
    alpha = 1,
    mu_true = NULL,      # vector for dim = 1, matrix (K x dim) for dim > 1
    Sigma_true = NULL         # vector for dim = 1, list of covariance matrices for dim > 1
) {
  #--- Component weights from Dirichlet(alpha,...,alpha) ---
  gamma_draws <- rgamma(K, shape = alpha, rate = 1)
  weights_true <- gamma_draws / sum(gamma_draws)
  
  #--- Assign clusters ---
  cluster_true <- sample(1:K, N, replace = TRUE, prob = weights_true)
  
  #--- Default mu_true ---
  if (dim == 1) {
    if (is.null(mu_true)) mu_true <- seq(-K, K, length.out = K)
  } else {
    if (is.null(mu_true)) {
      mu_true <- matrix(0, nrow = K, ncol = dim)
      for (i in 1:K) mu_true[i, ] <- rep(i, dim)  # simple separated defaults, consider making random
    }
  }
  
  #--- Default covariances ---
  if (dim == 1) {
    if (is.null(Sigma_true)) Sigma_true <- rep(1, K)
  } else {
    if (is.null(Sigma_true)) {
      Sigma_true <- vector("list", K)
      for (i in 1:K) Sigma_true[[i]] <- diag(dim)      # identity covariance
    }
  }
  
  #--- Initialize data ---
  y <- matrix(NA, nrow = N, ncol = dim)
  
  #--- Generate Gaussian data ---
  for (i in 1:K) {
    idx <- cluster_true == i
    ni  <- sum(idx)
    
    if (ni > 0) {
      if (dim == 1) {
        y[idx, ] <- rnorm(ni, mean = mu_true[i], sd = Sigma_true[i])
      } else {
        y[idx, ] <- MASS::mvrnorm(n = ni, mu = mu_true[i, ], Sigma = Sigma_true[[i]])
      }
    }
  }
  
  #--- Return vector for 1d case ---
  if (dim == 1) y <- as.numeric(y)
  
  return(list(
    data = y,
    cluster_true = cluster_true,
    weights_true = weights_true,
    mu_true = mu_true,
    Sigma_true = Sigma_true
  ))
}



#Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=sqrt(s)) )
mixture_cdf_1D <- function(x, weight, mean, sig2) sapply(x, function(xi) sum(weight * pnorm(xi, mean = mean, sd = sqrt(sig2))))

mixture_pdf_1D <- function(t, mean, sig2, weight) {
  mix_pdf <-   sapply(t, function(xi) sum(weight * dnorm(xi, mean = mean, sd = sqrt(sig2))))
  mix_pdf / sum(mix_pdf)
}

emp_pdf <- function(t, proj) {
  # Create histogram with length(t) bins over range of t
  h <- hist(proj,
            breaks = seq(min(t), max(t), length.out = length(t) + 1),
            plot = FALSE)
  
  density <- h$counts
  # Normalize
  density / sum(density)
}

KS_mixture_distance <- function(data, weight, mean, sig2) {
  as.numeric(ks.test(data, mixture_cdf_1D, weight, mean, sig2)$statistic)
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
  cdf_vals <- mixture_cdf_1D(bins, weight, mean, sig2)  # vectorized
  ff <- numeric(n_bins)
  ff[1] <- cdf_vals[1]
  ff[2:(n_bins-1)] <- diff(cdf_vals)
  ff[n_bins] <- 1 - cdf_vals[n_bins-1]
  
  # avoid zero division
  ff[ff == 0] <- 1e-50
  
  # Pearson χ² statistic
  sum((gg - ff)^2 / ff)
}



# Compute Distance Between Data and Mixture

compute_distance <- function(data, weight, mean, Sigma, method = "KS") {
  
  if (is.null(dim(data))) {
    D <- 1
  } else {
    D <- ncol(data)
  }
  
  if(D != 1 & method != "Wasserstein"){
    stop("only Sliced Wasserstein for higher dimensions")
  }
  
  K = length(weight)
  
  if (method == "KS") {
    KS_mixture_distance(data, weight, mean, Sigma)
  }else if (method == "Wasserstein") {
    # for now always p = 2, consider adapting
    M = 2e3
    t <- seq(
      from = -max(abs(y)) * sqrt(2 * D),
      to   =  max(abs(y)) * sqrt(2 * D),
      length.out = M
    )
    
    if(D == 1){
      
      I0 = mixture_pdf(t, mean, Sigma, weight)
      I1 = emp_pdf(t, data)
      Wasserstein_1D(I0, I1, p = 2)
      
    }else{ # Sliced Wasserstein
      
      L = 20 # for now fixed to low number
      
      theta = generateTheta(L, D)
      yproj = y %*% t(theta)
      
      projected_Sigma <- matrix(0, nrow = K, ncol = L)
      projected_Mu    <- matrix(0, nrow = K, ncol = L)
      
      for (k in 1:K) {
        for (l in 1:L) {
          projected_Sigma[k, l] <- sqrt( theta[l, , drop = FALSE] %*% Sigma[[k]] %*% t(theta[l, , drop = FALSE]) )
          projected_Mu[k, l]    <- theta[l, , drop = FALSE] %*% mean[[k]]
        }
      }
      
      SW = 0
      
      for (l in 1:L){
        RIx = mixture_pdf(t, projectedMu[,l], projectedSigma[,l], weights)
        RIy = emp_pdf(t, yproj[,l])
        SW = SW + pWasserstein(RIx, RIy,p=2) / L
      }
      
      SW
      
    }
    
  } else if (method == "Pearson") {
    Pearson_mixture_distance(data, weight, mean, Sigma)
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
  
  n <- nrow(data)
  D <- ncol(data)
  K <- max(c_alloc) # Number of clusters (assumed labeled 1...K)
  
  # Allocate storage
  alpha_post <- numeric(K)
  mu_post <- vector("list", K)
  Sigma_post <- vector("list", K)   # D>1: covariance matrices. D=1: sigma2.
  
  # ---------------------------------------------------------------------
  # CASE 1: UNIVARIATE (D = 1): Normal–Inverse-Gamma
  # ---------------------------------------------------------------------
  if (D == 1) {
    
    a_0     <- prior_list$a_0
    b_0     <- prior_list$b_0
    
    for (k in seq_len(K)) {
      y_k <- data[c_alloc == k]
      n_k <- length(y_k)
      
      # Posterior for weights
      alpha_post[k] <- alpha_0 + n_k
      
      ybar <- mean(y_k)
      S <- sum((y_k - ybar)^2)    # within-cluster sum of squares
      
      # posterior hyperparams
      kappa_post <- kappa_0 + n_k
      mu_n    <- (kappa_0 * mu_0 + n_k * ybar) / kappa_post
      a_post     <- a_0 + 0.5 * n_k
      b_post     <- b_0 + 0.5 * S + (kappa_0 * n_k) / (2 * kappa_post) * (ybar - mu_0)^2
      
      # sample variance then mean
      Sigma_post[[k]] <- 1 / rgamma(1, shape = a_post, rate = b_post)    # Inv-Gamma sample
      mu_post[[k]] <- rnorm(1, mean = mu_n, sd = sqrt(Sigma_post[[k]] / kappa_post))
      
    }
  
    # ---------------------------------------------------------------------
    # CASE 2: MULTIVARIATE (D > 1): Normal–Inverse-Wishart
    # ---------------------------------------------------------------------
    } else {
      
      nu_0 <- prior_list$nu_0
      S_0  <- prior_list$S_0
      
      for (k in seq_len(K)) {
        y_k <- data[c_alloc == k, , drop = FALSE]
        n_k <- nrow(y_k)
        
        # Posterior for weights
        alpha_post[k] <- alpha_0[k] + n_k
        
        # Means
        y_bar <- colMeans(y_k)
        diff <- sweep(y_k, 2, y_bar)
        S_k_data <- t(diff) %*% diff
        
        # Posterior mean of the mean
        kappa_post <- kappa_0 + n_k
        mu_post_k <- (kappa_0 * mu_0 + n_k * y_bar) / kappa_post
        
        nu_post <- nu_0 + n_k
        cross_term <- (kappa_0 * n_k / kappa_post) * tcrossprod(y_bar - mu_0) 
        S_k_post <- S_0 + S_data + cross_term
        
        # sample Sigma (Inv-Wishart) then mu | Sigma
        Sigma_post[[k]] <- MCMCpack::riwish(nu_post, S_k_post)
        mu_post[[k]] <- MASS::mvrnorm(1, mu_post_k, Sigma / kappa_post)
      }
    } 
    
  weights_post <- as.numeric(MCMCpack::rdirichlet(1, alpha_post))
  
  list(
    weights_post   = weights_post,
    mu_post   = mu_post,
    Sigma_post = Sigma_post
  )
}



#####################
# Sliced Wasserstein
#####################

generateTheta <- function(L, d) {
  theta <- matrix(0, nrow = L, ncol = d)
  
  # First vector: uniform on sphere
  th_l <- runif(d)
  th_l <- th_l / sqrt(sum(th_l^2))
  theta[1, ] <- th_l
  
  # Remaining vectors
  for (i in 2:L) {
    th_l <- rnorm(d)
    th_l <- th_l / sqrt(sum(th_l^2))
    
    m <- max(abs(theta[1:(i - 1), ] %*% th_l))
    
    while (m > 0.97) {
      th_l <- rnorm(d)
      th_l <- th_l / sqrt(sum(th_l^2))
      m <- max(abs(theta[1:(i - 1), ] %*% th_l))
    }
    theta[i, ] <- th_l
  }
  theta
}

Wasserstein_1D <- function(I0, I1, p = 2) {
  stopifnot(length(I0) == length(I1))
  
  eps <- 1e-7
  
  # Ensure strict positivity
  I0 <- I0 + eps
  I1 <- I1 + eps
  
  # Normalize to sum to 1 - both pdfs have already been normalized, so this step is redundant here
  I0 <- I0 / sum(I0)
  I1 <- I1 / sum(I1)
  
  # Compute CDFs
  J0 <- cumsum(I0)
  J1 <- cumsum(I1)
  
  # Grid
  x       <- seq_along(I0) - 1
  xtilde  <- seq(0, 1, length.out = length(I0))
  
  # Inverse CDF sampling (pseudo-quantiles)
  XI0 <- approx(J0, x, xout = xtilde, rule = 2)$y
  XI1 <- approx(J1, x, xout = xtilde, rule = 2)$y
  
  # Displacement field u(x)
  u <- approx(XI0, XI0 - XI1, xout = x, rule = 2)$y
  
  # # Transport map f(x) = x - u(x)
  # f <- x - u
  # 
  # # Potential phi(x)
  # phi <- cumsum(u / length(I0))
  # phi <- phi - mean(phi)
  
  # p-Wasserstein distance
  Wp <- (mean((abs(u)^p) * I0))^(1/p)
  
  Wp
}




########
# Search Algorithm
########



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
