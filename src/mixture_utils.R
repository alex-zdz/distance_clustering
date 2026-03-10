
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
      for (i in 1:K) Sigma_true[[i]] <- diag(dim) * 0.1     # identity covariance
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

# Prepare AntMAN output:
prepare_AM_posterior_params <- function(mix_post_draws) {
  
  D <- length(mix_post_draws$mu[[1]][[1]])
  
  weight_post = mix_post_draws$W
  
  if(D == 1){
    list(weight_post = weight_post,
         #weight = lapply(mix_post_draws$W, unlist),
         mu_post = lapply(mix_post_draws$mu, unlist),
         var_post = lapply(mix_post_draws$sig2, unlist)
    )
  }else{
    list(weight_post = weight_post,
         #weight = lapply(mix_post_draws$W, unlist),
         mu_post = mix_post_draws$mu,
         var_post = mix_post_draws$Sig
    )
  }
  
}


#Mixture_normal_cdf = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=sqrt(s)) )
mixture_cdf_1D <- function(x, weight, mean, sig2) sapply(x, function(xi) sum(weight * pnorm(xi, mean = mean, sd = sqrt(sig2))))

mixture_pdf_1D <- function(t, mean, sig2, weight) {
  mix_pdf <-   sapply(t, function(xi) sum(weight * dnorm(xi, mean = mean, sd = sqrt(sig2))))
  mix_pdf / sum(mix_pdf)
}


# Calculate the posterior expectations of a mixture model conditional on the chosen clustering
mixture_posterior <- function(c_alloc, data, prior_list) {
  
  use_post_mean <- TRUE
  
  # Extract priors
  alpha_0 <- prior_list$alpha_0
  kappa_0 <- prior_list$kappa_0
  mu_0    <- prior_list$mu_0
  
  n <- nrow(data)
  D <- if( is.null(ncol(data))) 1 else ncol(data)
  K <- max(c_alloc) # Number of clusters (assumed labeled 1...K)
  
  # ---------------------------------------------------------------------
  # CASE 1: UNIVARIATE (D = 1): Normal–Inverse-Gamma
  # ---------------------------------------------------------------------
  if (D == 1) {
    
    # Allocate storage
    alpha_post <- numeric(K)
    mu_post <- rep(NA, K)
    var_post <- rep(NA, K)   # D>1: covariance matrices. D=1: sigma2.
    
    a_0 <- prior_list$a_0
    b_0 <- prior_list$b_0
    
    for (k in seq_len(K)) {
      y_k <- data[c_alloc == k]
      n_k <- length(y_k)
      
      # Posterior for weights
      alpha_post[k] <- alpha_0 + n_k
      
      ybar <- mean(y_k)
      S <- sum((y_k - ybar)^2)    # within-cluster sum of squares
      
      # Posterior hyperparameters
      kappa_post <- kappa_0 + n_k
      mu_n <- (kappa_0 * mu_0 + n_k * ybar) / kappa_post
      a_post <- a_0 + 0.5 * n_k
      b_post <- b_0 + 0.5 * S +
        (kappa_0 * n_k) / (2 * kappa_post) * (ybar - mu_0)^2
      
      if (use_post_mean) {
        ## Plug-in posterior means
        var_post[[k]] <- b_post / (a_post - 1)   # mean of Inv-Gamma
        mu_post[[k]]  <- mu_n
      } else {
        ## Sample from posteriors
        var_post[[k]] <- 1 / rgamma(1, shape = a_post, rate = b_post)
        mu_post[[k]]  <- rnorm(1, mean = mu_n, sd = sqrt(var_post[[k]] / kappa_post))
      }
    }
  
    # ---------------------------------------------------------------------
    # CASE 2: MULTIVARIATE (D > 1): Normal–Inverse-Wishart
    # ---------------------------------------------------------------------
  } else {
    
    # Allocate storage
    alpha_post <- numeric(K)
    mu_post <- vector("list", K)
    var_post <- vector("list", K)   # D>1: covariance matrices. D=1: sigma2.
    
    nu_0 <- prior_list$nu_0
    S_0  <- prior_list$S_0
    
    for (k in seq_len(K)) {
      y_k <- data[c_alloc == k, , drop = FALSE]
      n_k <- nrow(y_k)
      
      # Posterior for weights
      alpha_post[k] <- alpha_0 + n_k
      
      # Means
      y_bar <- colMeans(y_k)
      diff <- sweep(y_k, 2, y_bar)
      S_k_data <- t(diff) %*% diff
      
      # Posterior parameters
      kappa_post <- kappa_0 + n_k
      mu_post_k <- (kappa_0 * mu_0 + n_k * y_bar) / kappa_post
      
      nu_post <- nu_0 + n_k
      cross_term <- (kappa_0 * n_k / kappa_post) * tcrossprod(y_bar - mu_0)
      S_k_post <- S_0 + S_k_data + cross_term
      
      if (use_post_mean) {
        ## Plug-in posterior means
        mu_post[[k]]  <- mu_post_k
        var_post[[k]] <- S_k_post / (nu_post - ncol(data) - 1)  # mean of Inv-Wishart
      } else {
        ## Sample from posteriors
        var_post[[k]] <- MCMCpack::riwish(nu_post, S_k_post)
        mu_post[[k]]  <- MASS::mvrnorm(1, mu_post_k, var_post[[k]] / kappa_post)
      }
    }
  }
  
  if (use_post_mean) {
    weight_post <- alpha_post / sum(alpha_post)
  } else {
    weight_post <- as.numeric(MCMCpack::rdirichlet(1, alpha_post))
  }
  
  list(
    weight_post   = weight_post,
    mu_post   = mu_post,
    var_post = var_post
  )
  
}

# Decides between MLE and posterior

get_mixture_params <- function(
    clustering,
    data,
    estimator = c("posterior", "mle"),
    prior_list = NULL,
    posterior_params = NULL,
    candidate_id = NULL,
    
    # --- variance prior hyperparameters ---
    alpha_0 = 2,   # 1D inverse-gamma shape
    nu_0    = NULL){
  
  estimator <- match.arg(estimator)
  D <- if (is.null(dim(data))) 1 else ncol(data)
  
  # =====================================================
  # POSTERIOR ESTIMATOR
  # =====================================================
  
  if (estimator == "posterior") {
    
    if (is.null(prior_list)) {
      stop("For estimator = 'posterior', prior_list must be supplied.")
    }
    
    post <- mixture_posterior(clustering, data, prior_list)
    
    return(list(
      weight = post$weight_post,
      mean   = post$mu_post,
      var    = post$var_post
    ))
    
  } else if (estimator == "mle") {
    
    # =====================================================
    # MLE with conjugate variance prior
    # =====================================================
    
    N <- if (D == 1) length(data) else nrow(data)
    K <- max(clustering)
    
    weights <- numeric(K)
    means   <- if (D == 1) numeric(K) else vector("list", K) #matrix(NA, nrow = K, ncol = D)
    vars    <- if (D == 1) numeric(K) else vector("list", K)
    
    # -------------------------------------------------
    # Data-driven prior scale parameters
    # -------------------------------------------------
    
    if (D == 1) {
      
      beta_0 <- var(data)
      
    } else {
      
      if (is.null(nu_0)) {
        nu_0 <- D + 2
      }
      
      Psi_0 <- cov(data)
    }
    
    # -------------------------------------------------
    
    for (k in seq_len(K)) {
      
      idx <- clustering == k
      n_k <- sum(idx)
      
      weights[k] <- n_k / N
      y_k <- if (D == 1) data[idx] else data[idx, , drop = FALSE]
      
      if (D == 1) {
        
        means[k] <- mean(y_k)
        
        # ---------------------------------------------
        # OLD VERSION (kept for reference)
        # vars[k] <- if (n_k > 1) var(y_k) else 1e-6
        # ---------------------------------------------
        
        if (n_k > 1) {
          S_k <- sum((y_k - means[k])^2)
        } else {
          S_k <- 0
        }
        
        alpha_post <- alpha_0 + n_k / 2
        beta_post  <- beta_0  + S_k / 2
        
        vars[k] <- beta_post / (alpha_post - 1)
        
      } else {
        
        means[[k]] <- colMeans(y_k)
        
        # ---------------------------------------------
        # OLD VERSION (kept for reference)
        # vars[[k]] <- if (n_k > 1) cov(y_k) else diag(D) * 1e-6
        # ---------------------------------------------
        
        if (n_k > 1) {
          S_k <- cov(y_k) * (n_k - 1)
        } else {
          S_k <- matrix(0, D, D)
        }
        
        nu_post  <- nu_0 + n_k
        Psi_post <- Psi_0 + S_k
        
        vars[[k]] <- Psi_post / (nu_post - D - 1)
      }
    }
    
    return(list(
      weight = weights,
      mean   = means,
      var    = vars
    ))
  }
}



