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
  if (is.null(means)) means <- seq(-2, 2, length.out = K)
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

#' Kolmogorov–Smirnov Distance for a Mixture Model
#'
#' Computes the KS distance between observed data and the CDF of a Gaussian mixture model.
#'
#' @param data Numeric vector. Observed data.
#' @param w Numeric vector of component weights.
#' @param mean Numeric vector of component means.
#' @param sig2 Numeric vector of component variances.
#'
#' @return Numeric scalar. KS distance.
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100)
#' KS_distance_mixture(data, w = c(0.5, 0.5), mean = c(0, 2), sig2 = c(1, 1))
#'
#' @export
KS_distance_mixture <- function(data, w, mean, sig2) {
  mixture_cdf <- function(x) sapply(x, function(xi) sum(w * pnorm(xi, mean = mean, sd = sqrt(sig2))))
  as.numeric(ks.test(data, mixture_cdf)$statistic)
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
compute_distance <- function(data, w, mean, sig2, method = "KS") {
  if (method == "KS") {
    KS_distance_mixture(data, w, mean, sig2)
  } else {
    stop("Unknown distance method")
  }
}


#' Initialize Clustering Based on Posterior Draws
#'
#' Selects the best posterior draw (clustering initialization) by minimizing
#' the chosen distance measure between observed data and the fitted mixture model.
#'
#' @param data Numeric vector. Observed data.
#' @param w_post List of numeric vectors. Posterior draws of mixture weights.
#' @param m_post List of numeric vectors. Posterior draws of component means.
#' @param sig2_post List of numeric vectors. Posterior draws of component variances.
#' @param method Character. Distance method to use. Currently only `"KS"` is implemented.
#'
#' @return Integer index of the posterior draw with the smallest distance.
#'
#' @examples
#' \dontrun{
#' # Suppose you have posterior draws:
#' data <- rnorm(100)
#' w_post <- list(c(0.4, 0.6), c(0.5, 0.5))
#' m_post <- list(c(0, 2), c(0, 1.8))
#' sig2_post <- list(c(1, 1), c(1, 1.2))
#' initialize_clustering(data, w_post, m_post, sig2_post, method = "KS")
#' }
#'
#' @export
initialize_clustering <- function(data, w_post, m_post, sig2_post, method = "KS") {
  n_mcmc <- length(w_post)          # number of posterior draws
  K <- length(w_post[[1]])          # number of mixture components
  distance <- numeric(n_mcmc)
  
  for (i in seq_len(n_mcmc)) {
    w_i    <- w_post[[i]]
    u_i    <- unlist(m_post[[i]])
    sig2_i <- unlist(sig2_post[[i]])
    distance[i] <- compute_distance(data, w_i, u_i, sig2_i, method)
  }
  
  which.min(distance)
}



Sweetening (add roxygen)
sweetening <- function(c_current, data, prior_list, D_current,
                       max_sweet_iter = 100, tol = 1e-10, method = "KS") {
  
  n_sweet <- 0
  
  while (n_sweet < max_sweet_iter) {
    n_sweet <- n_sweet + 1
    
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
                                      w = param_post$weights,
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


