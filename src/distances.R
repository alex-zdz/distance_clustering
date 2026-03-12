# distances.R
# Discrepancy measures between empirical distribution and Gaussian mixture model

# ──────────────────────────────────────────────────────────────────────────────
# Public interface
# ──────────────────────────────────────────────────────────────────────────────
#' Compute discrepancy between empirical distribution and Gaussian mixture
#'
#' @param data n × D numeric vector (D=1) or matrix (D>1)
#' @param weight length-K vector of mixture weights
#' @param mean K-vector (D=1) or list of K length-D vectors (D>1)
#' @param var K-vector of variances (D=1) or list of K D×D covariance matrices (D>1)
#' @param method one of `"KS"`, `"Pearson"`, `"Wasserstein"`
#' @param clusterwise logical; if `TRUE` returns weighted sum of per-component distances
#' @param c_alloc integer vector of cluster labels (length n); required if `clusterwise = TRUE`
#' @param L number of random projections for sliced Wasserstein (only when D>1 and theta = NULL)
#' @param M grid size or number of quantile points (only for Wasserstein)
#' @param theta L × D matrix of unit vectors for sliced Wasserstein (optional; if NULL, generated internally)
#'
#' @return scalar distance value  
#'   (when `clusterwise = TRUE` the return value carries attributes `by_cluster`, `weight`, `method`, `clusterwise`)
compute_distance <- function(data,
                             weight,
                             mean,
                             var,
                             method      = c("KS", "Pearson", "Wasserstein"),
                             clusterwise = FALSE,
                             c_alloc     = NULL,
                             L           = 1000,
                             M           = 2000,
                             theta       = NULL) {
  
  method <- match.arg(method)
  
  if (is.null(dim(data))) {
    D <- 1L
    data <- as.vector(data)
  } else {
    D <- ncol(data)
  }
  
  K <- length(weight)
  
  if (D > 1L && method != "Wasserstein") {
    stop("Only sliced Wasserstein is implemented for multivariate data (D > 1)")
  }
  
  if (clusterwise) {
    if (is.null(c_alloc)) stop("clusterwise = TRUE requires c_alloc")
    return(.clusterwise_distance(data, c_alloc, weight, mean, var,
                                 method = method, M = M, L = L, D = D, theta = theta))
  }
  
  # Global dispatch
  switch(method,
         "KS"      = .ks_global(data, weight, mean, var),
         "Pearson" = .pearson_global(data, weight, mean, var),
         "Wasserstein" = .wasserstein_global(data, weight, mean, var,
                                             D = D, L = L, M = M, theta = theta))
}


# ──────────────────────────────────────────────────────────────────────────────
# Method-specific global dispatch helpers
# ──────────────────────────────────────────────────────────────────────────────

.ks_global <- function(data, weight, mean, var) {
  KS_mixture_distance(data, weight, mean, var)
}

.pearson_global <- function(data, weight, mean, var) {
  Pearson_mixture_distance(data, weight, mean, var)
}

.wasserstein_global <- function(data, weight, mean, var, D, L, M, theta) {
  if (D == 1L) {
    .wasserstein_1d(data, weight, mean, var, M = M)
  } else {
    .sliced_wasserstein_global(data, weight, mean, var, L = L, M = M, theta = theta)
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Low-level Wasserstein implementations (univariate) — global always uses grid
# ──────────────────────────────────────────────────────────────────────────────

.wasserstein_1d <- function(y, weight, mu, sig2, M) {
  .wasserstein_1d_grid(y, weight, mu, sig2, M)
}

.wasserstein_1d_grid <- function(y, weight, mu, sig2, M) {
  y_range <- range(y)
  t <- seq(y_range[1] - 5 * sqrt(max(sig2)),
           y_range[2] + 5 * sqrt(max(sig2)),
           length.out = M)
  
  pdf_model <- mixture_pdf_1D(t, mean = mu, sig2 = sig2, weight = weight)
  pdf_emp   <- emp_pdf(t, y)
  
  .wasserstein_from_densities(pdf_model, pdf_emp)
}

.w2_quantile_1d_single <- function(y, mu, sig2, M) {
  u <- seq(0.0005, 0.9995, length.out = M)
  q_emp   <- quantile(y, probs = u, type = 8, names = FALSE)
  q_model <- qnorm(u, mean = mu, sd = sqrt(sig2))
  mean((q_emp - q_model)^2)
}


# ──────────────────────────────────────────────────────────────────────────────
# Low-level Wasserstein implementations (multivariate sliced)
# ──────────────────────────────────────────────────────────────────────────────

.sliced_wasserstein_global <- function(data, weight, mean, var, L, M, theta) {
  # Global mixture distance: always grid-based
  D <- ncol(data)
  
  if (is.null(theta)) {
    theta <- generate_theta_normal(L, D)
  } else {
    L <- nrow(theta)
  }
  
  yproj <- data %*% t(theta)
  
  total <- 0
  for (l in seq_len(L)) {
    proj_mu <- vapply(seq_along(weight),
                      function(k) as.numeric(theta[l, ] %*% mean[[k]]), 0)
    proj_sd <- vapply(seq_along(weight), function(k) {
      sqrt(theta[l, , drop = FALSE] %*% var[[k]] %*% t(theta[l, , drop = FALSE]))
    }, 0)
    
    y_range <- range(yproj[, l])
    t <- seq(y_range[1] - 3 * max(proj_sd),
             y_range[2] + 3 * max(proj_sd),
             length.out = M)
    
    pdf_model <- mixture_pdf_1D(t, proj_mu, sig2 = proj_sd^2, weight = weight)
    pdf_emp   <- emp_pdf(t, yproj[, l])
    w_l <- .wasserstein_from_densities(pdf_model, pdf_emp)
    
    total <- total + w_l / L
  }
  total
}

.sliced_wasserstein_single <- function(y, mu, Sigma, L, M, theta) {
  # Clusterwise single Gaussian: always quantile-based
  D_dim <- ncol(y)
  
  if (is.null(theta)) {
    theta <- generate_theta_normal(L, D_dim)
  } else {
    L <- nrow(theta)
  }
  
  yproj <- y %*% t(theta)
  
  total <- 0
  for (l in seq_len(L)) {
    proj_mu <- as.numeric(theta[l, ] %*% mu)
    proj_sd <- sqrt(theta[l, , drop = FALSE] %*% Sigma %*% t(theta[l, , drop = FALSE]))
    
    w_l <- .w2_quantile_1d_single(yproj[, l], proj_mu, proj_sd^2, M)
    total <- total + w_l / L
  }
  total
}

.wasserstein_from_densities <- function(pdf_model, pdf_emp) {
  eps <- 1e-8
  pdf_model <- pdf_model + eps; pdf_model <- pdf_model / sum(pdf_model)
  pdf_emp   <- pdf_emp   + eps; pdf_emp   <- pdf_emp   / sum(pdf_emp)
  
  cdf_model <- cumsum(pdf_model)
  cdf_emp   <- cumsum(pdf_emp)
  
  x <- seq_along(pdf_model) - 1
  u <- seq(0, 1, length.out = length(pdf_model))
  
  q_model <- approx(cdf_model, x, xout = u, rule = 2)$y
  q_emp   <- approx(cdf_emp,   x, xout = u, rule = 2)$y
  
  du <- approx(q_model, q_model - q_emp, xout = x, rule = 2)$y
  (mean(abs(du)^2 * pdf_model))^(1/2)
}


# ──────────────────────────────────────────────────────────────────────────────
# Cluster-wise wrapper
# ──────────────────────────────────────────────────────────────────────────────

.clusterwise_distance <- function(data, c_alloc, weight, mean, var,
                                  method, M, L, D, theta) {
  
  if (D > 1L && method != "Wasserstein") {
    stop("clusterwise for D > 1 is only supported with method = 'Wasserstein' (sliced)")
  }
  
  K <- length(weight)
  Dk <- numeric(K)
  names(Dk) <- paste0("cluster_", seq_len(K))
  
  n_per_k <- tabulate(c_alloc, nbins = K)
  
  for (k in seq_len(K)) {
    
    if (D == 1L) {
      yk <- data[c_alloc == k]
    } else {
      yk <- data[c_alloc == k, , drop = FALSE]
    }
    
    nk <- n_per_k[k]
    
    if (nk == 1L) {
      if (D == 1L) {
        sigma_k <- sqrt(var[[k]])
        Dk[k] <- sigma_k
      } else {
        diff_sq <- sum((yk[1, ] - mean[[k]])^2)
        tr_Sigma <- sum(diag(var[[k]]))
        Dk[k] <- sqrt(diff_sq + tr_Sigma)
      }
      next
    }
    
    if (method == "Wasserstein") {
      if (D == 1L) {
        Dk[k] <- .wasserstein_1d(yk, weight = 1, mu = mean[k],
                                 sig2 = var[[k]], M = M)
      } else {
        Dk[k] <- .sliced_wasserstein_single(yk, mean[[k]], var[[k]],
                                            L = L, M = M, theta = theta)
      }
    } else if (method == "KS") {
      Dk[k] <- KS_mixture_distance(yk, weight = 1, mean = mean[k], sig2 = var[[k]])
    } else if (method == "Pearson") {
      Dk[k] <- Pearson_mixture_distance(yk, weight = 1, mean = mean[k], sig2 = var[[k]])
    }
  }
  
  total <- sum(weight * Dk)
  
  structure(total,
            total       = total,
            by_cluster  = Dk,
            weight      = weight,
            method      = method,
            clusterwise = TRUE)
}


# ──────────────────────────────────────────────────────────────────────────────
# Auxiliary functions
# ──────────────────────────────────────────────────────────────────────────────

emp_pdf <- function(t, proj) {
  h <- hist(proj,
            breaks = seq(min(t), max(t), length.out = length(t) + 1),
            plot = FALSE)
  density <- h$counts
  density / sum(density)
}

KS_mixture_distance <- function(data, weight, mean, sig2) {
  as.numeric(ks.test(data, mixture_cdf_1D, weight, mean, sig2)$statistic)
}

Pearson_mixture_distance <- function(data, weight, mean, sig2) {
  n <- length(data)
  Fn <- ecdf(data)
  bins <- seq(-50, 50, by = 0.5)
  n_bins <- length(bins) + 1
  
  bin_counts <- tabulate(cut(knots(Fn), breaks = c(-Inf, bins, Inf)), nbins = n_bins)
  gg <- bin_counts / n
  
  cdf_vals <- mixture_cdf_1D(bins, weight, mean, sig2)
  ff <- numeric(n_bins)
  ff[1] <- cdf_vals[1]
  ff[2:(n_bins-1)] <- diff(cdf_vals)
  ff[n_bins] <- 1 - cdf_vals[n_bins-1]
  
  ff[ff == 0] <- 1e-50
  
  sum((gg - ff)^2 / ff)
}


# ──────────────────────────────────────────────────────────────────────────────
# Projection generators
# ──────────────────────────────────────────────────────────────────────────────

generate_theta_normal <- function(L, d) {
  theta <- matrix(rnorm(L * d), nrow = L, ncol = d)
  theta <- theta / sqrt(rowSums(theta^2))
  theta
}