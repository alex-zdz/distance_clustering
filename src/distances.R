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
#' @param L number of random projections for sliced Wasserstein (only when D>1)
#' @param version character; `"grid"` or `"quantile"` — used **only** when `method = "Wasserstein"`
#' @param M grid size or number of quantile points (only for Wasserstein)
#'
#' @return scalar distance value  
#'   (when `clusterwise = TRUE` the return value carries attributes `by_cluster`, `weight`, `method`, `version`, `clusterwise`)
compute_distance <- function(data,
                             weight,
                             mean,
                             var,
                             method      = c("KS", "Pearson", "Wasserstein"),
                             clusterwise = FALSE,
                             c_alloc     = NULL,
                             L           = 100,
                             version     = c("grid", "quantile"),
                             M           = 2000) {
  
  method  <- match.arg(method)
  version <- match.arg(version)
  
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
                                 method = method, version = version,
                                 M = M, L = L, D = D))
  }
  
  # Global dispatch
  switch(method,
         "KS"      = .ks_global(data, weight, mean, var),
         "Pearson" = .pearson_global(data, weight, mean, var),
         "Wasserstein" = .wasserstein_global(data, weight, mean, var,
                                             D = D, L = L, version = version, M = M)
  )
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

.wasserstein_global <- function(data, weight, mean, var, D, L, version, M) {
  if (D == 1L) {
    .wasserstein_1d(data, weight, mean, var, version = version, M = M)
  } else {
    .sliced_wasserstein_global(data, weight, mean, var, L = L, version = version, M = M)
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Low-level Wasserstein implementations (univariate)
# ──────────────────────────────────────────────────────────────────────────────

.wasserstein_1d <- function(y, weight, mu, sig2, version, M) {
  if (version == "grid") {
    .wasserstein_1d_grid(y, weight, mu, sig2, M)
  } else {
    .wasserstein_1d_quantile(y, weight, mu, sig2, M)
  }
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

.wasserstein_1d_quantile <- function(y, weight, mu, sig2, M) {
  # Weighted sum of squared 2-Wasserstein distances to each component
  sum_w2 <- 0
  for (k in seq_along(weight)) {
    sum_w2 <- sum_w2 + weight[k] * .w2_quantile_1d_single(y, mu[k], sig2[k], M = M)
  }
  sqrt(sum_w2)
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

.sliced_wasserstein_global <- function(data, weight, mean, var, L, version, M) {
  D <- ncol(data)
  theta <- generate_theta_normal(L, D)
  yproj <- data %*% t(theta)
  
  total <- 0
  for (l in seq_len(L)) {
    proj_mu <- vapply(seq_along(weight), function(k) as.numeric(theta[l, ] %*% mean[[k]]), 0)
    proj_sd <- vapply(seq_along(weight), function(k) {
      sqrt(theta[l, , drop = FALSE] %*% var[[k]] %*% t(theta[l, , drop = FALSE]))
    }, 0)
    
    if (version == "grid") {
      y_range <- range(yproj[, l])
      t <- seq(y_range[1] - 3 * max(proj_sd),
               y_range[2] + 3 * max(proj_sd),
               length.out = M)
      pdf_model <- mixture_pdf_1D(t, proj_mu, sig2 = proj_sd^2, weight = weight)
      pdf_emp   <- emp_pdf(t, yproj[, l])
      w_l <- .wasserstein_from_densities(pdf_model, pdf_emp)
    } else {
      w_l <- 0
      for (k in seq_along(weight)) {
        w_l <- w_l + weight[k] * .w2_quantile_1d_single(yproj[, l], proj_mu[k], proj_sd[k]^2, M = M)
      }
      w_l <- sqrt(w_l)
    }
    total <- total + w_l / L
  }
  total
}

.sliced_wasserstein_single <- function(y, mu, Sigma, L, version, M) {
  # y : n_k × D matrix (single cluster)
  # mu, Sigma : parameters of ONE Gaussian
  D_dim <- ncol(y)
  theta <- generate_theta_normal(L, D_dim)
  yproj <- y %*% t(theta)
  
  total <- 0
  for (l in seq_len(L)) {
    proj_mu <- as.numeric(theta[l, ] %*% mu)
    proj_sd <- sqrt(theta[l, , drop = FALSE] %*% Sigma %*% t(theta[l, , drop = FALSE]))
    
    if (version == "grid") {
      y_range <- range(yproj[, l])
      t <- seq(y_range[1] - 3 * proj_sd,
               y_range[2] + 3 * proj_sd,
               length.out = M)
      pdf_model <- mixture_pdf_1D(t, proj_mu, sig2 = proj_sd^2, weight = 1)
      pdf_emp   <- emp_pdf(t, yproj[, l])
      w_l <- .wasserstein_from_densities(pdf_model, pdf_emp)
    } else {
      w_l <- .w2_quantile_1d_single(yproj[, l], proj_mu, proj_sd^2, M)
    }
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
# Cluster-wise wrapper (now supports D > 1 via sliced Wasserstein)
# ──────────────────────────────────────────────────────────────────────────────


.clusterwise_distance <- function(data, c_alloc, weight, mean, var,
                                  method, version, M, L, D) {
  
  # Guard: only Wasserstein is supported for multivariate clusterwise
  if (D > 1L && method != "Wasserstein") {
    stop("clusterwise for D > 1 is only supported with method = 'Wasserstein' (sliced)")
  }
  
  K <- length(weight)
  Dk <- numeric(K)
  names(Dk) <- paste0("cluster_", seq_len(K))
  
  # Number of points per cluster
  n_per_k <- tabulate(c_alloc, nbins = K)
  
  for (k in seq_len(K)) {
    
    # Extract cluster data
    if (D == 1L) {
      yk <- data[c_alloc == k]
    } else {
      yk <- data[c_alloc == k, , drop = FALSE]
    }
    
    nk <- n_per_k[k]
    
    # ── Singleton case (nk == 1): use exact W₂ distance for all methods ──────
    if (nk == 1L) {
      if (D == 1L) {
        # Univariate: σ̂
        sigma_k <- sqrt(var[[k]])
        Dk[k] <- sigma_k
      } else {
        # Multivariate: √(‖x - μ‖² + trace(Σ))
        diff_sq <- sum((yk[1, ] - mean[[k]])^2)
        tr_Sigma <- sum(diag(var[[k]]))
        Dk[k] <- sqrt(diff_sq + tr_Sigma)
      }
      next
    }
    
    # ── Normal case (nk ≥ 2) ─────────────────────────────────────────────────
    if (method == "Wasserstein") {
      if (D == 1L) {
        Dk[k] <- .wasserstein_1d(yk, weight = 1, mu = mean[k],
                                 sig2 = var[[k]], version = version, M = M)
      } else {
        Dk[k] <- .sliced_wasserstein_single(yk, mean[[k]], var[[k]],
                                            L = L, version = version, M = M)
      }
    } else if (method == "KS") {
      Dk[k] <- KS_mixture_distance(yk, weight = 1, mean = mean[k], sig2 = var[[k]])
    } else if (method == "Pearson") {
      Dk[k] <- Pearson_mixture_distance(yk, weight = 1, mean = mean[k], sig2 = var[[k]])
    }
  }
  
  # Final weighted sum
  total <- sum(weight * Dk)
  
  structure(total,
            total       = total,
            by_cluster  = Dk,
            weight      = weight,
            method      = method,
            version     = if (method == "Wasserstein") version else NULL,
            clusterwise = TRUE)
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



# compute_distance <- function(data, weight, mean, var,
#                              method = "KS",
#                              clusterwise = FALSE,
#                              c_alloc = NULL,
#                              L = 100) {
#   
#   if (is.null(dim(data))) {
#     D <- 1
#   } else {
#     D <- ncol(data)
#   }
#   
#   if (D != 1 && method != "Wasserstein") {
#     stop("only Sliced Wasserstein for higher dimensions")
#   }
#   
#   K <- length(weight)
#   
#   ## =========================
#   ## CLUSTER-WISE DISTANCE
#   ## =========================
#   if (clusterwise) {
#     
#     if (is.null(c_alloc)) {
#       stop("clusterwise = TRUE requires c_alloc")
#     }
#     
#     if (!is.null(dim(data))) {
#       stop("cluster-wise distance currently only implemented for 1D")
#     }
#     
#     Dk <- numeric(K)
#     names(Dk) <- paste0("cluster_", seq_len(K))
#     
#     for (k in seq_len(K)) {
#       
#       yk <- data[c_alloc == k]
#       
#       # # skip empty / degenerate clusters
#       # if (length(yk) < 2) {
#       #   Dk[k] <- 0
#       #   next
#       # }
#       
#       if (method == "KS") {
#         Dk[k] <- KS_mixture_distance(
#           yk,
#           weight = 1,
#           mean   = mean[k],
#           sig2   = var[k]
#         )
#         
#       } else if (method == "Pearson") {
#         Dk[k] <- Pearson_mixture_distance(
#           yk,
#           weight = 1,
#           mean   = mean[k],
#           sig2   = var[k]
#         )
#         
#       } else if (method == "Wasserstein") {
#         
#         M <- 2e3
#         t <- seq(
#           from = min(yk) - 3 * sqrt(var[k]),
#           to   = max(yk) + 3 * sqrt(var[k]),
#           length.out = M
#         )
#         
#         I0 <- mixture_pdf_1D(t, mean[k], sig2 = var[k], weight = 1)
#         I1 <- emp_pdf(t, yk)
#         
#         Dk[k] <- Wasserstein_1D(I0, I1, p = 2)
#         
#       } else {
#         stop("Unknown distance method")
#       }
#     }
#     
#     total <- sum(weight * Dk)
#     
#     attr(total, "total")  <- total
#     attr(total, "weight")  <- weight
#     attr(total, "by_cluster")  <- Dk
#     attr(total, "method")      <- method
#     attr(total, "clusterwise") <- TRUE
#     
#     return(total)
#   }
#   
#   ## =========================
#   ## GLOBAL DISTANCE
#   ## =========================
#   if (method == "KS") {
#     
#     KS_mixture_distance(data, weight, mean, sig2 = var)
#     
#   } else if (method == "Wasserstein") {
#     
#     M <- 2e3
#     t <- seq(
#       from = -max(abs(data)) * sqrt(2 * D),
#       to   =  max(abs(data)) * sqrt(2 * D),
#       length.out = M
#     )
#     
#     if (D == 1) {
#       
#       I0 <- mixture_pdf_1D(t, mean, sig2 = var, weight)
#       I1 <- emp_pdf(t, data)
#       Wasserstein_1D(I0, I1, p = 2)
#       
#     } else {
#       # sliced wasserstein
#       Sigma <- var
#       theta <- generate_theta_normal(L, D)
#       yproj <- data %*% t(theta)
#       
#       projected_Sigma <- matrix(0, nrow = K, ncol = L)
#       projected_Mu    <- matrix(0, nrow = K, ncol = L)
#       
#       for (k in seq_len(K)) {
#         for (l in seq_len(L)) {
#           projected_Sigma[k, l] <-
#             sqrt(theta[l, , drop = FALSE] %*% Sigma[[k]] %*% t(theta[l, , drop = FALSE]))
#           projected_Mu[k, l] <-
#             theta[l, , drop = FALSE] %*% mean[[k]]
#         }
#       }
#       
#       SW <- 0
#       for (l in seq_len(L)) {
#         RIx <- mixture_pdf_1D(t, projected_Mu[, l], projected_Sigma[, l], weight)
#         RIy <- emp_pdf(t, yproj[, l])
#         SW  <- SW + Wasserstein_1D(RIx, RIy, p = 2) / L
#       }
#       
#       SW
#     }
#     
#   } else if (method == "Pearson") {
#     
#     Pearson_mixture_distance(data, weight, mean, sig2 = var)
#     
#   } else {
#     stop("Unknown distance method")
#   }
# }






#####################
# Sliced Wasserstein
#####################
# 
# generateTheta <- function(L, d) {
#   theta <- matrix(0, nrow = L, ncol = d)
#   
#   # First vector: uniform on sphere
#   th_l <- runif(d)
#   th_l <- th_l / sqrt(sum(th_l^2))
#   theta[1, ] <- th_l
#   
#   # Remaining vectors
#   for (i in 2:L) {
#     th_l <- rnorm(d)
#     th_l <- th_l / sqrt(sum(th_l^2))
#     
#     m <- max(abs(theta[1:(i - 1), ] %*% th_l))
#     
#     while (m > 0.97) {
#       th_l <- rnorm(d)
#       th_l <- th_l / sqrt(sum(th_l^2))
#       m <- max(abs(theta[1:(i - 1), ] %*% th_l))
#     }
#     theta[i, ] <- th_l
#   }
#   theta
# }
# 
# 
# generate_theta_normal <- function(L, d) {
#   theta <- matrix(rnorm(L * d), nrow = L, ncol = d)
#   # Normalize each row to have length 1
#   theta <- theta / sqrt(rowSums(theta^2))
#   theta
# }
# 
# Wasserstein_1D <- function(I0, I1, p = 2) {
#   stopifnot(length(I0) == length(I1))
#   
#   eps <- 1e-7
#   
#   # Ensure strict positivity
#   I0 <- I0 + eps
#   I1 <- I1 + eps
#   
#   # Normalize to sum to 1 - both pdfs have already been normalized, so this step is redundant here
#   I0 <- I0 / sum(I0)
#   I1 <- I1 / sum(I1)
#   
#   # Compute CDFs
#   J0 <- cumsum(I0)
#   J1 <- cumsum(I1)
#   
#   # Grid
#   x       <- seq_along(I0) - 1
#   xtilde  <- seq(0, 1, length.out = length(I0))
#   
#   # Inverse CDF sampling (pseudo-quantiles)
#   XI0 <- approx(J0, x, xout = xtilde, rule = 2)$y
#   XI1 <- approx(J1, x, xout = xtilde, rule = 2)$y
#   
#   # Displacement field u(x)
#   u <- approx(XI0, XI0 - XI1, xout = x, rule = 2)$y
#   
#   # # Transport map f(x) = x - u(x)
#   # f <- x - u
#   # 
#   # # Potential phi(x)
#   # phi <- cumsum(u / length(I0))
#   # phi <- phi - mean(phi)
#   
#   # p-Wasserstein distance
#   Wp <- (mean((abs(u)^p) * I0))^(1/p)
#   
#   Wp
# }
