Sequential_Allocation <- function(c_estimate_minus, c_sample_minus, c_sample_add){
  n = dim(c_sample_add)[2]
  permutation_subject_number = sample(1:n,n, replace=FALSE)
  result = rep(NA, n)
  allocated_label = c()
  for (j in 1:n) {
    i = permutation_subject_number[j]  # i: subject number
    c_sample_minus_add_i = as.matrix(cbind(c_sample_minus,c_sample_add[,permutation_subject_number[1:j]]))
    newlabel = argmin_loss(i,c_estimate_minus, c_sample_minus_add_i, loss_function=Binder_loss_posterior_expectation, SA=TRUE)
    c_estimate_minus = c(c_estimate_minus, newlabel)
    allocated_label = c(allocated_label, newlabel)
  }
  result[permutation_subject_number] = allocated_label
  return(result)
}



Initialization_phase <- function(c_sample, p_sa, k_d){
  # Input: mat I_mcmc*n c_sample
  #        double [0,1] p_sa: the probability to sequential allocation
  #        int k_d: the number of maximum clusters
  print("Start Initialization Phase")
  n = dim(c_sample)[2]
  result = rep(NA, n)
  u = runif(n=1, min = 0, max = 1)
  if(u < p_sa){
    c_estimate_minus = c()
    c_sample_minus = matrix(NA, nrow = dim(c_sample)[1], ncol = 0)
    result = Sequential_Allocation(c_estimate_minus, c_sample_minus,c_sample)
  }else{
    result = sample(1:k_d, n, replace=TRUE)
  }
  result = fulfill_gap_label(result)
  print("End Initialization Phase")
  return(result)
}


posterior_similarity_matrix <- function(c_sample){
  # c_sample: I_mcmc * n matrix
  n = dim(c_sample)[2]
  H = dim(c_sample)[1]
  pi_mat = matrix(0, nrow=n, ncol=n)
  for (i in 1:n){
    for(j in 1:n){
      pi_mat[i,j] = sum(c_sample[,i]==c_sample[,j])/H
    }
  }
  return(pi_mat)
}

mixture_model_cdf <- function(x, weight1, mean1, sd1, mean2, sd2) {
  f = weight1 * pnorm(x, mean1, sd1) + (1-weight1) * pnorm(x, mean2, sd2)
  return(f)
}

Wasserstein2_sq <- function(observed_data, mixture_model_cdf, param_list){
  weight1 = param_list$weight1
  mean1 = param_list$mean1
  sd1 = param_list$sd1
  mean2 = param_list$mean2
  sd2 = param_list$sd2
  
  Y_ordered = sort(observed_data)
  vec_FminusG = rep(0, n)
  for (i in 1:n) {
    F_i = mixture_model_cdf(Y_ordered[i], weight1, mean1, sd1, mean2, sd2)
    G_i = i/n
    vec_FminusG[i] = F_i - G_i
  }
  W2_sq = sum((vec_FminusG)^2)
  return(W2_sq)
}

add_back_c <- function(i, c_i, c_estimate_minus, SA=FALSE){
  # Input: int i: position, 
  #        int c_i: label, 
  #        vec c_estimate_minus: c_vec with c_i removed
  #        bool SA: TRUE when sequential allocation, otherwise FALSE
  # Output: vec c_new: c_vec with c_i add back
  
  n_minus = length(c_estimate_minus)
  c_new = rep(NA, n_minus+1)
  if(SA){
    c_new = c(c_estimate_minus,c_i)
  }else{
    if(i==1){
      c_new[i] = c_i
      c_new[(i+1):(n_minus+1)] = c_estimate_minus[i:n_minus]
    }else if(i==n_minus+1){
      c_new[1:(i-1)] = c_estimate_minus[1:(i-1)]
      c_new[i] = c_i
    }else{
      c_new[1:(i-1)] = c_estimate_minus[1:(i-1)]
      c_new[i] = c_i
      c_new[(i+1):(n_minus+1)] = c_estimate_minus[i:n_minus]
    }
  }
  
  return(c_new)
}

argmin_loss <- function(i, c_estimate_minus, c_sample_minus_add_i, loss_function, SA=FALSE){
  # Input: int i: the position of the current allocoting c_i, 
  #        vec c_estimate_minus: c_vec with c_i removed 
  #        mat c_sample_minus_add_i: c_sample with columns corresponding to c_minus plus i
  #        func loss function
  # Output: the cluster label that minimizes the posterior expected loss
  k_max = length(unique(c_estimate_minus))
  pi_mat = posterior_similarity_matrix(c_sample_minus_add_i)
  f = rep(0, k_max+1)
  if(k_max==0){
    newlabel = 1
  }else{
    for (k in 1:k_max) {
      f[k] = loss_function(add_back_c(i, c_i=k, c_estimate_minus,SA), c_sample_minus_add_i,pi_mat)
    }
    f[k_max+1] = loss_function(add_back_c(i, c_i=k_max+1, c_estimate_minus,SA), c_sample_minus_add_i,pi_mat)
    newlabel = min(which(f == min(f)))
  }
  return(newlabel)
}


Sweetening_phase <- function(c_old, c_sample){
  print("Start Sweetening Phase")
  n = length(c_old)
  c_last = c_old
  c_new = rep(NA, n)
  n_iter = 1
  while (TRUE) {
    print("Sweetening_phase: iter")
    print(n_iter)
    permutation_subject_number = sample(1:n,n, replace=FALSE)
    
    for (j in 1:n) {
      i = permutation_subject_number[j]  # i: subject number
      # remove c_i
      c_estimate_minus = c_last[-i]
      # compute loss
      newlabel = argmin_loss(i, c_estimate_minus, c_sample, 
                             loss_function=Binder_loss_posterior_expectation, SA=FALSE)
      # allocate
      # print("newlabel")
      # print(newlabel)
      c_new[i] = newlabel
    }
    
    if(sum(c_last==c_new)==n){
      break
    }else{
      c_last=c_new
      n_iter = n_iter + 1
    }
  }
  c_new = fulfill_gap_label(c_new)
  print("End Sweetening Phase")
  return(c_new)
}

fulfill_gap_label <- function(c_vec){
  n_unique_label = length(unique(c_vec))
  k_max = max(c_vec)
  if(n_unique_label==k_max){
    return(c_vec)
  }else{
    label_unique = unique(c_vec)
    c_vec_relabel = rep(NA,length(c_vec))
    for (i in 1:length(c_vec)) {
      c_vec_relabel[i] = which(label_unique==c_vec[i])
    }
    return(c_vec_relabel)
  }
}

Zealous_updates_phase <- function(c_old, c_sample, n_maxzea=10){
  print("Start Zealous Updates Phase")
  n = length(c_old)
  n_unique_label = length(unique(c_old))
  reallocate_cluster = sample(1:n_unique_label, min(n_maxzea, max(c_old)), replace=FALSE)
  pi_mat = posterior_similarity_matrix(c_sample)
  c_current = c_old
  for (j in 1:length(reallocate_cluster)) {
    current_loss = Binder_loss_posterior_expectation(c_current, c_sample, pi_mat, a=1, b=1)
    reallo_clu_j = reallocate_cluster[j]
    # removing cluster [reallo_clu_j]
    sub_j = which(c_current==reallo_clu_j)
    c_remove_j = c_current[-sub_j]
    c_sample_remove = as.matrix(c_sample[,-sub_j])
    c_sample_add =  as.matrix(c_sample[, sub_j])
    SA_result = Sequential_Allocation(c_remove_j, c_sample_remove, c_sample_add)
    c_candidate = c_current
    c_candidate[sub_j] = SA_result
    c_candidate = fulfill_gap_label(c_candidate)
    # print("c_sample")
    # print(dim(c_sample))
    candidate_loss = Binder_loss_posterior_expectation(c_candidate, c_sample, pi_mat, a=1, b=1)
    
    if(candidate_loss < current_loss){
      c_current = c_candidate
    }
  }
  print("End Zealous Updates Phase")
  return(c_current)
}

Salso_whole_procedure <- function(c_sample, n_runs){
  c_score = rep(0, n_runs)
  c_list = vector(mode="list", length=n_runs)
  pi_mat = posterior_similarity_matrix(c_sample)
  for (i in 1:n_runs){
    c_first = Initialization_phase(c_sample, p_sa=0.5,k_d=3)
    c_second = Sweetening_phase(c_first, c_sample)
    c_third = Zealous_updates_phase(c_second, c_sample, n_maxzea=10)
    c_list[[i]] = c_third
    c_score[i] = Binder_loss_posterior_expectation(c_third, c_sample, pi_mat, a=1, b=1)
  }
  # print(c_list)
  # print(c_score)
  select_no = min(which(c_score==min(c_score)))
  
  return(c_list[[select_no]])
}




# library(doParallel)
# registerDoParallel(cores=4)
# library(foreach)
# n_runs=4
# pi_mat = posterior_similarity_matrix(c_sample)

# Salso_ith_run <- function(c_sample,pi_mat){
#   c_first = Initialization_phase(c_sample, p_sa=0.5,k_d=3)
#   c_second = Sweetening_phase(c_first, c_sample)
#   c_third = Zealous_updates_phase(c_second, c_sample, n_maxzea=10)
#   c_score = Binder_loss_posterior_expectation(c_third, c_sample, pi_mat, a=1, b=1)
#   return(list(c_third,c_score))
# }
# 
# Salso_parallel <- function(c_sample,n_runs, Salso_ith_run){
#   pi_mat = posterior_similarity_matrix(c_sample)
#   x <- foreach(i=1:n_runs) %dopar% Salso_ith_run(c_sample,pi_mat)
#   c_score = rep(NA,n_runs)
#   for (i in 1:nruns) {
#     c_score[i] = x[[i]][[2]]
#   }
#   select_no = min(which(c_score==min(c_score)))
#   return(x[[select_no]][[1]])
# }



# Binder_loss <- function(c_true, c_estimated, a, b){
#   n = length(c_true)
#   loss = 0
#   for (i in 1:n) {
#     for (j in (i+1):n) {
#       I1 = (c_true[i]==c_true[j]) * 1.0
#       I2 = (c_estimated[i]!=c_estimated[j]) * 1.0
#       I3 = (c_true[i]!=c_true[j]) * 1.0
#       I4 = (c_estimated[i]==c_estimated[j]) * 1.0
#       
#       loss = loss + a * I1 * I2 + b * I3 * I4
#     }
#   }
#   return(loss)
# }
