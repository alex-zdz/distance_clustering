# functions for W2-Salso
library(LaplacesDemon)
library(Rcpp)

# source("W2-functions.R")
sourceCpp("W2-functions.cpp")

fulfill_gap_label <- function(c_vec){
  # relabel c_vec to make sure there are no empty clusters
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

# Initialization
## run W2 through all MCMC samples, select the one that minimizing W2
Initialization_phase <- function(c_sample, observed_data, prior_list){
  # Input: mat I_mcmc*n c_sample
  #        double [0,1] p_sa: the probability to sequential allocation
  #        int k_d: the number of maximum clusters
  print("Start Initialization Phase")
  result = min_W2(c_sample, observed_data, prior_list) + 1
  print("End Initialization Phase")
  return(result)
}

# Sweetening
add_back_c <- function(i, c_i, c_estimate_minus){
  # function to add label c_i to the i-th position in c_estimate_minus
  # Input: int i: position, 
  #        int c_i: label, 
  #        vec c_estimate_minus: c_vec with c_i removed
  # Output: vec c_new: c_vec with c_i add back
  
  n_minus = length(c_estimate_minus)
  c_new = rep(NA, n_minus+1)

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
  
  return(c_new)
}


argmin_W2 <- function(i, c_estimate_minus, observed_data, prior_list){
  # Input: int i: the position of the current allocoting c_i, 
  #        vec c_estimate_minus: c_vec with c_i removed 
  #        
  #        func loss function
  # Output: the cluster label that minimizes the posterior expected loss
  k_max = length(unique(c_estimate_minus))
  # pi_mat = posterior_similarity_matrix(c_sample_minus_add_i)
  f = rep(0, k_max+1)
  if(k_max==0){
    newlabel = 1
  }else{
    for (k in 1:(k_max+1)) {
      c_candidate = add_back_c(i, c_i=k, c_estimate_minus)
      f[k] = W2(data=observed_data, c_candidate, prior_list)
        # W2_function(add_back_c(i, c_i=k, c_estimate_minus,SA), c_sample_minus_add_i,pi_mat)
    }
    # c_candidate = add_back_c(i, c_i=k_max+1, c_estimate_minus,SA)
    # f[k_max+1] =  W2(data=observed_data, c_candidate, prior_list)
      # loss_function(add_back_c(i, c_i=k_max+1, c_estimate_minus,SA), c_sample_minus_add_i,pi_mat)
    newlabel = min(which(f == min(f)))
  }
  return(newlabel)
}

# Sweetening_phase <- function(c_old, observed_data, prior_list, tol=1e-3){
#   # update use c_last
#   print("Start Sweetening Phase")
#   n = length(c_old)
#   c_last = c_old
#   c_new = rep(NA, n)
#   n_iter = 1
#   # while (TRUE) {
#   # while (n_iter<5) {
#   while(n_iter<150){
#     
#     print("Sweetening_phase: iter")
#     print(n_iter)
#     W2_old =  W2(observed_data, c_last, prior_list)
#     permutation_subject_number = sample(1:n,n, replace=FALSE)
#     
#     for (j in 1:n) {
#       i = permutation_subject_number[j]  # i: subject number
#       # remove c_i
#       c_estimate_minus = c_last[-i]
#       # compute loss
#       newlabel = argmin_W2(i, c_estimate_minus, observed_data, prior_list)
#       # allocate
#       c_new[i] = newlabel
#     }
#     W2_new =  W2(observed_data, c_new, prior_list)
#     W2_diff = abs(W2_new - W2_old)
#     
#     c_last=c_new
#     # if(sum(c_last==c_new)==n){
#     if(W2_diff<tol){
#       break
#     }else{
#       n_iter = n_iter + 1
#     }
#   }
#   c_new = fulfill_gap_label(c_new)
#   print("End Sweetening Phase")
#   return(c_new)
# }


Sweetening_phase <- function(c_old, observed_data, prior_list, tol=1e-15){
  # update use c_new
  print("Start Sweetening Phase")
  n = length(c_old)
  c_last = c_old
  c_new = c_last
  n_iter = 1
  # while (TRUE) {
    # while (n_iter<5) {
  while(n_iter<150){

    print("Sweetening_phase: iter")
    print(n_iter)
    W2_old =  W2(observed_data, c_new, prior_list)
    print("W2_old")
    print(W2_old)
    permutation_subject_number = sample(1:n,n, replace=FALSE)

    for (j in 1:n) {
      i = permutation_subject_number[j]  # i: subject number
      # remove c_i
      c_estimate_minus = c_new[-i]
      # compute loss
      newlabel = argmin_W2(i, c_estimate_minus, observed_data, prior_list)
      # allocate
      c_new[i] = newlabel
    }
    W2_new =  W2(observed_data, c_new, prior_list)
    print("W2_new")
    print(W2_new)
    W2_diff = abs(W2_new - W2_old)

    # if(sum(c_last==c_new)==n){
    if(W2_diff<tol){
      break
    }else{
      # c_last=c_new
      n_iter = n_iter + 1
    }
  }
  c_new = fulfill_gap_label(c_new)
  print("End Sweetening Phase")
  return(c_new)
}


# Zealous updates phase

merge_cluster <- function(c_vec, clus1, clus2){
  # function to merge clus1 and clus2 in c_Vec
  sub_c1 = which(c_vec==clus1)
  c_vec[sub_c1] = clus2
  c_vec = fulfill_gap_label(c_vec)

  return(c_vec)
}



split_cluster <- function(c_vec,clus){
  # function to split a cluster into two clusters randomly in c_vec
  sub_c1 = which(c_vec==clus)
  n_unique_label = length(unique(c_vec))
  for (i in 1:length(sub_c1)) {
    U = runif(1,0,1)
    if(U<0.5){
      c_vec[sub_c1[i]] = n_unique_label+1
    }
  } 
  c_vec = fulfill_gap_label(c_vec)

  return(c_vec)
}


Zealous_updates_phase <- function(c_old, observed_data, prior_list, n_max=3, k=3){
  # merge and split
  
  print("Start Zealous Updates Phase")
  
  c_current = c_old
  n = length(c_current)
  for (i in 1:n_max) {
    print("n_zea=")
    print(i)
    
    n_unique_label = length(unique(c_current))
    current_loss = W2(data=observed_data, c_current, prior_list)
    
    n_merge_accept = 0
    if(n_unique_label>1){
      pairs = t(combn(1:n_unique_label, 2))
      for (j in 1:dim(pairs)[1]) {
        clusters_to_merge = pairs[j,]
        c_merge = merge_cluster(c_current, clusters_to_merge[1], clusters_to_merge[2])
        new_loss = W2(observed_data, c_merge, prior_list)
        if(new_loss<current_loss){
          c_current = c_merge
          # update current loss
          current_loss = W2(data=observed_data, c_current, prior_list)
          print("merge accept!")
          n_merge_accept = n_merge_accept+1
        }
      }
    }
    n_unique_label = length(unique(c_current))
    
    n_split_accept = 0
    for (i_k in 1:k) {
      # randomly split k times
      clus = sample(1:n_unique_label, 1)
      c_split = split_cluster(c_current,clus)
      new_loss = W2(observed_data, c_split, prior_list)
      if(new_loss<current_loss){
        c_current = c_split
        # update current loss
        current_loss = W2(data=observed_data, c_current, prior_list)
        print("split accept!")
        n_split_accept = n_split_accept+1
      }
    }
    if(n_merge_accept + n_split_accept == 0){
      break
    }
  }  
  print("End Zealous Updates Phase")
  return(c_current)
}


# Zealous_updates_phase <- function(c_old, observed_data, prior_list, n_maxzea=3, k){
#   # merge and split
#   
#   print("Start Zealous Updates Phase")
#   
#   c_current = c_old
#   
#   n = length(c_current)
#   
#   
#   for (i in 1:n_maxzea) {
#     n_unique_label = length(unique(c_current))
#     
#     if(n_unique_label>1){ # if there are more than 1 cluster, merge two and then split once
#       
#       # merge
#       current_loss = W2(data=observed_data, c_current, prior_list)
#       
#       for (i_k in 1:k) {
#         c_new_i = merge_cluster(c_current)
#         new_loss = W2(observed_data, c_new_i, prior_list)
#         if(new_loss<current_loss){
#           c_current = c_new_i
#           # update current loss
#           current_loss = W2(data=observed_data, c_current, prior_list)
#           print("merge!")
#         }
#       }
#       
#       # split
#       for (i_k in 1:k) {
#         c_new_i = split_cluster(c_current)
#         new_loss = W2(observed_data, c_new_i, prior_list)
#         if(new_loss<current_loss){
#           c_current = c_new_i
#           # update current loss
#           current_loss = W2(data=observed_data, c_current, prior_list)
#           print("split!")
#         }
#       }
#       
#       
#     }else{
#       # only split
#       for (i_k in 1:k) {
#         c_new_i = split_cluster(c_current)
#         new_loss = W2(observed_data, c_new_i, prior_list)
#         if(new_loss<current_loss){
#           c_current = c_new_i
#           # update current loss
#           current_loss = W2(data=observed_data, c_current, prior_list)
#           print("split!")
#         }
#       }
#       
#     }
#   }
#   print("End Zealous Updates Phase")
#   return(c_current)
# }




# reallocate_one_cluster <- function(reallo_label, c_vec, k){
#   # destory one cluster and reallocate the subjects to other clusters k times
#   
#   sub_j = which(c_vec==reallo_label)
#   n_subj = length(sub_j)
#   c_remove_j = c_vec[-sub_j]
#   cluster_left = unique(c_remove_j)
#   n_left = length(cluster_left)
#   # cluster_new = c(cluster_left, max(cluster_left)+1)
#   # n_new = length(cluster_new)
#   
#   new_c_vec_list = vector("list", k)
#   
#   for (i_k in 1:k) {
#     c_vec_k = c_vec
#     num_alloc = rcat(n_subj,p=rep(1/n_left, n_left))
#     c_vec_k[sub_j] = cluster_left[num_alloc]
#     c_vec_k = fulfill_gap_label(c_vec_k)
#     new_c_vec_list[[i_k]] = c_vec_k
#   }
#   
#   return(new_c_vec_list)
# }

# Zealous_updates_phase <- function(c_old, observed_data, prior_list, n_maxzea=3, k){
#   print("Start Zealous Updates Phase")
# 
#   c_current = c_old
#   for (i in 1:n_maxzea) {
#     if(length(unique(c_current))>1){
#     
#       n = length(c_current)
#       n_unique_label = length(unique(c_current))
#       reallocate_cluster = sample(1:n_unique_label, 1, replace=FALSE)
#     
#       current_loss = W2(data=observed_data, c_current, prior_list)
#       new_allocations_list = reallocate_one_cluster(reallocate_cluster, c_current, k)
#       
#       new_loss = rep(NA,k+1)
#       for (i_k in 1:k) {
#         c_new_i = new_allocations_list[[i_k]]
#         new_loss[i_k] = W2(observed_data, c_new_i, prior_list)
#       }
#       new_loss[k+1] = current_loss
#       # print("new_loss")
#       # print(new_loss)
#       min_partition_num = min(which(new_loss == min(new_loss)))
#       if (min_partition_num<=k){
#         print("update!")
#         c_current = new_allocations_list[[min_partition_num]]
#       }
#     }
#   }
#   print("End Zealous Updates Phase")
#   return(c_current)
# }



Salso_whole_procedure <- function(c_sample, observed_data,prior_list, n_runs, n_maxzea=0 ){
  c_score = rep(0, n_runs)
  c_list = vector(mode="list", length=n_runs)
  
  for (i in 1:n_runs){
    c_first = Initialization_phase(c_sample, observed_data, prior_list)
    c_second = Sweetening_phase(c_first, observed_data, prior_list)
    if(n_maxzea>0){
      c_third = Zealous_updates_phase(res_swt, observed_data, prior_list, n_maxzea, k=10)
    }else{
      c_third = c_second
    }
    c_list[[i]] = c_third
    c_score[i] = W2(observed_data, c_third, prior_list)
  }
  # print(c_list)
  # print(c_score)
  select_no = min(which(c_score==min(c_score)))
  
  return(c_list[[select_no]])
}



# result = reallocate_one_cluster(1, c(1,2,2,5,5,4,3,1), 20)
# Zealous_updates_phase(c_old, observed_data, prior_list, n_maxzea=3, k)