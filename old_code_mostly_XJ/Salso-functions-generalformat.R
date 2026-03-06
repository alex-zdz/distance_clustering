# functions for salso with distances
# library(LaplacesDemon)
# library(Rcpp)

# source("Pearson-functions.R")


min_distance <- function(c_samples, data, prior_list, distance_function){
  # c_samples: M*n matrix
  M = dim(c_samples)[1]
  result_vec = numeric(M)
  for (i in 1:M) {
    if(i%%100==0){
      print(i)
    }
    result_vec[i] = distance_function(data, c_samples[i,], prior_list)
  }
  min_m = which(result_vec==min(result_vec))[1]
  # print(min_m)
  return(c_samples[min_m,])
}


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
Initialization_phase <- function(c_sample, observed_data, prior_list, distance_function){
  # Input: mat I_mcmc*n c_sample
  #        double [0,1] p_sa: the probability to sequential allocation
  #        int k_d: the number of maximum clusters
  # print("Start Initialization Phase")
  result = min_distance(c_sample, observed_data, prior_list,distance_function) + 1
  # print("End Initialization Phase")
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


argmin_distance <- function(i, c_estimate_minus, observed_data, prior_list,distance_function){
  # Input: int i: the position of the current allocoting c_i, 
  #        vec c_estimate_minus: c_vec with c_i removed 
  #        
  #        func loss function
  # Output: the cluster label that minimizes the posterior expected loss
  k_max = length(unique(c_estimate_minus))
  f = rep(0, k_max+1)
  if(k_max==0){
    newlabel = 1
  }else{
    for (k in 1:(k_max+1)) {
      c_candidate = add_back_c(i, c_i=k, c_estimate_minus)
      f[k] = distance_function(observed_data, c_candidate, prior_list)
    }

    newlabel = min(which(f == min(f)))
  }
  return(newlabel)
}



Sweetening_phase <- function(c_old, observed_data, prior_list,distance_function, tol=1e-15){
  # update use c_new
  # print("Start Sweetening Phase")
  n = length(c_old)
  c_last = c_old
  c_new = c_last
  n_iter = 1
  # while (TRUE) {
  # while (n_iter<5) {
  # print("Sweetening_phase: iter")
  while(n_iter<150){
  
  # to delete:
    #while(n_iter<10){
      #if(n_iter%%10 == 0){
        print(n_iter)
      #}
  
    D_old =  distance_function(observed_data, c_new, prior_list)
    # print("D_old")
    # print(D_old)
    permutation_subject_number = sample(1:n,n, replace=FALSE)
    
    for (j in 1:n) {
      i = permutation_subject_number[j]  # i: subject number
      # remove c_i
      c_estimate_minus = c_new[-i]
      # compute loss
      newlabel = argmin_distance(i, c_estimate_minus, observed_data, prior_list,distance_function)
      # allocate
      c_new[i] = newlabel
    }
    D_new =  distance_function(observed_data, c_new, prior_list)
    D_diff = abs(D_new - D_old)
    
    if(D_diff<tol){
      break
    }else{
      n_iter = n_iter + 1
    }
  }
  c_new = fulfill_gap_label(c_new)
  # print("End Sweetening Phase")
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


Zealous_updates_phase <- function(c_old, observed_data, prior_list,distance_function, n_max=3, k=3){
  # merge and split
  
  # print("Start Zealous Updates Phase")
  
  c_current = c_old
  n = length(c_current)
  for (i in 1:n_max) {
    # print("n_zea=")
    # print(i)
    
    n_unique_label = length(unique(c_current))
    current_loss = distance_function(observed_data, c_current, prior_list)
    
    n_merge_accept = 0
    if(n_unique_label>1){
      pairs = t(combn(1:n_unique_label, 2))
      for (j in 1:dim(pairs)[1]) {
        clusters_to_merge = pairs[j,]
        c_merge = merge_cluster(c_current, clusters_to_merge[1], clusters_to_merge[2])
        new_loss = distance_function(observed_data, c_merge, prior_list)
        if(new_loss<current_loss){
          c_current = c_merge
          # update current loss
          current_loss = distance_function(observed_data, c_current, prior_list)#redundant just use new_loss!
          # print("merge accept!")
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
      new_loss = distance_function(observed_data, c_split, prior_list)
      if(new_loss<current_loss){
        c_current = c_split
        # update current loss
        current_loss = distance_function(observed_data, c_current, prior_list)
        print("split accept!")
        n_split_accept = n_split_accept+1
      }
    }
    if(n_merge_accept + n_split_accept == 0){
      break
    }
  }  
  # print("End Zealous Updates Phase")
  return(c_current)
}




Salso_whole_procedure <- function(c_sample, observed_data, prior_list,distance_function, n_runs=3, n_max=3 ){
  c_score = rep(0, n_runs)
  c_list = vector(mode="list", length=n_runs)
  
  for (i in 1:n_runs){
    c_first = Initialization_phase(c_sample, observed_data, prior_list,distance_function)
    c_second = Sweetening_phase(c_first, observed_data, prior_list,distance_function)
    if(n_max>0){
      c_third = Zealous_updates_phase(c_second, observed_data, prior_list,distance_function, n_max, k=3)
    }else{
      c_third = c_second
    }
    c_list[[i]] = c_third
    c_score[i] = distance_function(observed_data, c_third, prior_list)
  }
  # print(c_list)
  # print(c_score)
  select_no = min(which(c_score==min(c_score)))
  
  return(c_list[[select_no]])
}



# result = reallocate_one_cluster(1, c(1,2,2,5,5,4,3,1), 20)
# Zealous_updates_phase(c_old, observed_data, prior_list, n_maxzea=3, k)