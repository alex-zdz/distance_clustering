# entire algorithm in one go:

mde_main <- function(c_sample , data, prior_list, distance_function, n_runs = 2, n_max = 2, n_split = 3 ){
  
  
  # remember to add 1 to c_sample if it comes from AM
  c_score = rep(0, n_runs)
  c_list = vector(mode = "list", length = n_runs)
  n_mcmc = dim(c_sample)[1]
  n = dim(c_sample)[2]
  
  for (run in 1:n_runs){
    
    
    # 1.Initialization_phase
    result_vec = numeric(n_mcmc)
    
    for (i in 1:n_mcmc) {
      if(i%%100==0){
        print(i)
      }
      
      result_vec[i] = distance_function(data, c_sample[i,], prior_list) # KS_distance_mixture

    }
    min_m = which.min(result_vec)
    D_current <- result_vec[min_m]
    c_first = c_sample[min_m, ]  # + 1 since the input from AntMAN starts at 0
    
    # 2.Sweetening_phase
    # Initialize the new clustering configuration
    c_new <- c_first
    n_iter <- 1
    
    # Iterate until convergence or maximum iterations reached
    while (n_iter < 100) {
      
      if(n_iter%%1==0){
        print(n_iter)
      }
      
      # Generate a random permutation of subject indices
      permutation_subject_number <- sample(1:n, n, replace = FALSE)
      
      # Loop over all subjects in randomized order
      for (j in 1:n) {
        i <- permutation_subject_number[j]  # Current subject index
        
        # Temporarily remove the current item from clustering
        c_estimate_minus <- c_new[-i]
        n_cluster <- max(c_estimate_minus)
        
        # Vector to store losses for alternative allocations
        alternative_losses <- numeric(n_cluster + 1)
        
        # Evaluate loss for placing the removed item into each cluster
        # ToDo: this can be parallelized: create and evaluate all possible candidates at once
        for (k in 1:(n_cluster + 1)) {
          c_candidate <- append(c_estimate_minus, k, after = i - 1)
          alternative_losses[k] <- distance_function(data, c_candidate, prior_list)
        }
        newlabel <- which.min(alternative_losses)
        
        # Update clustering with the optimal assignment
        c_new[i] <- newlabel
      }
      
      # Compute the new loss after allocation
      D_new <- min(alternative_losses)
      
      #if (isTRUE(all.equal(D_new, D_current))) {
      if (abs(D_new - D_current) < 1e-15) {
        break
      } else {
        n_iter <- n_iter + 1
      }
      
      D_current <- D_new
      
    }
    
    D_current = D_new
    c_current = fulfill_gap_label(c_new)
    n_cluster = max(c_current)
  
    #3. Zealous_updates_phase 
    # merge and split
    if(n_max > 0){
      for (i in 1:n_max) {
        # Merging step:
          n_merge_accept = 0
          if(n_cluster > 1){
            pairs = t(combn(1:n_cluster, 2))
            for (j in 1:dim(pairs)[1]) {
              clusters_to_merge = pairs[j,]
              c_merge = c_current
              c_merge[which(c_merge==clusters_to_merge[1])] = clusters_to_merge[2]
              c_merge = fulfill_gap_label(c_current)
              D_new = distance_function(data, c_merge, prior_list)
              if(D_new < D_current){
                c_current = c_merge
                n_cluster = max(c_current)
                D_current = D_new
                n_merge_accept = n_merge_accept+1
              }
            }
          }
          
          # Split step:
          n_split_accept = 0
          c_split = c_current
          for (i_k in 1:n_split) {
            # randomly split n_split times
            cl_to_split = sample(1:n_cluster, 1)
            where_splitclust <- which(c_current==cl_to_split)
            size_of_splitclust <- length(where_splitclust)
            for (i in 1:size_of_splitclust) {
              U = runif(1,0,1)
              if(U < 0.5){
                c_split[where_splitclust[i]] = n_cluster + 1
              }
            } 
            c_split = fulfill_gap_label(c_split)
            D_new = distance_function(data, c_split, prior_list)
            if(D_new < D_current){
              c_current = c_split
              n_cluster = max(c_current)
              D_current = D_new
              n_split_accept = n_split_accept+1
            }
          }
          
          if(n_merge_accept + n_split_accept == 0){
            break
          }
      } # end of zealous update
    }
    
    c_list[[run]] = c_current
    c_score[run] = D_current
    
  } # end of n_runs
  
  select_no = which.min(c_score)
  
  return(c_list[[select_no]])
}

##################################################
# Helper functions
#################################################

fulfill_gap_label <- function(c_vec){
  # relabel c_vec to make sure there are no empty clusters
  n_unique_label = length(unique(c_vec))
  n_cluster = max(c_vec)
  if(n_unique_label==n_cluster){
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

data_generation <- function(N, type="Gaussian mixture"){
  
  if(type=="Gaussian mixture"){
    # data generation
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    # hist(c_alloc_true)
    
    y = rep(NA, N)
    # m = c(20,10,0,5,-2)
    m = c(-2,0,3)
    # m = c(7,10,0,5,-2)
    s = c(0.4, 1, 0.3)
    for (i in 1:K) {
      I = c_alloc_true==i
      y[I] = rnorm(sum(I), m[i], s[i])
    }
    # hist(y)
  }else if(type=="Skew Gaussian mixture"){
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    y = rep(NA, N)
    m = c(-2,0,3)
    s = c(0.4, 1, 0.3)
    a = c(1,10,4)
    
    for (i in 1:K) {
      I = c_alloc_true==i
      y[I] = rsn(sum(I),  m[i], s[i], a[i])
    }
    # plot(density(y))
  }else if(type=="Skew-symmetric mixture"){
    weights_true = c(0.45,0.25,0.3)
    K=3
    c_alloc_true = sample(1:K, N, replace = TRUE, prob=weights_true)
    y = rep(NA, N)
    
    I2 = c_alloc_true==2
    y[I2] = rsn(sum(I2),  0, 5, 4)
    I3 = c_alloc_true==3
    y[I3] = rnorm(sum(I3),  -4, sqrt(0.5))
    
    I1 = c_alloc_true==1
    set1 = which(I1)
    w11 = c(0.364,0.212,0.424)
    c11 = sample(1:3, sum(I1), replace = TRUE, prob=w11)
    I11 = c11==1
    y[set1[I11]] = rsn(sum(I11),  2.5, 1, -10)
    I12 = c11==2
    y[set1[I12]] = rnorm(sum(I12),  2.325, sqrt(0.2))
    I13 = c11==3
    y[set1[I13]] = rnorm(sum(I13),  1.085, sqrt(0.7))
    
    # plot(density(y))
  }
  
  
  return(list(data=y, c_true = c_alloc_true))
}









