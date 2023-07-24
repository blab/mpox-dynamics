library(dplyr)
library(doParallel)
library(foreach)
library(igraph)

## Distribution of cluster of identical sequences
cluster_alloc <- readRDS('../data/cluster_alloc.rds')

size_largest_polytomy <- cluster_alloc$df_size_distrib$cluster_size %>% max()
n_clust_identical_sequences <- sum(cluster_alloc$df_size_distrib$count)
n_seq <- sum(cluster_alloc$df_size_distrib$count * cluster_alloc$df_size_distrib$cluster_size)

## Probability that transmission occurs before mutation for mpox
p_trans_before_mut <- 0.66259

#######
get_vec_proba_obs_cluster_size <- function(vec_cluster_size, R0, k, p, p_detect, max_cluster_size){
  
  get_log_proba_size_cluster <- function(cluster_size, R0, k, p){
    
    log_proba <- 
      lgamma(k * cluster_size + cluster_size - 1) - lgamma(k * cluster_size) - lgamma(cluster_size + 1) +
      (cluster_size - 1)*log(p * R0 / k) - (k * cluster_size + cluster_size - 1) * log(1 + p * R0 / k)
    
    return(log_proba)
  }
  
  get_vec_log_proba_cluster_size <- function(vec_cluster_size, R0, k, p){
    vec_log_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
      get_log_proba_size_cluster(cluster_size, R0, k, p)
    })
    return(vec_log_proba)
  }
  
  vec_all_clust_size <- 1:max_cluster_size # l
  vec_proba_clust_size <- exp(get_vec_log_proba_cluster_size(vec_all_clust_size, R0, k, p)) # r_l
  
  vec_proba_obs_cluster_size <- sapply(vec_cluster_size, FUN = function(j){
    sum(sapply(j:max_cluster_size, FUN = function(l){
      if(p_detect < 1.0) {
        log_proba <- log(vec_proba_clust_size[l]) + lchoose(l, j) + j*log(p_detect) + (l-j)*log(1.0 - p_detect)
        proba <- exp(log_proba)
      }
      else{
        proba <- ifelse(j == l, vec_proba_clust_size[j], 0.)
      }
      
      #vec_proba_clust_size[l] * choose(l, j) * (p_detect ^ j) * ((1. - p_detect) ^ (l - j))
    }))
  })
  proba_obs_zero_element <- sum(sapply(1:max_cluster_size, FUN = function(l){
    vec_proba_clust_size[l] *  (1. - p_detect) ^ l
  }))
  
  return(vec_proba_obs_cluster_size/(1. - proba_obs_zero_element))
}
get_prop_01_offsprings <- function(R0, k){
  prop_indiv_01_offsprings <- dnbinom(x = c(0, 1), size = k, mu = R0)
  
  return(prop_indiv_01_offsprings)
}

max_cluster_size_sim <- 1e4
n_mpox_cases <- 85000
n_seq / n_mpox_cases
p_detect <- 0.055
p_detect_2 <- 0.5

vec_cluster_size <- 1:size_largest_polytomy

vec_R <- seq(0.1, 1.6, 0.01)
vec_log_k <- seq(-4.7, 2.4, 0.1)

df_prop_indiv_R_greater_than_1 <- expand.grid(R = vec_R,
                                              log_k = vec_log_k) %>% 
  mutate(k = exp(log_k))

df_prop_indiv_R_greater_than_1$prop_0_offsprings <- sapply(1:nrow(df_prop_indiv_R_greater_than_1), FUN = function(i_row){
  get_prop_01_offsprings(R0 = df_prop_indiv_R_greater_than_1$R[i_row], k = df_prop_indiv_R_greater_than_1$k[i_row])[1]
})

df_prop_indiv_R_greater_than_1$prop_01_offsprings <- sapply(1:nrow(df_prop_indiv_R_greater_than_1), FUN = function(i_row){
  sum(get_prop_01_offsprings(R0 = df_prop_indiv_R_greater_than_1$R[i_row], k = df_prop_indiv_R_greater_than_1$k[i_row])[1:2])
})

df_proba_obs_large_clust <- expand.grid(R = vec_R,
                                        log_k = vec_log_k) %>% 
  mutate(k = exp(log_k))


n_cores <- 12

cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
t0 <- Sys.time()
df_proba_obs_large_clust$proba_at_least_size_largest_polytomy <-
  Reduce('c',foreach(i_row = 1:nrow(df_proba_obs_large_clust), .packages = c('dplyr')) %dopar% {

  curr_R <- df_proba_obs_large_clust[i_row, 'R']
  curr_k <- df_proba_obs_large_clust[i_row, 'k']

  ## Probability to observe a cluster of vec_cluster_size
  curr_vec_prob <- get_vec_proba_obs_cluster_size(vec_cluster_size, curr_R, curr_k, p_trans_before_mut, p_detect, max_cluster_size_sim)

  proba_cluster_greater_than_cluster <- 1.0 - sum(curr_vec_prob[1:(size_largest_polytomy - 1)])
})

df_proba_obs_large_clust$proba_at_least_size_largest_polytomy_2 <-
  Reduce('c',foreach(i_row = 1:nrow(df_proba_obs_large_clust), .packages = c('dplyr')) %dopar% {
    
    curr_R <- df_proba_obs_large_clust[i_row, 'R']
    curr_k <- df_proba_obs_large_clust[i_row, 'k']
    
    ## Probability to observe a cluster of vec_cluster_size
    curr_vec_prob <- get_vec_proba_obs_cluster_size(vec_cluster_size, curr_R, curr_k, p_trans_before_mut, p_detect_2, max_cluster_size_sim)
    
    proba_cluster_greater_than_cluster <- 1.0 - sum(curr_vec_prob[1:(size_largest_polytomy - 1)])
  })

stopCluster(cl)
t1 <- Sys.time()
print(t1 - t0)

df_proba_obs_large_clust <- df_proba_obs_large_clust %>%
  mutate(proba_obs_cluster_size = 1.0 - exp(n_clust_identical_sequences * log(1.0 - proba_at_least_size_largest_polytomy)),
         proba_obs_cluster_size_2 = 1.0 - exp(n_clust_identical_sequences * log(1.0 - proba_at_least_size_largest_polytomy_2)))

#saveRDS(df_proba_obs_large_clust, '../results/proba_large_polytomy/df_proba_obs_large_clust.rds')
