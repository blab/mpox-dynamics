library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(parallel)
source('utils_inference.R')
source('config_runs.R')

## Build the model
create_full_model <- function(list_config, prop_inf_det = 1.0, max_cluster_size_inference = 10000){
  
  n_locations <- length(list_config)
  
  ## Probability that transmission occurs before mutation for mpox 
  p_trans_before_mut <- 0.66259
  
  ## Setting up the bounds for the grid search
  R_min <- 0.01
  R_max <- 10.0
  k_min <- 0.001
  k_max <- 10.0
  
  ## Proportion of infections sequenced by location
  vec_prop_inf_sequenced <- sapply(1:length(list_config), FUN = function(i_location){
    curr_config <- list_config[[i_location]]
    prop_inf_det * curr_config$n_seq_during_period / curr_config$n_cases_during_period
  })
  
  ## 
  minus_log_lik_location <- function(i_location, R, k){
    curr_config <- list_config[[i_location]]
    curr_prop_infections_sequenced <- vec_prop_inf_sequenced[i_location]
    
    ll <- get_loglik_obs(vec_cluster_size = curr_config$df_clust_dist$n_in_clique,
                         vec_n_clusters = curr_config$df_clust_dist$n_clusters,
                         R0 = R, k = k,
                         p = p_trans_before_mut,
                         p_detect = curr_prop_infections_sequenced,
                         max_cluster_size = max_cluster_size_inference)
    
    return(-ll)
  }
  
  minus_log_lik_joint <- function(param){
    vec_R <- param[1:n_locations]
    k <- param[n_locations + 1]
    
    minus_ll <- Reduce('+', lapply(1:n_locations, FUN = function(i_location){
      minus_log_lik_location(i_location, vec_R[i_location], k)
    }))
    return(minus_ll)
  }
  
  get_mle_joint <- function(){
    mle_estim <- optim(par = c(rep(1.0, n_locations), 0.1),
                       fn = minus_log_lik_joint,
                       method = 'L-BFGS-B',
                       lower = c(rep(R_min, n_locations), k_min), 
                       upper = c(rep(R_max, n_locations), k_max))
    
    
    return(mle_estim)
  }
  
  return(list(
    n_locations = n_locations,
    p_trans_before_mut = p_trans_before_mut,
    minus_log_lik_location = minus_log_lik_location,
    minus_log_lik_joint = minus_log_lik_joint,
    get_mle_joint = get_mle_joint
  ))
}

## Configuration for the different countries
vec_name_locations <- c('Portugal', 'United Kingdom', 'Washington', 'California', 'New York')

list_config <- lapply(vec_name_locations, get_config)

## Models with different assumptions regarding the fraction of infections detected
full_mod_10percent <- create_full_model(list_config, prop_inf_det = 0.1, max_cluster_size_inference = 10000)
full_mod_50percent <- create_full_model(list_config, prop_inf_det = 0.5, max_cluster_size_inference = 10000)
full_mod_100percent <- create_full_model(list_config, prop_inf_det = 1.0, max_cluster_size_inference = 10000)

## MLEs for the different models with CIs
vec_R_CI <- seq(0.1, 1.65, 0.01) # Vectors with the resolution used to compute the likelihood profiles CIs
vec_k_CI <- c(seq(0.001, 0.1, 0.001), seq(0.11, 1.0, 0.01), seq(1.1, 10.0, 0.1))

list_mod <- list(full_mod_10percent, full_mod_50percent, full_mod_100percent)

vec_p_seq_inf <- c(0.1, 0.5, 1.0)
n_cores <- 3
cl <- makeForkCluster(n_cores)
t0 <- Sys.time()
df_inference_joint <- Reduce('bind_rows', parLapply(cl, 1:length(list_mod), fun = function(i_mod){
  curr_mod <- list_mod[[i_mod]]
  n_locations <- curr_mod$n_locations
  
  ## MLE
  curr_mle <- curr_mod$get_mle_joint()
  
  ## Confidence intervals
  CI_k <- get_CI(n_locations + 1, curr_mod$minus_log_lik_joint, curr_mle, vec_param = vec_k_CI)
  
  mat_CI_R <- Reduce('bind_rows', lapply(1:n_locations, FUN = function(i_location){
    CI_R <- get_CI(i_location, curr_mod$minus_log_lik_joint, curr_mle, vec_param = vec_R_CI)
  }))
  
  tibble(param = c(rep('R', n_locations), 'k'),
         location = c(vec_name_locations, 'all'),
         mle_estim = curr_mle$par) %>% 
    bind_cols(bind_rows(mat_CI_R, CI_k)) %>% 
    mutate(prop_inf_det = vec_p_seq_inf[i_mod])
}))
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

#saveRDS(df_inference_joint, '../results/estim_dispersion/df_inference_joint_k.rds')

