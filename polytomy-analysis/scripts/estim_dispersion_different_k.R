library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(parallel)
source('utils_inference.R')
source('config_runs.R')

## Build the model
create_model <- function(curr_config, prop_inf_det = 1.0, max_cluster_size_inference = 10000){
  
  ## Probability that transmission occurs before mutation for mpox 
  p_trans_before_mut <- 0.66259
  
  ## Setting up the bounds for the grid search
  R_min <- 0.01
  R_max <- 10.0
  k_min <- 0.001
  k_max <- 10.0
  
  ## Proportion of infections sequenced
  prop_inf_sequenced <- prop_inf_det * curr_config$n_seq_during_period / curr_config$n_cases_during_period
  
  minus_log_lik <- function(R, k){
    
    ll <- get_loglik_obs(vec_cluster_size = curr_config$df_clust_dist$n_in_clique,
                         vec_n_clusters = curr_config$df_clust_dist$n_clusters,
                         R0 = R, k = k,
                         p = p_trans_before_mut,
                         p_detect = prop_inf_sequenced,
                         max_cluster_size = max_cluster_size_inference)
    
    return(-ll)
  }
  
  minus_log_lik_country <- function(param){
    R <- param[1]
    k <- param[2]
    
    minus_ll <- minus_log_lik(R, k)
    
    return(minus_ll)
  }
  
  get_mle_location <- function(){
    mle_estim <- optim(par = c(1.0, 0.1),
                       fn = minus_log_lik_country,
                       method = 'L-BFGS-B',
                       lower = c(R_min, k_min), 
                       upper = c(R_max, k_max))
    
    
    return(mle_estim)
  }
  
  return(list(
    p_trans_before_mut = p_trans_before_mut,
    minus_log_lik_country = minus_log_lik_country,
    get_mle_location = get_mle_location
  ))
}

## Configuration for the different countries
vec_name_locations <- c('Portugal', 'United Kingdom', 'Washington', 'California', 'New York')

list_config <- lapply(vec_name_locations, get_config)

## Models with different assumptions regarding the fraction of infections detected
list_mod_10percent <- lapply(list_config, FUN = function(curr_config){
  create_model(curr_config, prop_inf_det = 0.1, max_cluster_size_inference = 10000)
})
list_mod_50percent <- lapply(list_config, FUN = function(curr_config){
  create_model(curr_config, prop_inf_det = 0.5, max_cluster_size_inference = 10000)
})
list_mod_100percent <- lapply(list_config, FUN = function(curr_config){
  create_model(curr_config, prop_inf_det = 1.0, max_cluster_size_inference = 10000)
})

## MLEs for the different models with CIs
vec_R_CI <- seq(0.1, 1.65, 0.01) # Vectors with the resolution used to compute the likelihood profiles CIs
vec_k_CI <- c(seq(0.001, 0.1, 0.001), seq(0.11, 1.0, 0.01), seq(1.1, 10.0, 0.1))

n_cores <- 5
cl <- makeForkCluster(n_cores)
t0 <- Sys.time()
df_inference_10percent <- Reduce('bind_rows', parLapply(cl, 1:length(list_mod_10percent), fun = function(i_location){
  curr_mod <- list_mod_10percent[[i_location]]
  curr_location <- vec_name_locations[i_location]
  ## MLE
  curr_mle <- curr_mod$get_mle_location()
  
  ## Confidence intervals
  CI_R <- get_CI(1, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_R_CI)
  CI_k <- get_CI(2, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_k_CI)
  
  
  tibble(param = c('R', 'k'),
         location = rep(curr_location, 2),
         mle_estim = curr_mle$par) %>% 
    bind_cols(bind_rows(CI_R, CI_k)) %>% 
    mutate(prop_inf_det = 0.1)
}))
df_inference_50percent <- Reduce('bind_rows', parLapply(cl, 1:length(list_mod_50percent), fun = function(i_location){
  curr_mod <- list_mod_50percent[[i_location]]
  curr_location <- vec_name_locations[i_location]
  ## MLE
  curr_mle <- curr_mod$get_mle_location()
  
  ## Confidence intervals
  CI_R <- get_CI(1, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_R_CI)
  CI_k <- get_CI(2, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_k_CI)
  
  
  tibble(param = c('R', 'k'),
         location = rep(curr_location, 2),
         mle_estim = curr_mle$par) %>% 
    bind_cols(bind_rows(CI_R, CI_k)) %>% 
    mutate(prop_inf_det = 0.5)
}))
df_inference_100percent <- Reduce('bind_rows', parLapply(cl, 1:length(list_mod_100percent), fun = function(i_location){
  curr_mod <- list_mod_100percent[[i_location]]
  curr_location <- vec_name_locations[i_location]
  ## MLE
  curr_mle <- curr_mod$get_mle_location()
  
  ## Confidence intervals
  CI_R <- get_CI(1, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_R_CI)
  CI_k <- get_CI(2, curr_mod$minus_log_lik_country, curr_mle, vec_param = vec_k_CI)
  
  
  tibble(param = c('R', 'k'),
         location = rep(curr_location, 2),
         mle_estim = curr_mle$par) %>% 
    bind_cols(bind_rows(CI_R, CI_k)) %>% 
    mutate(prop_inf_det = 1.0)
}))
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

df_inference_by_location <- bind_rows(df_inference_10percent, df_inference_50percent, df_inference_100percent)
#saveRDS(df_inference_by_location, '../results/estim_dispersion/df_inference_different_k.rds')


get_prop_0_offspring <- function(R0, k){
  prop_indiv_0_offspring <- dnbinom(x = c(0), size = k, mu = R0)
  return(prop_indiv_0_offspring)
}


df_inference_by_location %>% 
  mutate(location = factor(location, levels = c('Portugal', 'United Kingdom', 'California', 'New York', 'Washington'))) %>% 
  mutate(param_char = paste0(round(mle_estim, 3), ' (', lower_95, ' - ', upper_95, ')')) %>% 
  select(param, location, prop_inf_det, param_char, mle_estim) %>% 
  pivot_wider(names_from = param, values_from = c(param_char, mle_estim)) %>% 
  mutate(prop_0 = get_prop_0_offspring(mle_estim_R, mle_estim_k)*100) %>% 
  select(location, prop_inf_det, param_char_R, param_char_k, prop_0) %>% 
  arrange(location, prop_inf_det)
  
