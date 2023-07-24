get_proba_size_cluster <- function(cluster_size, R0, k, p){
  
  proba <-
    gamma(k * cluster_size + cluster_size - 1) /
    ( gamma(k * cluster_size) * gamma(cluster_size + 1) ) *
    (p * R0 / k) ^ (cluster_size - 1) * (1 + p * R0 / k) ^ (1 - k * cluster_size - cluster_size) 
  
  return(proba)
}

get_log_proba_size_cluster <- function(cluster_size, R0, k, p){
  
  log_proba <- 
    lgamma(k * cluster_size + cluster_size - 1) - lgamma(k * cluster_size) - lgamma(cluster_size + 1) +
    (cluster_size - 1)*log(p * R0 / k) - (k * cluster_size + cluster_size - 1) * log(1 + p * R0 / k)
  
  return(log_proba)
}

get_vec_proba_cluster_size <- function(vec_cluster_size, R0, k, p){
  vec_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_proba_size_cluster(cluster_size, R0, k, p)
  })
  return(vec_proba)
}

get_vec_log_proba_cluster_size <- function(vec_cluster_size, R0, k, p){
  vec_log_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_log_proba_size_cluster(cluster_size, R0, k, p)
  })
  return(vec_log_proba)
}

get_loglik <- function(vec_cluster_size, vec_n_clusters, R0, k, p){
  loglik <- sum(
    vec_n_clusters * get_vec_log_proba_cluster_size(vec_cluster_size, R0, k, p)
  )
  return(loglik)
}


######## Same as above but accounting for imperfect observation
get_vec_proba_obs_cluster_size <- function(vec_cluster_size, R0, k, p, p_detect, max_cluster_size){
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

get_vec_proba_obs_cluster_size_dependent <- function(vec_cluster_size, R0, k, p, p_sent, max_cluster_size){
  vec_all_clust_size <- 1:max_cluster_size # l
  vec_proba_clust_size <- exp(get_vec_log_proba_cluster_size(vec_all_clust_size, R0, k, p)) # r_l
  
  vec_proba_obs_cluster_size <- sapply(vec_cluster_size, FUN = function(j){
    vec_proba_clust_size[j] * (1. - (1. - p_sent) ^ j)
  })
  proba_obs_zero_element <- sum(sapply(1:max_cluster_size, FUN = function(l){
    vec_proba_clust_size[l] * (1. - p_sent) ^ l
  }))
  
  return(vec_proba_obs_cluster_size/(1. - proba_obs_zero_element))
}

get_loglik_obs <- function(vec_cluster_size, vec_n_clusters, R0, k, p, p_detect, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R0, k, p, p_detect, max_cluster_size))
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  loglik <- sum(
    vec_n_clusters * log_proba_obs
  )
  return(loglik)
}

get_loglik_obs_no_sum <- function(vec_cluster_size, R0, k, p, p_detect, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R0, k, p, p_detect, max_cluster_size))
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  return(log_proba_obs)
}

get_loglik_obs_size_dep <- function(vec_cluster_size, vec_n_clusters, R0, k, p, p_sent, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size_dependent(vec_cluster_size, R0, k, p, p_sent, max_cluster_size))
  
  loglik <- sum(
    vec_n_clusters * log_proba_obs
  )
  return(loglik)
}


######### Running a grid search
run_grid_search <- function(vec_R0, vec_k, p_detect, max_cluster_size_inference,
                            cluster_alloc, p_trans_before_mut){
  df_grid <- expand.grid(R0 = vec_R0, k = vec_k) %>% 
    mutate(p_seq_infection = p_detect)
  
  df_grid$log_lik <- sapply(1:nrow(df_grid), FUN = function(i){
    get_loglik_obs(cluster_alloc$df_dist_clique_size$cluster_size,
                   cluster_alloc$df_dist_clique_size$count,
                   df_grid[i, 'R0'], df_grid[i, 'k'],
                   p = p_trans_before_mut, p_detect = p_detect,
                   max_cluster_size = max_cluster_size_inference)
  })
  
  return(df_grid)
  
}

run_grid_search_from_vector <- function(vec_R0, vec_k, p_detect,
                                        max_cluster_size_inference,
                                        vec_clusters_sizes, vec_counts_cluster_sizes,
                                        p_trans_before_mut){
  df_grid <- expand.grid(R0 = vec_R0, k = vec_k) %>% 
    mutate(p_seq_infection = p_detect)
  
  df_grid$log_lik <- sapply(1:nrow(df_grid), FUN = function(i){
    get_loglik_obs(vec_clusters_sizes,
                   vec_counts_cluster_sizes,
                   df_grid[i, 'R0'], df_grid[i, 'k'],
                   p = p_trans_before_mut, p_detect = p_detect,
                   max_cluster_size = max_cluster_size_inference)
  })
  
  return(df_grid)
  
}

######### Obtaining likelihood profile CI from the results of a grid search
get_profile_likelihood_CI <- function(df_grid, col_loglik = 'log_lik_full'){
  # https://web.stat.tamu.edu/~suhasini/teaching613/chapter3.pdf
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  df_profile_R0 <- df_grid %>%
    rename(loglik = col_loglik) %>%
    group_by(R0) %>% 
    summarise(loglik_profile = max(loglik)) %>% 
    mutate(max_loglik = max(loglik_profile)) %>% 
    mutate(LR = 2 * (max_loglik - loglik_profile),
           is_within_95 = (LR < chisq_stat_95),
           is_within_90 = (LR < chisq_stat_90),
           is_within_50 = (LR < chisq_stat_50))
  
  df_profile_k <- df_grid %>%
    rename(loglik = col_loglik) %>%
    group_by(k) %>% 
    summarise(loglik_profile = max(loglik)) %>% 
    mutate(max_loglik = max(loglik_profile)) %>% 
    mutate(LR = 2 * (max_loglik - loglik_profile),
           is_within_95 = (LR < chisq_stat_95),
           is_within_90 = (LR < chisq_stat_90),
           is_within_50 = (LR < chisq_stat_50))
  
  CI_R0 <- df_profile_R0 %>% 
    summarise(mle_estim = R0[loglik_profile == max_loglik],
              lower_95 = min(R0[is_within_95]),
              upper_95 = max(R0[is_within_95]),
              lower_90 = min(R0[is_within_90]),
              upper_90 = max(R0[is_within_90]),
              lower_50 = min(R0[is_within_50]),
              upper_50 = max(R0[is_within_50])) %>% 
    mutate(param = 'R0')
  
  CI_k <- df_profile_k %>% 
    summarise(mle_estim = k[loglik_profile == max_loglik],
              lower_95 = min(k[is_within_95]),
              upper_95 = max(k[is_within_95]),
              lower_90 = min(k[is_within_90]),
              upper_90 = max(k[is_within_90]),
              lower_50 = min(k[is_within_50]),
              upper_50 = max(k[is_within_50])) %>% 
    mutate(param = 'k')
  
  df_CI <- bind_rows(CI_R0, CI_k)
  return(df_CI)
}


get_profile_likelihood_CI_multiple_params <- function(df_grid, 
                                                      names_param = c('R0', 'k'),
                                                      col_loglik = 'loglik'){
  
  # https://web.stat.tamu.edu/~suhasini/teaching613/chapter3.pdf
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  df_all_CI <- Reduce('bind_rows', lapply(1:length(names_param), FUN = function(i_param){
    df_profile <- df_grid %>% 
      rename(loglik = col_loglik) %>% 
      rename(curr_param = names_param[i_param]) %>% 
      group_by(curr_param, p_detect_cases) %>% 
      summarise(loglik_profile = max(loglik)) %>% 
      group_by(p_detect_cases) %>% 
      mutate(max_loglik = max(loglik_profile)) %>% 
      mutate(LR = 2 * (max_loglik - loglik_profile),
             is_within_95 = (LR < chisq_stat_95),
             is_within_90 = (LR < chisq_stat_90),
             is_within_50 = (LR < chisq_stat_50))
    
    df_CI <- df_profile %>% 
      group_by(p_detect_cases) %>% 
      summarise(mle_estim = curr_param[loglik_profile == max_loglik],
                lower_95 = min(curr_param[is_within_95]),
                upper_95 = max(curr_param[is_within_95]),
                lower_90 = min(curr_param[is_within_90]),
                upper_90 = max(curr_param[is_within_90]),
                lower_50 = min(curr_param[is_within_50]),
                upper_50 = max(curr_param[is_within_50])) %>% 
      mutate(param = names_param[i_param])
    
    df_CI
  }))
  
  return(df_all_CI)
}



get_CI <- function(i_param,
                   minus_log_lik, mle_estim,
                   vec_param){
  
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)

  
  mle_param <- mle_estim$par
  
  ## Log-likelihood profile
  vec_log_lik <- sapply(vec_param, FUN = function(curr_param){
    mle_param[i_param] <- curr_param
    - minus_log_lik(mle_param)
  })
  
  ## Likelihood ratio
  vec_LR <- 2 * (- mle_estim$value - vec_log_lik)
  
  param_within_95 <- vec_param[vec_LR < chisq_stat_95]
  param_within_90 <- vec_param[vec_LR < chisq_stat_90]
  param_within_50 <- vec_param[vec_LR < chisq_stat_50]
  
  return(
    c('lower_95' = min(param_within_95),
      'upper_95' = max(param_within_95),
      'lower_90' = min(param_within_90),
      'upper_90' = max(param_within_90),
      'lower_50' = min(param_within_50),
      'upper_50' = max(param_within_50))
  )
  
  
}

