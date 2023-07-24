get_df_timing_clusters <- function(cluster_alloc_timing){
  df_timing <- cluster_alloc_timing$df_clusters_timing %>% 
    mutate(year = substr(date, 1, 4),
           month = substr(date, 6, 7)) %>% 
    group_by(clique_id) %>% 
    #mutate(n_month_NA_within_clique = sum(month == 'XX')) %>% 
    #filter(n_month_NA_within_clique == 0) %>% 
    mutate(imputed_month_for_min = case_when(month == 'XX' ~ '13',
                                             TRUE ~ month),
           imputed_month_for_max = case_when(month == 'XX' ~ '00',
                                             TRUE ~ month)) %>% 
    group_by(clique_id) %>% 
    summarise(min_month = min(as.numeric(paste0(year, imputed_month_for_min))),
              max_month = max(as.numeric(paste0(year, imputed_month_for_max))),
              n_in_clique = n()) %>% 
    group_by(n_in_clique, min_month, max_month) %>% 
    summarise(n_clusters = n())
  
  return(df_timing)
}

get_cluster_distrib_during_period <- function(df_timing_clusters, vec_month_beginning, month_end){
  df_clust_dist <- df_timing_clusters %>% 
    filter(min_month %in% vec_month_beginning, max_month <= month_end) %>% 
    group_by(n_in_clique) %>% 
    summarise(n_clusters = sum(n_clusters)) 
  
  return(df_clust_dist)
}

get_cluster_distrib_over_time <- function(df_timing_clusters,
                                          vec_month_beginning = c(202205, 202206, 202207, 202208, 202209, 202210, 202211, 202212)){
  
  df_cluster_distrib_over_time <- Reduce('bind_rows', lapply(vec_month_beginning, FUN = function(curr_month_beginning){
    df_timing_clusters %>% 
      filter(min_month == curr_month_beginning) %>% 
      group_by(n_in_clique) %>% 
      summarise(n_clusters = sum(n_clusters)) %>% 
      mutate(month = curr_month_beginning)
  }))
}

get_config <- function(name_location){
  
  file_path_clean_owid <- '../data/owid-monkeypox.rds' # File path for cases from OWID
  
  if(name_location == 'World'){
    file_path_cluster_alloc <- paste0('../data/cluster_alloc_timing.rds')
    file_path_metadata <- paste0('../data/metadata_seq.rds')
  } else{
    file_path_cluster_alloc <- paste0('../data/cluster_alloc_timing_', name_location, '.rds') # File path for cluster alloc file
    file_path_metadata <- paste0('../data/metadata_for_analysis_', name_location, '.rds') # File path for sequence metadata
  }
  
  vec_month_beginning <- c(202208, 202209, 202210, 202211, 202212)
  month_end <- 202303
  
  ## Load data
  cluster_alloc_timing <- readRDS(file_path_cluster_alloc)
  df_timing_clusters <- get_df_timing_clusters(cluster_alloc_timing)
  
  ## Cluster distribution over time
  df_cluster_distrib_over_time <- get_cluster_distrib_over_time(df_timing_clusters)
  
  df_seq_per_month <- readRDS(file_path_metadata) %>%  filter(month != 'XX') %>% 
    group_by(month, year) %>% summarise(n_seq = n()) %>% 
    mutate(full_month = as.numeric(paste0(year, month))) %>% 
    ungroup()
  
  ## Load case data
  if(name_location %in% c( 'California', 'New York')){
    file_path_state_case <- paste0('../data/daily_cases_', name_location, '.rds')
    df_cases_per_month <- readRDS(file_path_state_case) %>% 
      mutate(date_char = as.character(date),
             year = substr(date_char, 1, 4), month = substr(date_char, 6, 7),
             full_month = as.numeric(paste0(year, month))) %>% 
      rename(new_cases = n_cases) %>% 
      select(full_month, new_cases) %>% 
      group_by(full_month) %>% 
      summarise(new_cases = sum(new_cases))
    
  } else if(name_location == 'Washington'){
    file_path_state_case <- paste0('../data/weekly_cases_', name_location, '.rds')
    df_cases_per_month <- readRDS(file_path_state_case) %>% 
      mutate(date_char = as.character(date_week),
             year = substr(date_char, 1, 4), month = substr(date_char, 6, 7),
             full_month = as.numeric(paste0(year, month))) %>% 
      rename(new_cases = n_cases) %>% 
      select(full_month, new_cases) %>% 
      group_by(full_month) %>% 
      summarise(new_cases = sum(new_cases))
  } else{
    df_cases_per_month <- readRDS(file_path_clean_owid) %>%
      filter(country == name_location) %>% 
      mutate(date_char = as.character(date),
             year = substr(date_char, 1, 4), month = substr(date_char, 6, 7),
             full_month = as.numeric(paste0(year, month))) %>% 
      select(full_month, new_cases) %>% 
      group_by(full_month) %>% 
      summarise(new_cases = sum(new_cases))
  }
  
  ## 
  df_clust_dist <- get_cluster_distrib_during_period(df_timing_clusters, vec_month_beginning, month_end)
  
  n_seq_during_period <- df_seq_per_month %>%
    filter(full_month >= min(vec_month_beginning), full_month <= month_end) %>% 
    summarise(n_seq = sum(n_seq)) %>% as.numeric()
  
  n_cases_during_period <- df_cases_per_month %>%
    filter(full_month >= min(vec_month_beginning), full_month <= month_end) %>% 
    summarise(n_cases = sum(new_cases)) %>% as.numeric()
  
  return(
    list('df_clust_dist' = df_clust_dist,
         'n_seq_during_period' = n_seq_during_period,
         'n_cases_during_period' = n_cases_during_period,
         'df_cluster_distrib_over_time' = df_cluster_distrib_over_time)
  )
}
  
  
