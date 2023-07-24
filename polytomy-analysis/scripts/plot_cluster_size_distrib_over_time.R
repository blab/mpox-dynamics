library(dplyr)
library(RColorBrewer)
source('utils_inference.R')
source('config_runs.R')

## Mean cluster size over time for the locations of interest
df_mean_clust_size <- Reduce('bind_rows', lapply(c('World', vec_name_locations), FUN = function(curr_location){
  
  get_config(curr_location)$df_cluster_distrib_over_time %>% 
    group_by(month) %>% 
    summarise(mean_clust_size = sum(n_in_clique * n_clusters) / sum(n_clusters),
              n_clusters = sum(n_clusters)) %>% 
    mutate(location = curr_location)
  
}))

col_world <- 'brown'
col_portugal <- '#EEA160'
col_UK <- '#356D4C'
col_CA <- '#aad7ff'
col_NY <- '#2664a5'
col_WA <- '#133253'

plt_mean_cluster_size <- df_mean_clust_size %>% 
  filter(n_clusters >= 10) %>% 
  ggplot(aes(x = as.factor(month), y = mean_clust_size, 
             group = location, colour = as.factor(location))) +
  geom_point() + geom_line() +
  scale_x_discrete(name = 'Month of first cluster detection',
                   breaks = 202205:202212,
                   labels = c('May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')) +
  scale_y_continuous(name = 'Mean size of clusters\nof identical sequences') +
  scale_colour_manual(name = 'Location',
                      breaks = c('World','Portugal', 'United Kingdom', 
                                 'California', 'New York', 'Washington'),
                      labels = c('World','Portugal', 'United Kingdom', 
                                 'California', 'New York', 'Washington'),
                      values = c(col_world, col_portugal, col_UK, col_CA, col_NY, col_WA)) +
  theme_classic()
