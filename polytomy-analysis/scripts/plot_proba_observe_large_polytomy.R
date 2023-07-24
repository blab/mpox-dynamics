library(dplyr)
library(RColorBrewer)

df_proba_obs_large_clust <- readRDS('../results/proba_large_polytomy/df_proba_obs_large_clust.rds')


plot_proba_large_clust <- function(df_proba_obs_large_clust){
  vec_log_k <- seq(-4.7, 2.4, 0.1)
  
  breaks_k <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
  log_breaks_k <- sapply(log(breaks_k), FUN = function(curr_log_break){
    vec_log_k[vec_log_k > curr_log_break][1]
  })
  
  ## Previous estimates from the analysis of epidemiological cluster
  ## during previous outbreaks (Blumberg and Lloyd-Smith, 2013)
  R_mle <- 0.30
  R_lower <- 0.21
  R_upper <- 0.42
  k_mle <- 0.36
  k_lower <- 0.14
  k_upper <- 2.57
  
  df_proba_obs_large_clust %>% 
    mutate(k = exp(log_k)) %>% 
    filter(k > 0.09, k < 0.12, proba_obs_cluster_size > 0.05) %>% 
    group_by(k) %>% 
    filter(proba_obs_cluster_size == min(proba_obs_cluster_size)) %>% 
    print()
  
  plt_proba <- df_proba_obs_large_clust %>% 
    mutate(crop_prob = ifelse(proba_obs_cluster_size < 1e-8, 1e-8, proba_obs_cluster_size),
           log10_prob = log10(crop_prob)) %>% 
    ggplot(aes(y = log_k, x = R)) +
    geom_raster(aes(fill = log10_prob)) +
    geom_contour(aes(z = log10_prob),
                 breaks = c(-4, -2, -1),
                 col = 'white', linetype = 'dashed') +
    scale_fill_gradientn(name = paste0('Probability to observe a polytomy of size ', size_largest_polytomy),
                         colours = brewer.pal(9, 'OrRd'),
                         breaks = c(-8, -6, -4, -2, 0),
                         labels = c(expression(paste('\u2264 ', 10^{-8})), expression(10^{-6}),
                                    expression(10^{-4}), expression(10^{-2}), expression(10^{0}))) +
    scale_y_continuous(name = 'Dispersion parameter k',
                       breaks = log_breaks_k,
                       labels = breaks_k,
                       expand = expansion(mult = c(0., 0.))) +
    scale_x_continuous(name = 'Reproduction number R',
                       breaks = seq(0.1, 1.6, 0.2),
                       expand = expansion(mult = c(0., 0.))) +
    theme_classic() +
    coord_cartesian(ylim = c(log(0.009), log(4))) +
    theme(axis.line = element_blank(),
          legend.position = 'top',
          legend.key.width = unit(1.5, "cm")) +
    geom_point(aes(x = R_mle, y = log(k_mle))) +
    geom_segment(aes(x = R_lower, xend = R_upper,
                     y = log(k_mle), yend = log(k_mle))) +
    geom_segment(aes(x = R_mle, xend = R_mle,
                     y = log(k_lower), yend = log(k_upper))) +
    guides(fill = guide_colorbar(title.position = 'top',
                                 title.hjust = 0.5)) +
    annotate(geom = 'text', colour = 'white', angle = 30, 
             size = 2,
             x = 1.03, y = log(0.012), label = expression(10^{-1})) +
    annotate(geom = 'text', colour = 'white', angle = 30, 
             size = 2,
             x = 0.92, y = log(0.01205), label = expression(10^{-2})) +
    annotate(geom = 'text', colour = 'white', angle = 30, 
             size = 2,
             x = 0.83, y = log(0.014), label = expression(10^{-4}))
  
  return(plt_proba)
}

size_largest_polytomy <- 118

plot_proba_large_clust(df_proba_obs_large_clust)
