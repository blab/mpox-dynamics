library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(parallel)
source('utils_inference.R')
source('config_runs.R')

df_inference_joint <- readRDS('../results/estim_dispersion/df_inference_joint_k.rds')

### Joint inference
df_inference_joint <- readRDS('../results/estim_dispersion/df_inference_joint_k.rds')

plt_R_estim_joint_inference <- df_inference_joint %>%
  filter(param == 'R') %>% 
  mutate(location = factor(location, 
                           levels = c('Portugal', 'United Kingdom',
                                      'California', 'New York', 'Washington'))) %>% 
  mutate(is_US = (location %in% c('California', 'New York', 'Washington'))) %>% 
  ggplot(aes(x = as.factor(location), y = mle_estim,
             group = interaction(prop_inf_det, param, location),
             colour = as.factor(prop_inf_det))) +
  geom_rect(data = NULL,
            alpha = 0.1, fill = 'darkgrey', color = NA,
            aes(xmin = -Inf, xmax = Inf, ymin = 0.21, ymax = 0.42)) +
  geom_hline(yintercept = 0.30, linetype = 'dashed', color = 'black') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.4)) +
  scale_x_discrete(name = '') +
  scale_y_continuous(limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05)),
                     name = 'R estimate',
                     breaks = seq(0.0, 1.5, 0.25)) +
  scale_colour_manual(name = 'Proportion of\ninfections\ndetected',
                      breaks = c(0.1, 0.5, 1.0),
                      labels = c('10%', '50%', '100%'),
                      values = c('#b8aebf', '#7d6dc1', '#483d8b')) +
  facet_wrap(. ~ is_US, scales = 'free_x') +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank()) 

plt_k_estim_joint_inference <- df_inference_joint %>%
  filter(param == 'k') %>% 
  ggplot(aes(x = as.factor(location), y = mle_estim,
             group = interaction(prop_inf_det, param, location),
             colour = as.factor(prop_inf_det))) +
  geom_rect(data = NULL,
            alpha = 0.1, fill = 'darkgrey', color = NA,
            aes(xmin = -Inf, xmax = Inf, ymin = 0.21, ymax = 0.42)) +
  geom_hline(yintercept = 0.30, linetype = 'dashed', color = 'black') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.4)) +
  scale_x_discrete(name = '') +
  scale_y_continuous(limits = c(0.01, 10.),
                     trans = 'log',
                     expand = expansion(mult = c(0., 0.0)),
                     breaks = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05,
                                0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                     labels = c('0.001', '0.002', '0.005', '0.01', '0.02', '0.05',
                                '0.1', '0.2', '0.5', '1', '2', '5', '10'),
                     name = 'k estimate') +
  scale_colour_manual(name = 'Proportion of\ninfections\ndetected',
                      breaks = c(0.1, 0.5, 1.0),
                      labels = c('10%', '50%', '100%'),
                      values = c('#b8aebf', '#7d6dc1', '#483d8b')) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank()) 

plt_joint_inference <- ggarrange(plt_R_estim_joint_inference + theme(legend.position = 'none'),
                                 plt_k_estim_joint_inference,
                                 ncol = 2, widths = c(5., 3.0),
                                 labels = c('C', 'D'))

plot(plt_joint_inference)
