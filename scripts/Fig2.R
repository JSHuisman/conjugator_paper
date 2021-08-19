###########################################################
# Conjugation rate estimation simulations
# Author: Jana Huisman
###########################################################

library(conjugator)

library(deSolve)
library(ggplot2)
library(reshape2)
#library(ggrepel)
#library(gridExtra) # to arrange multiple ggplots

library(tidyverse) # config df
library(patchwork)

source('./R/model_functions.R')
source('./R/plotting_functions.R')
source('./R/utils.R')

col_palet <- c('#332288','#88CCEE',
               '#44AA99','#117733',
               '#999933','#DDCC77',
               '#CC6677','#882255',
               '#AA4499','#DDDDDD',
               '#000000')

###########################################################
## Varying Dsize ----
Dsize_vec = 10**seq(4, 10, 0.2)
vars <- cbind(r=Dsize_vec, t=rep(0, length(Dsize_vec)),
              d=Dsize_vec, c = rep(1E14, length(Dsize_vec)))

# equal rates
growth_rates = c(psi.D = 1, psi.R = 1, psi.T = 1)
gammas = c(gamma.t = 10^-13, gamma.d = 10^-13)
config_df_eq <- create_config_df(vars, growth_rates, gammas)
# add high donor growth
growth_rates = c(psi.D = 2, psi.R = 1, psi.T = 1)
gammas = c(gamma.t = 10^-13, gamma.d = 10^-13)
config_df <- bind_rows(config_df_eq, create_config_df(vars, growth_rates, gammas))

## compute outputs
t.vec <- seq(0, 8, length = 81)
model_out <- compute_model_scan(config_df, t.vec, output_timepoints = c(21, 41, 81))

full_out <- compute_estimates(model_out$full_out)

test <- full_out[,c('t', 'D.0', 'TD', 'psi.D')]
molten_test <- melt(test, id.vars = c('D.0', 't', 'psi.D'))
molten_test$time <- as.factor(molten_test$t)

Dsize_fc_plot <- ggplot(molten_test) +
  geom_point(aes(x = D.0,
                 y = value,
                 shape = factor(t),
                 colour = factor(psi.D))) +
  scale_y_continuous(trans = 'log', breaks = .base_breaks()) +
  scale_x_continuous(trans = 'log', breaks = .base_breaks()) +
  labs(x = bquote('Initial population size'),
       y = bquote('T/D'),
       shape = 'Time (h)',
       colour = 'Donor\ngrowth\nrate (1/h)') +
  scale_color_manual(values = col_palet[c(3,9)]) +
  theme_minimal() + theme(text = element_text(size=15))

Dsize_fc_plot
ggsave(filename = paste0('../figures/',
                         'TD_afo_Dsize_time_log.png'),
       plot = Dsize_fc_plot, width = 6, height = 4)

###########################################################
## Varying D:R ----
# [9:1, 4:1, 1:1, 1:4, 1:9]
# total = 1E7
DRratio_vec = 9:1 #c(9, 8, 5, 2, 1)
vars <- cbind(r=1E6*(10-DRratio_vec), t=rep(0, 9),
              d=1E6*DRratio_vec, c = rep(1E14, 9))

# equal rates
growth_rates = c(psi.D = 1, psi.R = 1, psi.T = 1)
gammas = c(gamma.t = 10^-13, gamma.d = 10^-13)
config_df_eq <- create_config_df(vars, growth_rates, gammas)
# add high donor growth
growth_rates = c(psi.D = 2, psi.R = 1, psi.T = 1)
gammas = c(gamma.t = 10^-13, gamma.d = 10^-13)
config_df <- bind_rows(config_df_eq, create_config_df(vars, growth_rates, gammas))

## compute outputs
t.vec <- seq(0, 8, length = 81)
model_out <- compute_model_scan(config_df, t.vec, output_timepoints = c(21, 41, 81))
#apply(model_out$crit_times, 1, min)
full_out <- compute_estimates(model_out$full_out)

test <- full_out[,c('t', 'D.0', 'TD', 'psi.D')]
molten_test <- melt(test, id.vars = c('D.0', 't', 'psi.D'))
molten_test$time <- as.factor(molten_test$t)

DRratio_fc_plot <- ggplot(molten_test) +
  geom_point(aes(x = D.0,
                 y = value,
                 shape = factor(t),
                 colour = factor(psi.D))) +
  scale_x_continuous(breaks = 1E6*DRratio_vec,
                     labels = c('9:1', '8:2', '7:3', '6:4',
                                '1:1', '4:6', '3:7', '1:8', '1:9')) +
  scale_y_continuous(trans = 'log', breaks = .base_breaks()) +
  labs(x = 'Initial donor to recipient ratio',
       y = 'T/D',
       shape = 'Time (h)',
       colour = 'Donor\ngrowth\nrate (1/h)') +
  scale_color_manual(values = col_palet[c(3,9)]) +
  theme_minimal() + theme(text = element_text(size=15))

DRratio_fc_plot
ggsave(filename = paste0('../figures/',
                         'TD_afo_DRratio_time_log.png'),
       plot = DRratio_fc_plot, width = 6, height = 4)

###########################################################

Dsize_fc_plot + DRratio_fc_plot +
  plot_layout(ncol = 2, guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'right')

ggsave(filename = '../figures/TD_afo_DR_log.png', 
       width = 10, height = 4)
