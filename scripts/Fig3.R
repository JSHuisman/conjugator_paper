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

source('../R/model_functions.R')
source('../R/plotting_functions.R')
source('../R/utils.R')

###########################################################
# Fold change - growth rate ----
vars <- c(r=5E6, t=0, d=5E6, c = 1E14) # initial population size
growth_rates = cbind(psi.D = seq(0.5, 1.5, 0.01),
                     psi.R = rep(1, 101), psi.T = rep(1, 101))
gammas = c(gamma.t = 10^-13, gamma.d = 10^-13)
config_df <- create_config_df(vars, growth_rates, gammas)
#
t.vec <- seq(0, 8, length = 81)
#
model_out_gr <- compute_model_scan(config_df, t.vec, output_timepoints = c(11, 81))
full_out_gr <- compute_estimates(model_out_gr$full_out)
full_out_gr <- compute_foldchange(full_out_gr)

test_gr <- full_out_gr[,c('t', 'psi.D', 'gamma_max_fc', 'gammaD_fc')]
molten_test_gr <- melt(test_gr, id.vars = c('psi.D', 't'))
molten_test_gr$time <- as.factor(molten_test_gr$t)

growth_rate_fc_plot <- ggplot(molten_test_gr) +
  geom_point(aes(x = psi.D,
                 y = value,
                 shape = time,
                 color = variable)) +
  #scale_y_continuous(trans = 'log') +
  scale_y_continuous(limits = c(0.75, 2.5)) +
  scale_color_manual(name = 'Method',
                     breaks = c('gamma_max_fc', 'gammaD_fc'),
                     labels = c(bquote(gamma[max]), bquote(gamma[D])),
                     values = c('lightblue', 'blue')) +
  scale_shape_manual(breaks = c(1, 8),
                     labels = c('1', '8'),
                     values = c(16, 3)) +
  labs(x = bquote('Donor growth rate'~psi[Dmax]),
       y = bquote('Conjugation rate foldchange'),
       shape = 'Time (h)',
       color = 'Method') +
  theme_minimal() + theme(text = element_text(size=20))

growth_rate_fc_plot

ggsave(filename = paste0('../../figures/',
                         'Growth_rate_fc.png'),
       plot = growth_rate_fc_plot, width = 6, height = 4)



growth_rate_crittimes_plot <- ggplot(NULL) +
  geom_point(aes(x = model_out_gr$full_out$psi.D,
                 y = model_out_gr$crit_times[['min_tcrit']])) +
  labs(x = bquote('Donor growth rate'~psi[Dmax]),
       y = 'Critical time (h)') +
  scale_y_continuous(limits = c(6, 13)) +
  theme_minimal() + theme(text = element_text(size=20))

growth_rate_crittimes_plot
ggsave(filename = paste0('../../figures/',
                         'Growth_rate_crit_times.png'),
       plot = growth_rate_crittimes_plot, width = 3, height = 4)

# higher initial densities, more fold-change impact of
# low Donor growth; strong time-dependence of measurement


# reduce the number of x-axis options and make boxplot
subset_out_gr <- full_out_gr[full_out_gr$psi.D %in% seq(0.5, 1.5, 0.2), ]
boxplot_scan(subset_out_gr, 'psi.D', 'TD', 't',
             plotfile = '../../figures/TD_afo_psiD_time_log.png',
             logplot = TRUE, add_means = TRUE)

######################################
# Fold change - conjugation rate ----
vars <- c(r=5E6, t=0, d=5E6, c = 1E14) # initial population size
growth_rates = c(psi.D = 1, psi.R = 1, psi.T = 1)

#gammas = cbind(gamma.t = 10**seq(-13, -8, 0.1), gamma.d = rep(10^-11, 51))
gammas = cbind(gamma.t = 10**seq(-15, -10, 0.1), gamma.d = rep(10^-13, 51))
config_df <- create_config_df(vars, growth_rates, gammas)
#
t.vec <- seq(0, 8, length = 81)
#
model_out <- compute_model_scan(config_df, t.vec, output_timepoints = c(11, 41, 61, 81))
full_out <- compute_estimates(model_out$full_out)
full_out <- compute_foldchange(full_out)

test <- full_out[,c('t', 'gamma.T', 'gamma_max_fc', 'gammaD_fc')]
test['crit_time'] <- model_out$crit_times[['min_tcrit']]
molten_test <- melt(test, id.vars = c('gamma.T', 't', 'crit_time'))
molten_test$time <- as.factor(molten_test$t)

conj_rate_fc_plot <- ggplot(molten_test, aes(x = gamma.T,
                                             y = value)) +
  geom_point(aes(x = gamma.T,
                 y = value,
                 shape = factor(t),
                 color = variable), position = position_jitter()) +
  #scale_y_continuous(trans = 'log', breaks = .base_breaks()) +
  scale_y_continuous(limits = c(0.75, 2.5)) +
  scale_x_continuous(trans = 'log', breaks = .base_breaks()) +
  scale_color_manual(name = 'Method',
                     breaks = c('gamma_max_fc', 'gammaD_fc'),
                     labels = c(bquote(gamma[max]), bquote(gamma[Dmax])),
                     values = c('lightblue', 'blue')) +
  labs(x = bquote('Transconjugant conjugation rate'~gamma[Tmax]),
       y = bquote('Conjugation rate foldchange'),
       shape = 'Time (h)',
       color = 'Method') +
  theme_minimal() + theme(text = element_text(size=20))

conj_rate_fc_plot

ggsave(filename = paste0('../../figures/',
                         #'conj_rate_fc_high_gamma.png'),
                         'conj_rate_fc.png'),
       plot = conj_rate_fc_plot, width = 6, height = 4)


conj_rate_crittimes_plot <- ggplot(molten_test, aes(x = gamma.T,
                                                    y = value)) +
  geom_point(aes(x = gamma.T,
                 y = crit_time)) +
  labs(x = bquote('Transconjugant conjugation rate'~gamma[Tmax]),
       y = 'Critical time (h)') +
  #scale_y_continuous(trans = 'log', breaks = .base_breaks()) +
  scale_y_continuous(limits = c(6, 13)) +
  scale_x_continuous(trans = 'log', breaks = .base_breaks()) +
  theme_minimal() + theme(text = element_text(size=20))

conj_rate_crittimes_plot
ggsave(filename = paste0('../../figures/',
                         #'conj_rate_crit_times_high_gamma.png'),
                         'conj_rate_crit_times.png'),
       plot = conj_rate_crittimes_plot, width = 5, height = 4)

## Combine plots ----
library('cowplot')

all_plots <- plot_grid( growth_rate_fc_plot + theme(legend.position = "none"),
                        conj_rate_fc_plot + theme(legend.position = c(0.2, 0.6)),
                        growth_rate_crittimes_plot,
                        conj_rate_crittimes_plot,
                        align = 'v', axis = 'l', nrow = 2, ncol = 2,
                        labels = "AUTO", label_size = 20)
all_plots
ggsave(filename = paste0('../../figures/',
                         'conj_growth_rate_combined.png'),
       plot = all_plots, width = 11, height = 10)

#####
# reduce the number of x-axis options and make boxplot
full_out <- full_out[full_out$gamma.T %in% 10**seq(-15, -10, 1), ]
boxplot_scan(full_out, 'gamma.T', 'TD', 't',
             plotfile = '../../figures/TD_afo_gammat_time_log.png',
             logplot = TRUE, add_means = TRUE)
