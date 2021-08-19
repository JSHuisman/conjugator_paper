###########################################################
# Figure 4
# Author: Jana Huisman
###########################################################

library(conjugator)
library(tidyverse)

library(deSolve)
library(ggplot2)
# library(reshape2)
# library(ggrepel)
library(gridExtra) # to arrange multiple ggplots
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

ggsave(filename = paste0('../figures/',
                         'conj_rate_fc.png'),
       plot = conj_rate_fc_plot, width = 6, height = 4)


conj_rate_crittimes_plot <- ggplot(molten_test, aes(x = gamma.T,
                                                    y = value)) +
  geom_point(aes(x = gamma.T,
                 y = crit_time)) +
  labs(x = bquote('Transconjugant conjugation rate'~gamma[Tmax]),
       y = 'Critical time (h)') +
  scale_y_continuous(limits = c(6, 13)) +
  scale_x_continuous(trans = 'log', breaks = .base_breaks()) +
  theme_minimal() + theme(text = element_text(size=20))

conj_rate_crittimes_plot
ggsave(filename = paste0('../figures/',
                         'conj_rate_crit_times.png'),
       plot = conj_rate_crittimes_plot, width = 5, height = 4)

###########################################################
# varying conjugation rate functions ----

sim_as_function_of_conj <- function(vars, t.vec, psi, gammaD_vec, gammaT_vec,
                                    psi_window = 30){
  tmin_conjugation = array(data=NA, dim = c(length(gammaD_vec),
                                            length(gammaT_vec)), dimnames = list(
                                              gammaD = gammaD_vec, gammaT = gammaT_vec))
  gamma_max_conjugation = array(data=NA, dim = c(length(gammaD_vec),
                                                 length(gammaT_vec)), dimnames = list(
                                                   gammaD = gammaD_vec, gammaT = gammaT_vec))
  gammaD_conjugation = array(data=NA, dim = c(length(gammaD_vec),
                                              length(gammaT_vec)), dimnames = list(
                                                gammaD = gammaD_vec, gammaT = gammaT_vec))
  
  for (gd_ind in 1:length(gammaD_vec)){
    for (gt_ind in 1:length(gammaT_vec)){
      gammaD = gammaD_vec[gd_ind]
      gammaT = gammaT_vec[gt_ind]*gammaD
      
      parms <- c(psi.R = psi, psi.T = psi, psi.D = psi,
                 gamma.t = gammaT, gamma.d = gammaD,
                 q = 1E2)
      # calculate critical times
      tmin_conjugation[gd_ind, gt_ind] = estimate_crittime_from_sim(parms, vars, tol_factor = 10)[['min_tcrit']]
      # simulate model with these values
      out <- integrate_model(vars, t.vec, "model_ESM", parms)
      conj_result_SM <- estimate_conj_from_sim(out, method = 'SM', parms, exp_window = 20)
      conj_result_ASM <- estimate_conj_from_sim(out, method = 'ASM', parms, exp_window = 20)
      gamma_max_conjugation[gd_ind, gt_ind] <- conj_result_SM$estimate[psi_window]
      gammaD_conjugation[gd_ind, gt_ind] <- conj_result_ASM$estimate[psi_window]
    }
  }
  
  return(list(tmin_est = tmin_conjugation, gamma_max_est = gamma_max_conjugation, gammaD_est = gammaD_conjugation))
}

plot_tmin_afo_conj <- function(tmin_est, t_meas, return_plot = TRUE){
  tmin_conj_molten <- melt(tmin_est)
  min(tmin_conj_molten$value)
  
  p_tmin<- ggplot(tmin_conj_molten, aes(x=gammaD, y = gammaT, fill = value)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red", midpoint = t_meas) +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    ggplot2::labs(x = bquote("Conjugation rate"~gamma[Dmax]),
                  y = bquote("Ratio"~gamma[Tmax]/gamma[Dmax]),
                  fill = bquote(atop("Minimal", t[crit]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'tmin_afo_conj','.png'),
         plot = p_tmin, width = 6, height = 4)
  
  if (return_plot){
    return(p_tmin)
  }
}

plot_gamma_max_afo_conj <- function(gamma_max_est, return_plot = TRUE, colour_lims = c(0,8)){
  gamma_max_conj_molten <- melt(gamma_max_est)
  
  gamma_max_conj_molten[gamma_max_conj_molten == Inf] <- NA
  
  gamma_max_conj_molten['error'] <- 100*(gamma_max_conj_molten$value - gamma_max_conj_molten$gammaD)/gamma_max_conj_molten$gammaD
  
  gamma_max_conj_molten['fold_change'] <- gamma_max_conj_molten$value/gamma_max_conj_molten$gammaD
  p_foldchange<- ggplot(gamma_max_conj_molten, aes(x=gammaD, y = gammaT, fill = log(fold_change))) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = colour_lims) +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    ggplot2::labs(x = bquote("Conjugation rate"~gamma[Dmax]),
                  y = bquote("Ratio"~gamma[Tmax]/gamma[Dmax]),
                  fill = bquote(atop("Fold change", gamma[max]))) +
    ggplot2::theme_minimal() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gamma_max_afo_conj_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}

plot_gammaD_afo_conj <- function(gammaD_est, return_plot = TRUE, colour_lims = c(0,8)){
  gammaD_conj_molten <- melt(gammaD_est)
  
  gammaD_conj_molten[gammaD_conj_molten < 0] <- NA
  
  gammaD_conj_molten['error'] <- 100*(gammaD_conj_molten$value - gammaD_conj_molten$gammaD)/gammaD_conj_molten$gammaD
  
  gammaD_conj_molten['fold_change'] <- gammaD_conj_molten$value/gammaD_conj_molten$gammaD
  p_foldchange <- ggplot(gammaD_conj_molten, aes(x=gammaD, y = gammaT, fill = log(fold_change))) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = colour_lims) +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    ggplot2::labs(x = bquote("Conjugation rate"~gamma[Dmax]),
                  y = bquote("Ratio"~gamma[Tmax]/gamma[Dmax]),
                  fill = bquote(atop("Fold change", gamma[Dmax]))) +
    ggplot2::theme_minimal() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gammaD_afo_conj_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}

##########################################################################
# varying conj rate ----

#####
vars <- c(r=5E6, t=0, d=5E6, c = 1E14) # initial population size
psi <- 1.0
t_meas = 8
t.vec <- seq(0, t_meas, length = 100)
#
gammaD_vec = 10**seq(-15, -10, 0.5)
gammaT_vec = 10**seq(-2, 4, 1)

conj_scan <- sim_as_function_of_conj(vars, t.vec, psi, gammaD_vec, gammaT_vec)
tmin_conj <- plot_tmin_afo_conj(conj_scan$tmin_est, t_meas)
# The minimal critical time depends very strongly on the
# relative conjugation rates, as well as absolute magnitude

gmax_conj <- plot_gamma_max_afo_conj(conj_scan$gamma_max_est, colour_lims = c(-0.1, 10))
gD_conj <- plot_gammaD_afo_conj(conj_scan$gammaD_est, colour_lims = c(-0.1, 10))

################
## Combine plots ----

conj_rate_fc_plot + conj_rate_crittimes_plot +
  gmax_conj + tmin_conj + 
  gD_conj + plot_spacer() +
  plot_layout(ncol = 2, guides = 'keep') + 
  plot_annotation(tag_levels = 'A') &
  theme_minimal() &
  theme(legend.position = 'right',
        text = element_text(size=15)) 

ggsave(filename = paste0('../figures/',
                         'conj_rate_combined.png'),
       width = 11, height = 10)



