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

ggsave(filename = paste0('../figures/',
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
ggsave(filename = paste0('../figures/',
                         'Growth_rate_crit_times.png'),
       plot = growth_rate_crittimes_plot, width = 3, height = 4)

# higher initial densities, more fold-change impact of
# low Donor growth; strong time-dependence of measurement

##########################################################################
# varying growth rate functions ----

sim_as_function_of_growth <- function(vars, t.vec, gamma, psiD_options, psiR_options,
                                      psi_window = 30){
  tmin_growth = array(data=NA,
                      dim = c(length(psiD_options), length(psiR_options)),
                      dimnames = list(psiD = psiD_options,
                                      psiR = psiR_options))
  gamma_max_growth = array(data=NA,
                           dim = c(length(psiD_options), length(psiR_options)),
                           dimnames = list(psiD = psiD_options,
                                           psiR = psiR_options))
  gammaD_growth = array(data=NA,
                        dim = c(length(psiD_options), length(psiR_options)),
                        dimnames = list(psiD = psiD_options,
                                        psiR = psiR_options))
  
  for (psiD_ind in 1:length(psiD_options)){
    for (psiR_ind in 1:length(psiR_options)){
      
      parms <- c(psi.R = psiR_options[psiR_ind],
                 psi.T = 1.0,
                 psi.D = psiD_options[psiD_ind],
                 gamma.t = gamma, gamma.d = gamma, q = 1E2)
      
      # calculate critical times
      tmin_growth[psiD_ind, psiR_ind] = estimate_crittime_from_sim(parms, vars, tol_factor = 10)[['min_tcrit']]
      # simulate model with these values
      out <- integrate_model(vars, t.vec, "model_ESM", parms)
      
      conj_result_SM <- estimate_conj_from_sim(out, c('SM'), parms, exp_window = 100)
      conj_result_ASM <- estimate_conj_from_sim(out, c('ASM'), parms, exp_window = 100)
      gamma_max_growth[psiD_ind, psiR_ind] <- conj_result_SM$estimate[psi_window]
      gammaD_growth[psiD_ind, psiR_ind] <- conj_result_ASM$estimate[psi_window]
    }
  }
  return(list(tmin_est = tmin_growth,
              gamma_max_est = gamma_max_growth,
              gammaD_est = gammaD_growth))
}

plot_tmin_afo_growth <- function(tmin_est, t_meas, return_plot = TRUE){
  tmin_growth_molten <- melt(tmin_est)
  
  min(tmin_growth_molten$value)
  
  p_tmin <- ggplot(tmin_growth_molten, aes(x = psiD, y = psiR, fill = value)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red", midpoint = t_meas) +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote(atop("Minimal", t[crit]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'tmin_afo_growth','.png'),
         plot = p_tmin, width = 6, height = 4)
  
  if (return_plot){
    return(p_tmin)
  }
}

plot_gamma_max_afo_growth <- function(gamma_max_est, return_plot = TRUE, colour_lims = c(0, 2)){
  gamma_max_growth_molten <- melt(gamma_max_est)
  
  gamma_max_growth_molten['perc_error'] <- 100*(gamma_max_growth_molten$value - gamma)/gamma
  
  gamma_max_growth_molten['fold_change'] <- gamma_max_growth_molten$value/gamma
  p_foldchange<- ggplot(gamma_max_growth_molten, aes(x= psiD, y = psiR, fill = fold_change)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 1, limits = colour_lims) +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote(atop("Fold change", ~gamma[max]))) +
    ggplot2::theme_minimal() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gamma_max_afo_growth_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
  
}

plot_gammaD_afo_growth <- function(gammaD_est, return_plot = TRUE, colour_lims = c(0,2)){
  gammaD_growth_molten <- melt(gammaD_est)
  
  gammaD_growth_molten['error'] <- 100*(gammaD_growth_molten$value - gamma)/gamma
  
  gammaD_growth_molten['fold_change'] <- gammaD_growth_molten$value/gamma
  p_foldchange <- ggplot(gammaD_growth_molten, aes(x = psiD, y = psiR, fill = fold_change)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 1, limits = colour_lims) +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote(atop("Fold change",gamma[D]))) +
    ggplot2::theme_minimal() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gammaD_afo_growth_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}


##########################################################################
# varying growth rate ----
###
vars <- c(r=5E6, t=0, d=5E6, c = 1E14) # initial population size
gamma <- 10**-13
t_meas = 8
t.vec <- seq(0, t_meas, length = 100)
#
psiD_options = seq(0.5, 1.5, 0.1)
psiR_options = seq(0.5, 1.5, 0.1) # assumption that T is like R
#
growth_scan <- sim_as_function_of_growth(vars, t.vec, gamma, psiD_options, psiR_options)
tmin_growth <- plot_tmin_afo_growth(growth_scan$tmin_est, t_meas)
# The minimal critical time does not depend very strongly
# on the relative growth rate, just the overal growth

gmax_growth <- plot_gamma_max_afo_growth(growth_scan$gamma_max_est)
gD_growth <- plot_gammaD_afo_growth(growth_scan$gammaD_est)
# For high growth rates of D/R with respect to T ( psiT set to 1)
# the Simonsen method underestimates the conjugation rate
# for low growth rates the conjugation rate is overestimated
# this is symmetric in D and R

# but it is most likely that D is a different species than T/R
# and R would grow slower because of selective pressure for plasmid?
# or virulence factor on plasmid?!


## Combine plots ----

growth_rate_fc_plot + growth_rate_crittimes_plot +
  gmax_growth + tmin_growth + 
 gD_growth + plot_spacer() +
  plot_layout(ncol = 2, guides = 'keep') + 
  plot_annotation(tag_levels = 'A') &
  theme_minimal() &
  theme(legend.position = 'right',
        text = element_text(size=15)) 

ggsave(filename = paste0('../figures/',
                         'growth_rate_combined.png'),
       width = 11, height = 10)


