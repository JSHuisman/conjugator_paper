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
library(cowplot)

source('../R/model_functions.R')
source('../R/plotting_functions.R')
source('../R/utils.R')

col_palet <- c('#332288','#88CCEE',
               '#44AA99','#117733',
               '#999933','#DDCC77',
               '#CC6677','#882255',
               '#AA4499','#DDDDDD',
               '#000000')

###########################################################
# varying conjugation rate ----

sim_as_function_of_conj <- function(vars, t.vec, psi, gammaD_vec, gammaT_vec){
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
      gamma_max_conjugation[gd_ind, gt_ind] <- conj_result_SM$estimate[100]
      gammaD_conjugation[gd_ind, gt_ind] <- conj_result_ASM$estimate[100]
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

plot_gamma_max_afo_conj <- function(gamma_max_est, return_plot = TRUE){
  gamma_max_conj_molten <- melt(gamma_max_est)
  
  gamma_max_conj_molten[gamma_max_conj_molten == Inf] <- NA
  
  gamma_max_conj_molten['error'] <- 100*(gamma_max_conj_molten$value - gamma_max_conj_molten$gammaD)/gamma_max_conj_molten$gammaD
  # p_perc <- ggplot(gamma_max_conj_molten, aes(x=gammaD, y = gammaT, fill = log(error))) +
  #   geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
  #   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
  #   ggplot2::labs(x = bquote("Conjugation rate"~gamma[D]),
  #                 y = bquote("Ratio of conjugation rates"~gamma[T]/gamma[D]),
  #                 fill = bquote(log("error")~gamma[max])) +
  #   ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  # ggsave(filename = paste0('../../conjugationrateestimation/figures/',
  #                          'gamma_max_afo_conj_perc','.png'),
  #        plot = p_perc, width = 6, height = 4)
  
  gamma_max_conj_molten['fold_change'] <- gamma_max_conj_molten$value/gamma_max_conj_molten$gammaD
  p_foldchange<- ggplot(gamma_max_conj_molten, aes(x=gammaD, y = gammaT, fill = log(fold_change))) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    ggplot2::labs(x = bquote("Conjugation rate"~gamma[Dmax]),
                  y = bquote("Ratio"~gamma[Tmax]/gamma[Dmax]),
                  fill = bquote(atop("Fold change", gamma[max]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gamma_max_afo_conj_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}

plot_gammaD_afo_conj <- function(gammaD_est, return_plot = TRUE){
  gammaD_conj_molten <- melt(gammaD_est)
  
  gammaD_conj_molten[gammaD_conj_molten < 0] <- NA
  
  gammaD_conj_molten['error'] <- 100*(gammaD_conj_molten$value - gammaD_conj_molten$gammaD)/gammaD_conj_molten$gammaD
  # p_perc <- ggplot(gammaD_conj_molten, aes(x=gammaD, y = gammaT, fill = log(error))) +
  #   geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
  #   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
  #   ggplot2::labs(x = bquote("Conjugation rate"~gamma[D]),
  #                 y = bquote("Ratio of conjugation rates"~gamma[T]/gamma[D]),
  #                 fill = bquote(log("error")~gamma[D])) +
  #   ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  # ggsave(filename = paste0('../../conjugationrateestimation/figures/',
  #                          'gammaD_afo_conj_perc','.png'),
  #        plot = p_perc, width = 6, height = 4)
  
  gammaD_conj_molten['fold_change'] <- gammaD_conj_molten$value/gammaD_conj_molten$gammaD
  p_foldchange <- ggplot(gammaD_conj_molten, aes(x=gammaD, y = gammaT, fill = log(fold_change))) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    ggplot2::labs(x = bquote("Conjugation rate"~gamma[Dmax]),
                  y = bquote("Ratio"~gamma[Tmax]/gamma[Dmax]),
                  fill = bquote(atop("Fold change", gamma[Dmax]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gammaD_afo_conj_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}

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

gmax_conj <- plot_gamma_max_afo_conj(conj_scan$gamma_max_est)
gD_conj <- plot_gammaD_afo_conj(conj_scan$gammaD_est)

# for t = 8
# It seems there is some kind of maximal error
# at high gammaT, and intermediate gammaD
# presumably because for low gammaD there are not enough transconjugants
# generated to be a problem, whereas at high gammaD, this value
# can be estimated with less error/problem???

# for t_meas = 12
# Simonsen method returns very bad results

##########################################################################
# varying growth rate ----

sim_as_function_of_growth <- function(vars, t.vec, gamma, psiD_options, psiR_options){
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
      gamma_max_growth[psiD_ind, psiR_ind] <- conj_result_SM$estimate[100]
      gammaD_growth[psiD_ind, psiR_ind] <- conj_result_ASM$estimate[100]
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

plot_gamma_max_afo_growth <- function(gamma_max_est, return_plot = TRUE){
  gamma_max_growth_molten <- melt(gamma_max_est)
  
  gamma_max_growth_molten['perc_error'] <- 100*(gamma_max_growth_molten$value - gamma)/gamma
  p_perc <- ggplot(gamma_max_growth_molten, aes(x= psiD, y = psiR, fill = perc_error)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote("Percent error"~gamma[max])) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gamma_max_afo_growth_perc','.png'),
         plot = p_perc, width = 6, height = 4)
  
  gamma_max_growth_molten['fold_change'] <- gamma_max_growth_molten$value/gamma
  p_foldchange<- ggplot(gamma_max_growth_molten, aes(x= psiD, y = psiR, fill = fold_change)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red", midpoint = 1) +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote(atop("Fold change", ~gamma[max]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gamma_max_afo_growth_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
  
}

plot_gammaD_afo_growth <- function(gammaD_est, return_plot = TRUE){
  gammaD_growth_molten <- melt(gammaD_est)
  
  gammaD_growth_molten['error'] <- 100*(gammaD_growth_molten$value - gamma)/gamma
  # p_perc <- ggplot(gammaD_growth_molten, aes(x = psiD, y = psiR, fill = error)) +
  #   geom_tile() + scale_fill_gradient2(low = "blue", high = "red") +
  #   ggplot2::labs(x = bquote("Growth rate"~psi[D]), y = bquote("Growth rate"~psi[R]),
  #                 fill = bquote("Percent error"~gamma[D])) +
  #   ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  # ggsave(filename = paste0('../../conjugationrateestimation/figures/',
  #                          'gammaD_afo_growth_foldchange','.png'),
  #        plot = p_perc, width = 6, height = 4)
  
  gammaD_growth_molten['fold_change'] <- gammaD_growth_molten$value/gamma
  p_foldchange <- ggplot(gammaD_growth_molten, aes(x = psiD, y = psiR, fill = fold_change)) +
    geom_tile() + scale_fill_gradient2(low = "blue", high = "red", midpoint = 1) +
    ggplot2::labs(x = bquote("Growth rate"~psi[Dmax]), y = bquote("Growth rate"~psi[Rmax]),
                  fill = bquote(atop("Fold change",gamma[Dmax]))) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15))
  ggsave(filename = paste0('../figures/',
                           'gammaD_afo_growth_foldchange','.png'),
         plot = p_foldchange, width = 6, height = 4)
  
  if (return_plot){
    return(p_foldchange)
  }
}

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

################
# Combine plots for paper ----

tmin_growth_legend <- get_legend(tmin_growth +
                                   theme(legend.position = "right", text = element_text(size=15)))
tmin_conj_legend <- get_legend(tmin_conj +
                                 theme(legend.position = "right", text = element_text(size=15)))
gmax_growth_legend <- get_legend(gmax_growth +
                                   theme(legend.position = "right", text = element_text(size=15)))
gmax_conj_legend <- get_legend(gmax_conj +
                                 theme(legend.position = "right", text = element_text(size=15)))
gD_growth_legend <- get_legend(gD_growth +
                                 theme(legend.position = "right", text = element_text(size=15)))
gD_conj_legend <- get_legend(gD_conj +
                               theme(legend.position = "right", text = element_text(size=15)))



all_plots <- plot_grid( tmin_growth + theme(legend.position = "none"),
                        tmin_growth_legend,
                        tmin_conj + theme(legend.position = "none"),
                        tmin_conj_legend,
                        gmax_growth + theme(legend.position = "none"),
                        gmax_growth_legend,
                        gmax_conj + theme(legend.position = "none"),
                        gmax_conj_legend,
                        gD_growth + theme(legend.position = "none"),
                        gD_growth_legend,
                        gD_conj + theme(legend.position = "none"),
                        gD_conj_legend,
                        nrow = 3, ncol = 4,
                        labels = c("A", "", "B", "", "C", "", "D",
                                   "", "E", "", "F", ""),
                        rel_widths = c(1, .3, 1, .3,
                                       1, .3, 1, .3,
                                       1, .3, 1, .3),
                        label_size = 15)

all_plots

ggsave(filename = paste0('../figures/',
                         'conj_growth_rate_scans_combined.pdf'),
       plot = all_plots, width = 11, height = 9)




