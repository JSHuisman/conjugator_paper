###########################################################
# Figure 1
# Author: Jana Huisman
###########################################################

library(conjugator)
library(tidyverse)

library(deSolve)
library(ggplot2)
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
# Plots for in the paper
vars <- c(r=1E6, t=0, d=5E6, c = 1E12) # initial population size
# growth rates, conjugation rates
t_meas = 15
t.vec <- seq(0, t_meas, length = 100)
tol_factor <- 10

parms_cost <- c(psi.R=1.0, psi.T=0.7, psi.D=0.7,
                gamma.t=10^-11, gamma.d=10^-14, q = 1E2)

titles = c('cost')
parms = rbind(parms_cost)
i = 1

  out <- integrate_model(vars, t.vec, "model_ESM", parms[i,])
  crit_times <- estimate_crittime_from_sim(parms[i,], vars, tol_factor)
  
  simple_dyn_plot <- plot_dynamics_ggplot(out, crit_times) +
    scale_color_manual(breaks = c('r', 'd', 't'),
                       labels = c('Recipients', 'Donors', 'Transconjugants'),
                       values = col_palet[c(3,6,9)])
  ggsave(filename = paste0('../figures/',
                           'dyn_plot_',titles[i],'.png'),
         plot = simple_dyn_plot, width = 6, height = 4)
  
  conj_result <- estimate_conj_from_sim(out, method = c('TD', 'T_DR', 'T_RT', 
                                                        'Gama', 'ASM', 'SM'), 
                                        parms[i,], id_cols = 't')
  
  # take absolute values
  conj_result <- conj_result %>%
    mutate(t = as.numeric(levels(t))[t],
           estimate = ifelse(method %in% c('Gama'), abs(estimate), estimate))
  
  # ALL CONJ RATE MEASURES
  # needs general ggplot with NULL; otherwise tcrit can't be added
  conj_rate_plot <- ggplot(NULL) +
    geom_line(data = conj_result, aes(x = t, y = estimate, color = method)) +
    scale_y_continuous(trans = 'log', breaks = .base_breaks()) +
    coord_cartesian(xlim = c(0, t_meas)) +
    scale_color_manual(name = 'Estimation method',
                       breaks = c('TD', 'T_RT', 'Gama', 'T_DR','ASM', 'SM'),
                       #values = c('red', 'orange', 'pink', 'purple','blue', 'lightblue'),
                       values = col_palet[c(1,4,7,5,8,2)],
                       labels = c('T/D', 'T/(T+R)', bquote(abs(log('T'/sqrt('DR')))),
                                  bquote('T/DR'),
                                  bquote(gamma[Dmax]~'in mL/(CFUxHour)'),
                                  bquote(gamma[max]~'in mL/(CFUxHour)'))) +
    labs(x = "Time (h)", y = "Conjugation estimate") +
    theme_light() + 
    theme(text = element_text(size=15))
  
  ggsave(filename = paste0('../../conjugationrateestimation/figures/',
                           'all_conj_measures_',titles[i],'.png'),
         plot = conj_rate_plot, width = 6, height = 4)
  
  conj_rate_plot <- add_tcrit_to_plot(conj_rate_plot, crit_times,
                                      y_pos = min(parms[i, "psi.R"],
                                                  parms[i, "psi.T"],
                                                  parms[i, "psi.D"]), id_cols = 'ID')
  
  
  dyn_legend <- get_legend(simple_dyn_plot + theme(legend.position = "right"))
  conj_rate_legend <- get_legend(conj_rate_plot + theme(legend.position = "right"))
  
  all_plots <- plot_grid( simple_dyn_plot + theme(legend.position = "none"),
                          dyn_legend,
                          conj_rate_plot + theme(legend.position = "none"),
                          conj_rate_legend,
                          nrow = 2, ncol = 2,
                          labels = c("A", "", "B", ""),
                          rel_widths = c(1, .4, 1, .4),
                          label_size = 20)
  
  all_plots
  ggsave(filename = paste0('../figures/',
                           'dyn_plot_conj_cost.png'),
         plot = all_plots, width = 8, height = 7)
