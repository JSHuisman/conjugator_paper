###########################################################
# Figure 1
# Author: Jana Huisman
###########################################################

library(conjugator)
library(tidyverse)

library(deSolve)
library(ggplot2)
library(cowplot)
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
                     values = col_palet[c(3,6,9)]) +
  geom_hline(yintercept = 1, linetype = 'dashed')
simple_dyn_plot

# ggsave(filename = paste0('../figures/',
#                          'dyn_plot_',titles[i],'.png'),
#        plot = simple_dyn_plot, width = 6, height = 4)

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
                                bquote(gamma[D]~'in mL/(CFUxHour)'),
                                bquote(gamma[max]~'in mL/(CFUxHour)'))) +
  labs(x = "Time (h)", y = "Conjugation rate") +
  theme_minimal() + 
  theme(text = element_text(size=15))

conj_rate_plot

# ggsave(filename = paste0('../../conjugationrateestimation/figures/',
#                          'all_conj_measures_',titles[i],'.png'),
#        plot = conj_rate_plot, width = 6, height = 4)

conj_rate_plot <- add_tcrit_to_plot(conj_rate_plot, crit_times,
                                    y_pos = min(parms[i, "psi.R"],
                                                parms[i, "psi.T"],
                                                parms[i, "psi.D"]), id_cols = 'ID')
###############################
# Experimental data
data_timecourse <- read_csv('../data/Time_course_DRT.csv')

growth_cols_S = list(psi.D = 1.46, psi.R = 1.43, psi.T = 1.40)
data_full <- cbind(data_timecourse, growth_cols_S)
result <- estimate_conj_rate(data_full, c('SM', 'TD', 'T_DR', 'T_RT', 'Gama'), id_cols = c('Replicate', 't'))

result <- result %>%
  mutate(Replicate = factor(Replicate),
         t = as.numeric(levels(t))[t]) %>%
  mutate(estimate = ifelse(method == 'Gama', abs(estimate), estimate))

exp_plot <- ggplot(result, aes(x = t, y = estimate, fill = method)) +
  geom_boxplot(aes(group = interaction(t, method)), position = "identity", 
               width = 0.75, show.legend = F) +
  geom_line(aes(colour = method), show.legend = F) +
  scale_y_continuous(trans = 'log10', breaks = .base_breaks()) +
  scale_fill_manual(name = 'Estimation method',
                      breaks = c('TD', 'T_RT', 'Gama', 'T_DR','SM'),
                      values = col_palet[c(1,4,7,5,2)],
                      labels = c('T/D', 'T/(T+R)', bquote(abs(log('T'/sqrt('DR')))),
                                 bquote('T/DR'),
                                 bquote(gamma[D]~'in mL/(CFUxHour)'),
                                 bquote(gamma[max]~'in mL/(CFUxHour)'))) +
  scale_colour_manual(name = 'Estimation method',
                     breaks = c('TD', 'T_RT', 'Gama', 'T_DR','SM'),
                     values = col_palet[c(1,4,7,5,2)],
                     labels = c('T/D', 'T/(T+R)', bquote(abs(log('T'/sqrt('DR')))),
                                bquote('T/DR'),
                                bquote(gamma[D]~'in mL/(CFUxHour)'),
                                bquote(gamma[max]~'in mL/(CFUxHour)'))) +
  labs(y = 'Conjugation rate', x = 'Time (h)') +
  theme_minimal() + 
  theme(text = element_text(size=15))

exp_plot

#########
simple_dyn_plot + guide_area() + #plot_spacer() +
  conj_rate_plot + exp_plot +
  plot_layout(ncol = 2, guides = 'collect') + 
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'right')

ggsave(filename = paste0('../figures/dyn_plot_conj_cost.png'), 
       width = 10, height = 8)
