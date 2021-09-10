###########################################################
# Experimentally observed differences
# Author: Jana Huisman
###########################################################

library(conjugator)

library(deSolve)

library(tidyverse) # config df
library(ggplot2)
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
## Different treatments; measure at 4, 8, 24 hours

nrep = 5
ntreat = 5
base_growth = 1.2
growth_sd = base_growth*5e-2

psiR_treat = c(rnorm(nrep, base_growth, growth_sd), 
               rnorm(nrep, base_growth, growth_sd),
               rnorm(nrep, 1.5*base_growth, 1.5*growth_sd),
               rnorm(nrep, 0.5*base_growth, 0.5*growth_sd),
               rnorm(nrep, 0.5*base_growth, 0.5*growth_sd))
psiD_treat = c(rnorm(nrep, base_growth, growth_sd), 
               rnorm(nrep, 0.75*base_growth, 0.75*growth_sd),
               rnorm(nrep, base_growth, growth_sd),
               rnorm(nrep, 0.5*base_growth, 0.5*growth_sd),
               rnorm(nrep, 0.5*base_growth, 0.5*growth_sd))
psiT_treat = c(rnorm(nrep, base_growth, growth_sd), 
               rnorm(nrep, 0.75*base_growth, 0.75*growth_sd),
               rnorm(nrep, 1.5*base_growth, 1.5*growth_sd),
               rnorm(nrep, base_growth, growth_sd),
               rnorm(nrep, 0.5*base_growth, 0.5*growth_sd))

fixed_vars <- cbind(r = rnorm(nrep*ntreat, mean = 1e6, sd = 5e4), 
              t = rep(0, nrep*ntreat),
              d = rnorm(nrep*ntreat, mean = 1e6, sd = 5e4), 
              c = rep(1E14, nrep*ntreat),
              q = rep(1E2, nrep*ntreat),
              gamma.t = rnorm(nrep*ntreat, mean = 1e-13, sd = 5e-15), 
              gamma.d = rnorm(nrep*ntreat, mean = 1e-13, sd = 5e-15),
              psi.R = psiR_treat,
              psi.D = psiD_treat, 
              psi.T = psiT_treat,
              treatment = unlist(lapply(1:ntreat, function(x){rep(x, nrep)})) )
config_df <- as_tibble(fixed_vars)

## compute outputs ####
t.vec <- seq(0, 24, length = 241)
model_out <- compute_model_scan(config_df, t.vec, output_timepoints = c(41, 81, 241), 
                                psi_time = 20) 
# a psi measurement in stationary phase (psi_time = 240) 
#reduces gamma estimates for first 3 treatments

full_out <- .rename_fullout(model_out$full_out)

results <- conjugator::estimate_conj_rate(full_out, c("T_RT", "TD", "Dionisio", "T_DR", "ASM", "SM"),
                                          id_cols = c('t', 'treatment')) %>%
  select(t, treatment, method, estimate) %>%
  mutate(t = factor(t),
         method = factor(method, levels = c("TD", "T_RT", "Dionisio", "T_DR", "SM", "ASM"),
                         labels = c("T/D", "T/(R+T)", "log(T/sqrt(DR))", "T/DR", 
                                    "SM endpoint", "ASM endpoint")),
         treatment = factor(treatment))

conjugator::scan_crit_time(full_out %>% filter(treatment == 1, t == 4), 
                           mult_seq = 1, id_cols = c('t', 'treatment'))

## Plotting the results ####

## We need to split the plot to fix the y-axis
top_row <- ggplot(results %>% filter(method %in% c("T/D", "T/(R+T)", "log(T/sqrt(DR))"))) +
  geom_boxplot(aes(y = estimate, x = t, 
                   colour = treatment, group = interaction(t, treatment))) +
  facet_wrap(vars(method)) +
  scale_y_continuous(trans = 'log', breaks = .base_breaks(), limits = c(1e-6, 1e2)) +
  labs(x = 'Time', y = 'Conjugation proficiency',
       colour = 'Treatment') +
  scale_color_manual(name = 'Treatment',
                     values = col_palet[c(1,4,7,5,8,2)],
                     breaks = 1:5,
                     labels = c('Equal growth', 'Plasmid cost', 'Fast recipient',
                                'Slow Donor/Recipient', 'All slow')) 

## We need to split the plot to fix the y-axis
bottom_row <- ggplot(results %>% filter(method %in% c("T/DR", 
                                                      "SM endpoint", "ASM endpoint"))) +
  geom_boxplot(aes(y = estimate, x = t, 
                   colour = treatment, group = interaction(t, treatment))) +
  facet_wrap(vars(method)) +
  scale_y_continuous(trans = 'log', breaks = .base_breaks(), limits = c(5e-14, 1e-10)) +
  labs(x = 'Time', y = 'Conjugation rate',
       colour = 'Treatment') +
  scale_color_manual(name = 'Treatment',
                     values = col_palet[c(1,4,7,5,8,2)],
                     breaks = 1:5,
                     labels = c('Equal growth', 'Plasmid cost', 'Fast recipient',
                                'Slow Donor/Recipient', 'All slow')) 

top_row + bottom_row +
  plot_layout(ncol = 1, guide = 'collect') &
  theme_minimal() & 
  theme(text = element_text(size=15),
        legend.position = 'bottom')


ggsave(filename = paste0('../figures/',
                         'Growth_experiments.png'),
       width = 12, height = 8)

