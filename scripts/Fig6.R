###########################################################
# Figure 5
# Author: Jana Huisman
###########################################################

library(conjugator)
library(tidyverse)
 
library(ggplot2)
library(patchwork)

source('./R/plotting_functions.R')

###########################################################
col_palet <- c('#332288','#88CCEE',
               '#44AA99','#117733',
               '#999933','#DDCC77',
               '#CC6677','#882255',
               '#AA4499','#DDDDDD',
               '#000000')
###########################################################

## Estimate conjugation rates ####
data_full <- read_csv(paste0('../data/', 'Full_protocol_DRT_TRT.csv'))
result <- estimate_conj_rate(data_full, c("SM", "ASM"), id_cols = c('Transfer', 'Replicate', 't'))

result <- result %>%
  mutate(estimate = ifelse(estimate <= 0, NA, estimate),
         #log_estimate = log10(estimate),
         Replicate = factor(Replicate))

ggplot(result, aes(x = Transfer, y = estimate, colour = method)) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  geom_point(position = position_dodge(width = 0.75), aes(group = method)) +
  scale_y_continuous(trans = 'log10', breaks = .base_breaks(), limits = c(1e-12, 1e-9)) +
  facet_wrap(vars(t), labeller = labeller(t = c('4' = '4 hours', '24' = '24 hours'))) +
  labs(y = 'Conjugation Rate Estimate', colour = 'Method') +
  scale_color_manual(values = col_palet[c(3,9)]) +
  theme_minimal() + theme(text = element_text(size=15))

ggsave(paste0('../figures/', 'Full_protocol_conj_estimate.pdf'), width = 10, height = 4)

## Crit Time Estimation ####
DRT_data <- data_full %>%
  filter(Transfer == 'DRT')

TRT_data <- data_full %>%
  filter(Transfer == 'DRT')

estimate_crit_time(DRT_data, TRT_data, id_cols = c('Replicate', 't')) %>%
  group_by(t) %>%
  summarise(min_tcrit = min(min_tcrit))

########
#New fig with points for revisions
left_plot = ggplot(result %>% filter(t == '4'), aes(x = Transfer, y = estimate, colour = method)) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  geom_point(position = position_dodge(width = 0.75), aes(group = method)) +
  scale_y_continuous(trans = 'log10', breaks = .base_breaks(), limits = c(1e-12, 1e-9)) +
  facet_wrap(vars(t), labeller = labeller(t = c('4' = '4 hours', '24' = '24 hours'))) +
  labs(y = 'Conjugation Rate Estimate', colour = 'Method') +
  scale_color_manual(values = col_palet[c(3,9)]) +
  theme_minimal() + 
  theme(text = element_text(size=15),
        axis.title.x = element_blank())

right_plot = ggplot(result %>% filter(t == '24'), aes(x = Transfer, y = estimate, colour = method)) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  geom_point() +
  scale_y_continuous(trans = 'log10', breaks = .base_breaks(), limits = c(1e-12, 1e-9)) +
  facet_wrap(vars(t), labeller = labeller(t = c('4' = '4 hours', '24' = '24 hours'))) +
  labs(y = 'Conjugation Rate Estimate', colour = 'Method') +
  scale_color_manual(values = col_palet[c(3,9)]) +
  theme_minimal() + 
  theme(text = element_text(size=15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

left_plot + right_plot +
  plot_layout(ncol = 2, guides = 'collect') + 
  theme(legend.position = 'right')

ggsave(paste0('../figures/', 'Full_protocol_conj_estimate.pdf'), width = 10, height = 4)
