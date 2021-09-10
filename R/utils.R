###########################################################
##### Utility Functions (just for internal use)
##### Authors: Jana S. Huisman
###########################################################

#' Change column type
#'
#' \code{.change_column_type} is a helper function
#'   to change the type of a dataframe column
#'
#' @param columns Which column
#' @param dataframe The dataframe
#' @param type Type to change to, default as.numeric
#'
#' @aliases change_col_type
.change_column_type <- function(columns, dataframe, type=as.numeric){
  dataframe[, columns] <- sapply(dataframe[, columns], type)
  return(dataframe)
}

#' Transform simulation output to estimation input
#'
#' \code{.sim_output_to_data} is a helper function
#'   to transform simulation output to estimation input
#'
#' D.t, R.t, T.t are final population sizes;
#' D.0, R.0 are initial pop sizes;
#'
#' @param out Dataframe. Output of lsoda
#' @param parms List. Should contain psi.D, psi.R, psi.T if included
#' @param exp_window Integer. If the growth rate should be calculated from
#' the population sizes, which time step is still in exponential phase.
#'
.sim_output_to_data <- function(out, parms = NULL, exp_window = 20){

  n_sim_steps = dim(out)[1]

  df = data.frame(out[1:4])
  names(df) <- c('t', 'R.t', 'T.t', 'D.t')
  df[c('R.0', 'T.0', 'D.0')] <- out[1, 2:4]

  psi.max = log(out$tot[exp_window]/out$tot[1])/(out$time[exp_window]-out$time[1])
  df['psi.max'] <- psi.max
  
  if(is.null(parms)){
    df[c('psi.R', 'psi.T', 'psi.D')] <- psi.max
  } else {
    df[c('psi.R', 'psi.T', 'psi.D')] <- cbind(rep(parms["psi.R"], n_sim_steps),
                                              rep(parms["psi.T"], n_sim_steps),
                                              rep(parms["psi.D"], n_sim_steps))
  }

  return(df)
}

estimate_conj_from_sim <- function(out, method, parms = NULL, exp_window = 20, 
                                   id_cols = 't', verbose = T){
  data <- .sim_output_to_data(out, parms, exp_window)

  result <- conjugator::estimate_conj_rate(data, method, id_cols, verbose = verbose)
  return(result)
}


estimate_crittime_from_sim <- function(parms, vars, tol_factor = 10, id_cols = 'ID', verbose = T){
  data = data.frame(as.list(c(parms, vars)))
  colnames(data) <- gsub('gamma.t', 'gamma.T', colnames(data)) %>% 
                 gsub('gamma.d', 'gamma.D', .) %>%
                 gsub('r', 'R.0', .) %>%
                 gsub('d', 'D.0', .) %>%
                 gsub('t', 'T.0', .)

  result <- conjugator::estimate_crit_time(data, TRT = NULL, tol_factor, id_cols, verbose = verbose)
  return(result)
}

.rename_fullout <- function(full_out){
  
  colnames(full_out) <- gsub('gamma.t', 'gamma.T', colnames(full_out)) %>% 
    gsub('gamma.d', 'gamma.D', .) %>%
    gsub('r0', 'R.0', .) %>%
    gsub('d0', 'D.0', .) %>%
    gsub('^t$', 'T.t', .) %>%
    gsub('r$', 'R.t', .) %>%
    gsub('d$', 'D.t', .) %>%
    gsub('time', 't', .)
  
  return(full_out)
}
###########################################################
# Fig 2 and Fig 3 of the paper

config_rand_gamma <- function(vars, gamma_mean = 10^-11, gamma_sd = 10^-12){
  config_df <- data.frame()
  for (gamma in rnorm(5, mean = gamma_mean, sd = gamma_sd)){
    parms <- c(psi.R = 1, psi.T = 1, psi.D = 1,
               gamma.t = gamma, gamma.d = gamma,
               q = 1E2)
    parms_rep <- matrix(parms, nrow = dim(vars)[1], ncol=length(parms), byrow=TRUE)
    colnames(parms_rep) <- names(parms)
    
    config_df <- rbind(config_df, cbind(vars, parms_rep))
  }
  return(config_df)
}

create_config_df <- function(vars, growth_rates, gammas){
  
  if(!is.matrix(vars)){
    vars <- as.data.frame(t(vars))
  } else {
    vars <- as.data.frame(vars)
  }
  
  if(!is.matrix(growth_rates)){
    growth_rates <- as.data.frame(t(growth_rates))
  } else {
    growth_rates <- as.data.frame(growth_rates)
  }
  
  if(!is.matrix(gammas)){
    gammas <- as.data.frame(t(gammas))
  } else {
    gammas <- as.data.frame(gammas)
  }
  
  config_df <- tidyr::expand_grid(growth_rates, vars, gammas)
  config_df['q'] = rep(1E2, dim(config_df)[1])
  
  return(config_df)
}

compute_model_scan <- function(config_df, t.vec, output_timepoints = c(21, 41, 81), psi_time = 20, verbose = T){
  full_out <- data.frame()
  crit_times <- data.frame()
  for (row_ind in 1:dim(config_df)[1]){
    # initial population size
    vars <- as.numeric(config_df[row_ind, c('r', 't', 'd', 'c')])
    names(vars) <- c('r', 't', 'd', 'c')
    parms <- as.numeric(config_df[row_ind, c('psi.R', 'psi.T', 'psi.D',
                                             'gamma.t', 'gamma.d', 'q')])
    names(parms) <- c('psi.R', 'psi.T', 'psi.D',
                      'gamma.t', 'gamma.d', 'q')
    other_cols <- config_df[row_ind, setdiff(colnames(config_df), c('r', 't', 'd', 'c',
                                        'psi.R', 'psi.T', 'psi.D',
                                       'gamma.t', 'gamma.d', 'q') )]
    
    out <- integrate_model(vars, t.vec, "model_ESM", parms)
    for (timepoint in output_timepoints){
      full_out <- rbind(full_out, c(parms, r0 = as.numeric(vars['r']),
                                    d0 = as.numeric(vars['d']), out[timepoint,],
                                    psi.max = log(out$tot[psi_time]/out$tot[1])/(out$time[psi_time]-out$time[1]),
                                    other_cols)
      )
      crit_times <- rbind(crit_times, estimate_crittime_from_sim(parms, vars, tol_factor = 10, verbose = verbose))
    }
    
  }
  
  return(list(full_out = full_out, crit_times = crit_times))
}

compute_estimates <- function(full_out, verbose = T){
  full_out <- .rename_fullout(full_out)
  estimates <- conjugator::estimate_conj_rate(full_out, c("T_RT", "TD", "Gama", "ASM", "SM"),
                                              verbose = verbose) %>%
    pivot_wider(id_cols = 'ID', 
                names_from = 'method', values_from = 'estimate')
  
  full_out <- cbind(full_out, estimates)
  return(full_out)
}

compute_foldchange <- function(full_out, verbose = T){
  
  if (!'ASM' %in% colnames(full_out)){
    full_out <- compute_estimates(full_out, verbose = verbose)
  }
  
  full_out['gammaD_fc'] <- full_out['ASM']/full_out['gamma.D']
  full_out['gamma_max_fc'] <- full_out['SM']/full_out['gamma.D']
  
  ref_TD <- exp((full_out['psi.T']- full_out['psi.D'])*full_out['t'])/full_out['D.0']
  
  full_out['TD_fc'] <- full_out['TD']/ref_TD
  
  return(full_out)
}

boxplot_scan <- function(full_out, x_var, y_var, fill_var,
                         logplot = TRUE, add_means = FALSE, plotfile = NULL){
  
  means <- aggregate(full_out[, c(x_var, y_var, fill_var)],
                     by = list(full_out[, x_var], full_out[, fill_var]), mean)
  means[, x_var] <- as.factor(means[, x_var])
  means[, fill_var] <- as.factor(means[, fill_var])
  
  full_out[, x_var] <- as.factor(full_out[, x_var])
  full_out[, fill_var] <- as.factor(full_out[, fill_var])
  
  box_plot <- ggplot(NULL) +
    geom_boxplot(aes(x = full_out[, x_var],
                     y = full_out[, y_var],
                     fill = full_out[, fill_var]),
                 position=position_dodge(1)) +
    labs(x = x_var, y = y_var, fill = fill_var) +
    theme_minimal() + 
    theme(text = element_text(size=15))
  
  if (logplot){
    box_plot <- box_plot +
      scale_y_continuous(trans = 'log', breaks = .base_breaks())
  }
  if (add_means){
    box_plot <- box_plot +
      geom_point(aes(x = means[, x_var], y = means[, y_var], color = means[, fill_var])) +
      labs(x = x_var, y = y_var, color = fill_var)
  }
  
  if (!is.null(plotfile)){
    ggsave(filename = plotfile,
           plot = box_plot, width = 6, height = 4)
  }
  box_plot
}

