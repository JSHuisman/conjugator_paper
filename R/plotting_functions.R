###########################################################
## FUNCTIONS FOR PLOTTING
## Authors: Jana S. Huisman & Sebastian Bonhoeffer
###########################################################

#' Plotting overall dynamics with ggplot
#'
#' \code{plot_dynamics_ggplot} plots the overall dynamics
#'
#' @param out Output of a simulated model
#' @param crit_times Named list of numerics.
#'   Critical times, output of function time_crit()
#'
#' @family plot functions
plot_dynamics_ggplot <- function(out, crit_times = NULL){

  # this is to avoid the  "transformation introduced
  # infinite values in continuous y-axis" error
  # Think of it as "detection limit"
  lowest_y_axis_val <- 1E-6
  out[out == 0] <- lowest_y_axis_val

  out_molten <- reshape2::melt(out[, c('time', 'r', 'd', 't')], id.vars = 'time')

  # plot the dynamics
  dyn_plot <- ggplot2::ggplot(NULL, ggplot2::aes(x = out_molten$time,
                                        y = out_molten$value)) +
              ggplot2::geom_line(ggplot2::aes(color = out_molten$variable))
  # add critical times to the plot
  if(is.list(crit_times)){
    dyn_plot <- add_tcrit_to_plot(dyn_plot, crit_times,
                                  y_pos = lowest_y_axis_val)
  }
  # add scale, labels, legend
  dyn_plot  +
    ggplot2::scale_y_continuous(trans = 'log', 
                                breaks = .base_breaks()) +
    coord_cartesian(ylim = c(lowest_y_axis_val, max(out$tot)),
                    xlim = range(out$time)) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15),
                                            legend.position = c(0.7, 0.2)) +
    #annotation_logticks(base = 10, sides='l') +
    ggplot2::labs(x = "Time (h)", y = "Population size (CFU/ml)", color = "Population") +
    ggplot2::scale_colour_manual(values = c('red', 'blue', 'green'),
                        breaks = levels(out_molten$variable),
                        labels = c('Recipients', 'Donors', 'Transconjugants'))
}

.base_breaks <- function(n = 5){
  function(x) {
    grDevices::axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    #pretty(gamma_max_growth, n = 5)
  }
}
###########################################################
#' Add critical times to a ggplot object
#'
#' \code{add_tcrit_to_plot} plots the overall dynamics
#'
#' @param plot_obj A plot object, that the plotting
#'   of the critical times will be added to
#' @param crit_times Named list of numerics.
#'   Critical times, output of function time_crit()
#' @param y_pos A numeric. Position of the tcrit label
#'
#' @family plot functions
add_tcrit_to_plot <- function(plot_obj, crit_times, y_pos = 0.001, id_cols = 'ID'){
  crit_df <- reshape2::melt(crit_times[c(id_cols, 'tcrit1', 'tcrit2', 'tcrit3')], id.vars = id_cols)

  plot_obj <- plot_obj +
    ggplot2::geom_vline(xintercept = crit_df$value, linetype = 'dotted') +
    ggrepel::geom_text_repel(ggplot2::aes(x = crit_df$value, y = y_pos,
                                 label = c('t[c1]', 't[c2]', 't[c3]')),
                             parse = TRUE)
  return(plot_obj)
}

###########################################################
#' Plotting estimated growth together with the model
#'
#' \code{plot_growth_rate_ggplot} plots the growth rate
#'   estimate together with input values of growth
#'   rates of D,R and T
#'
#' @param out Simulation output
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T;
#'   conjugation rates gamma.t, gamma.d)
#' @param crit_times A named list. Output of time_crit()
#'
#' @family plot functions
plot_growth_rate_ggplot <- function(out, parms, crit_times = NULL){
  sim_data <- .sim_output_to_data(out, parms, exp_window = max(out$time))

  # psi.max will be NaN at time = 0
  growth_plot <- ggplot2::ggplot(NULL, ggplot2::aes(x = out$time,
                                           y = sim_data$psi.D)) +
    ggplot2::geom_point(ggplot2::aes(color = 'c1')) +
    ggplot2::labs(x = "Time (h)", y = "Growth rate (1/h)") +
    coord_cartesian(xlim = range(out$time)) +
    ggplot2::scale_color_manual(name = 'Growth rate', breaks = 'c1',
                       values = 'red', labels = expression(psi[max])) +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15),
                          legend.position = c(0.8, 0.2)) +
    ggplot2::geom_hline(yintercept = c(parms["psi.R"], parms["psi.T"], parms["psi.D"]),
                        linetype = 'dashed') +
    ggrepel::geom_text_repel(ggplot2::aes(x = 20,
                                 y = c(parms["psi.R"], parms["psi.T"], parms["psi.D"]),
                                 label = c('psi[R]', 'psi[T]', 'psi[D]')),
                             parse = TRUE)
  # add critical times to this plot
  if(is.list(crit_times)){
    growth_plot <- add_tcrit_to_plot(growth_plot, crit_times,
                                     y_pos = min(parms["psi.R"],
                                                 parms["psi.T"],
                                                 parms["psi.D"]))
  }
  # plot the final result
  growth_plot
}

###########################################################
#' Plotting estimated conjugation rates together with the model
#'
#' \code{plot_conj_rate_ggplot} plots the conjugation rate
#'   estimate together with input values of conjugation rates
#'
#' This function sometimes gives warnings about missing
#' values when the R population becomes 0 at some point,
#' and correspondingly the Bonhoeffer approximation becomes
#' negative
#'
#' @param out Simulation output
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T;
#'   conjugation rates gamma.t, gamma.d)
#' @param crit_times A named list. Output of time_crit()
#'
#' @family plot functions
plot_conj_rate_ggplot <- function(out, parms, crit_times = NULL){
  SM_conj_result <- estimate_conj_from_sim(out, "SM", parms, exp_window = max(out$time))
  ASM_conj_result <- estimate_conj_from_sim(out, "ASM", parms, exp_window = max(out$time))

  lowest_y_axis_val <- min(c(parms["gamma.d"], parms["gamma.t"]))*10**-2
  highest_y_axis_val <- max(c(parms["gamma.d"], parms["gamma.t"]))*10**2
  
  gamma_plot <- ggplot2::ggplot(NULL) +
    ggplot2::geom_point(ggplot2::aes(x = out$time[!is.na(SM_conj_result$estimate)],
                            y = SM_conj_result$estimate[!is.na(SM_conj_result$estimate)],
                            color = 'c1')) +
    ggplot2::geom_point(ggplot2::aes(x = out$time[!is.na(ASM_conj_result$estimate)],
                            y = ASM_conj_result$estimate[!is.na(ASM_conj_result$estimate)],
                            color = 'c2')) +
    ggplot2::scale_y_continuous(trans = 'log',
                       breaks = .base_breaks()) +
    coord_cartesian(ylim = c(lowest_y_axis_val, highest_y_axis_val),
                    xlim = range(out$time)) +
    ggplot2::scale_color_manual(name = 'Conjugation rate',
                       breaks = c('c1', 'c2'),
                       values = c('blue', 'lightblue'),
                       labels = c(expression(gamma[max]), expression(gamma[D]))) +
    ggplot2::labs(x = "Time (h)", y = "Conjugation rate (mL/(CFU*h))") +
    ggplot2::theme_light() + ggplot2::theme(text = ggplot2::element_text(size=15),
                         legend.position = c(0.7, 0.8)) +
    ggplot2::geom_hline(yintercept = c(parms["gamma.d"], parms["gamma.t"]),
                        linetype = 'dashed') +
    ggrepel::geom_text_repel(ggplot2::aes(x = max(out$time)-2,
                                          y = c(parms["gamma.d"], parms["gamma.t"])),
                             label = c('gamma[D]', 'gamma[T]'), parse = TRUE)
  # add critical times
  if(is.list(crit_times)){
    gamma_plot <- add_tcrit_to_plot(gamma_plot, crit_times,
                                    y_pos = 10**-2*min(c(parms["gamma.d"],
                                                         parms["gamma.t"])))
  }
  # plot the final result
  gamma_plot
}


#################

plot_crittimes_sweep <- function(data, id_cols = "ID"){
  # apply this to the summarised data [per groupby variable]
  crit_result = scan_crit_time(data, tol_factor = 10, id_cols = id_cols, mult_seq = 10**seq(-2, 4, 0.1))

  crit_result <- crit_result %>%
    mutate_at(c("gamma.T", "min_tcrit"), as.double) %>%
    mutate_at(c(id_cols), as.factor)
  
  if(length(id_cols)>1){
    col_id <- id_cols[1]
  } else {col_id <- id_cols}
  
  crittimes_plot <- ggplot(NULL) +
    geom_line(aes(x = crit_result[["gamma.T"]]/crit_result[["gamma.D"]],
                     y = crit_result[["min_tcrit"]],
                     color = crit_result[[col_id]] )) +
    scale_x_continuous(trans = 'log', breaks = .base_breaks()) +
    coord_cartesian(ylim = c(0, ceiling(max(crit_result[["min_tcrit"]])))) +
    labs(x = bquote("Ratio"~gamma[T]/gamma[D]), y = "Critical time", color = col_id) +
    theme_light() + theme(text = element_text(size=15))

  return(crittimes_plot)
}


###########################################################
#' Create boxplot of the SM and ASM estimates
#'
#' \code{create_boxplot}
#'
#' @param data Dataframe with estimates of SM
#'   and ASM gamma rates
#'
#' @family plot functions
create_boxplot <- function(data){
  plot_data <- reshape2::melt(
    as.data.frame(data[, c('ID', 'gamma.max', 'gamma.D')]),
    id.vars='ID')

  p <- ggplot2::ggplot(NULL) +
    ggplot2::geom_boxplot(ggplot2::aes(x=plot_data$variable,
                                       y=plot_data$value))
  return(p)
}


boxplot_scan <- function(full_out, x_var, y_var, fill_var,
                         logplot = TRUE, add_means = FALSE, plotfile = NULL){
  
  if(!all(c(x_var, y_var, fill_var) %in% colnames(full_out))){
    return(ggplot(NULL))
  }
  
  if(add_means){
    means <- aggregate(full_out[, c(x_var, y_var, fill_var)],
                     by = list(full_out[, x_var], full_out[, fill_var]), mean)
    means[, x_var] <- as.factor(means[, x_var])
    means[, fill_var] <- as.factor(means[, fill_var])
  }

  full_out[[fill_var]] <- as.factor(signif(full_out[[fill_var]], digits = 3))

  box_plot <- ggplot(NULL) +
    geom_boxplot(aes(x = full_out[[x_var]],
                     y = full_out[[y_var]],
                     fill = full_out[[fill_var]]),
                 position=position_dodge(1)) +
    labs(x = x_var, y = y_var, fill = fill_var) +
    theme_light() + theme(text = element_text(size=15))

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
  return(box_plot)
}
