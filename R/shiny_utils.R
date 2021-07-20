###########################################################
## Shiny Helper Functions
## Authors: Jana S. Huisman
###########################################################

###########################################################
if.na <- function(value, value.if.na){
  ifelse(is.na(value), value.if.na, value)
}

join_and_overwrite_df <- function(orig_df, add_df, id_cols = "ID"){

  shared_col_all <- intersect(colnames(orig_df), colnames(add_df))
  if (length(id_cols) < 1){
    if (length(shared_col_all)>1){
      id_cols = shared_col_all[1]
    } else {
      stop('No shared columns.')
    }
  }
  shared_col <- setdiff(shared_col_all, id_cols)
  
  joint_df <- dplyr::full_join(orig_df, add_df, by = id_cols, stringsAsFactors = FALSE)
  for (col in shared_col){
    joint_df[, col] <- if.na(joint_df[, paste0(col, '.y')],
                             joint_df[, paste0(col, '.x')])
  }
  joint_df <- joint_df %>%
    select(-contains(".x")) %>%
    select(-contains(".y")) %>%
    select(all_of(id_cols), setdiff(colnames(orig_df), id_cols), setdiff(colnames(add_df), colnames(orig_df)) )
  
  return(joint_df)
}

###########################################################
get.vars <- function(data, id = NULL){
  vars <- cbind(r = data$R.0,
                t = ifelse("T.0" %in% colnames(data), data$T.0, 0),
                d = data$D.0)

  vars[is.na(vars[,"t"]), "t"] <- 0

  if ("C.0" %in% colnames(data)){
    vars <- cbind(vars, data$C.0)
  }

  if (!is.null(id)){
    vars <- vars[id,]
  }
  return(vars)
}

get.growth_rates <- function(data, id = NULL){
  growth_rates <- cbind(psi.R = data$psi.R,
                        psi.T = data$psi.T,
                        psi.D = data$psi.D)
  if (!is.null(id)){
    growth_rates <- growth_rates[id,]
  }
  return(growth_rates)
}

get.pars <-function(data_drt, data_trt, id, input){
  c(get.growth_rates(data_drt, id),
    gamma.t=estimate_conjugation(data_trt, "ASM")[id],
    gamma.d=estimate_conjugation(data_drt, "ASM")[id],
    q = 1E2)
}

get.vars_and_pars <-function(data_drt, data_trt, id, input){
  list(vars = get.vars(data_drt, id), pars = get.pars(data_drt, data_trt, id, input))
}

###########################################################
# Used in Data Analysis Tab
###########################################################

turn_into_df <- function(data_obj){
  if(!is.matrix(data_obj)){
    data_obj <- as.data.frame(t(data_obj))
  } else {
    data_obj <- as.data.frame(data_obj)
  }
  return(data_obj)
}



create_config_df <- function(vars, growth_rates, gammas){
  vars <- turn_into_df(vars)
  growth_rates <- turn_into_df(growth_rates)
  gammas <- turn_into_df(gammas)

  config_df <- tidyr::expand_grid(growth_rates, vars, gammas)
  config_df['q'] = rep(1E2, dim(config_df)[1])

  return(config_df)
}
