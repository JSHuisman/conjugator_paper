###########################################################
## FUNCTIONS FOR SIMULATION MODELS
## Authors: Jana S. Huisman & Sebastian Bonhoeffer
###########################################################

#' The approximate extended Simonsen model (ASM)
#'
#' \code{model_ASM} returns the derivatives
#'   corresponding to the approximate extended Simonsen model
#'
#' @param time Timepoint of evaluation, necessary
#'   for integration using lsoda
#' @param vars Named list of numerics.
#'   Initial population sizes (r, t, d)
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T;
#'   conj rates gamma.t, gamma.d)
#'
#' @family model functions
#' @aliases bonhoeffer_model
#'
#' @export
model_ASM <- function(time, vars, parms) {
  with(as.list(c(parms,vars)),{
    # calculate derivatives
    if (r>=0){
    deriv.r <- psi.R * r - (gamma.t * t + gamma.d * d) *r
    } else {deriv.r <- 0}
    if (t>=0){
    deriv.t <- psi.T * t + (gamma.t * t + gamma.d * d) *r
    } else {deriv.t <- 0}
    if (d>=0){
    deriv.d <- psi.D * d
    } else {deriv.d <- 0}
    # return differentials
    list(c(deriv.r, deriv.t, deriv.d));
  })
}

#' The Simonsen model (SM)
#'
#' \code{model_simonsen} returns the derivatives
#'   corresponding to the Simonsen model
#'
#' @param time Timepoint of evaluation, necessary
#'   for integration using lsoda
#' @param vars Named list of numerics.
#'   Initial population sizes (r, t, d, c)
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T OR
#'   psi.max; conjugation rates gamma.t, gamma.d OR
#'   gamma.max)
#'
#' @family model functions
#' @aliases simonsen_model
#'
#' @export
model_SM <- function(time, vars, parms) {
  with(as.list(c(parms,vars)),{
    if (!"psi.max" %in% names(parms)){
      psi.max = max(c(psi.D, psi.R, psi.T))
    }
    if (!"gamma.max" %in% names(parms)){
      gamma.max = max(c(gamma.d, gamma.t))
    }
    e = 1
    
    # calculate derivatives
    if (r>=0){
    deriv.r <- psi.max * c/(c + q) * r - gamma.max * c/(c + q) * (t + d) * r
    } else {deriv.r <- 0}
    if (t>=0){
    deriv.t <- psi.max * c/(c + q) * t + gamma.max * c/(c + q) * (t + d) * r
    } else {deriv.t <- 0}
    if (d>=0){
    deriv.d <- psi.max * c/(c + q) * d
    } else {deriv.d <- 0}
    if (c>0){
      deriv.c <- -psi.max * c/(c + q) * (r + d + t)*e
    } else{
      deriv.c <- 0
    }
    # return differentials
    list(c(deriv.r, deriv.t, deriv.d, deriv.c));
  })
}

#' The extended Simonsen model (ESM)
#'
#' \code{model_simonsen_extended} returns the derivatives
#'   corresponding to the extended Simonsen model
#'
#' @param time Timepoint of evaluation, necessary
#'   for integration using lsoda
#' @param vars Named list of numerics.
#'   Initial population sizes (r, t, d, c)
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T;
#'   conjugation rates gamma.t, gamma.d)
#'
#' @family model functions
#' @aliases extended_simonsen_model
#'
#' @export
model_ESM<- function(time, vars, parms) {
  with(as.list(c(parms,vars)),{
    # calculate derivatives
    if (r>=0){
      deriv.r <- psi.R * c/(q +c) * r - (gamma.d* c/(q +c) * d  +  gamma.t * c/(q +c) * t) * r
    } else {
      deriv.r <- 0
    }
    if (t>=0){
    deriv.t <- psi.T * c/(q +c) * t + (gamma.d* c/(q +c) * d  +  gamma.t * c/(q +c) * t) * r
    } else {deriv.t <- 0}
    if (d>=0){
    deriv.d <- psi.D * c/(q +c) * d
    } else {deriv.d <- 0}
    if (c>=0){
      deriv.c <- -(psi.R * r + psi.T * t + psi.D * d)* c/(q +c)
    } else {
      deriv.c <- 0
    }

    # return differentials
    list(c(deriv.r, deriv.t, deriv.d, deriv.c));
  })
}

#' Integrating the models
#'
#' \code{integrate_model} solves the ODE model
#'   numerically
#'
#' @param vars Named list of numerics.
#'   Initial population sizes; depending on the model
#'   ASM: (r, t, d); SM: (r, t, d, c)
#' @param t.vec Vector of evaluation timepoints
#' @param model The selected model as a string
#'   Choice from model_ASM, model_SM,
#'   model_ESM
#' @param parms Named list of numerics.
#'   Parameters (growth rates psi.R, psi.D, psi.T;
#'   conjugation rates gamma.t, gamma.d)
#'
#' @aliases simulate_model
#'
#' @export
integrate_model <- function(vars, t.vec, model, parms){

  if(model == "model_ASM"){
    stopifnot(c('r', 't', 'd') %in% names(vars))
    vars <- vars[c('r', 't', 'd')]
  } else if (model == "model_SM"){
    stopifnot(c('r', 't', 'd', 'c') %in% names(vars))
    vars <- vars[c('r', 't', 'd', 'c')]
  }

  out <- as.data.frame(deSolve::lsoda(y = vars, times = t.vec,
                                      func = get(model), parms = parms))

  # remove negative entries from the recipient population
  #out <- .clean_negative_Rpop(out)
  # add total population size to the output
  out$tot <- apply(out[, 2:4], 1, sum)
  return(out)
}

