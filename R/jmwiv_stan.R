#' Bayesian multivariate joint model with within-individual variability-based association
#'
#' @param formulaLong1 Formula for the first longitudinal marker.
#' @param formulaLong2 Formula for the second longitudinal marker.
#' @param dataLong Longitudinal dataset.
#' @param id_var Character string for name of ID variable in `dataLong`.
#' @param time_varLong Character string for name of time variable in `dataLong`.
#' @param formulaEvent Formula for the time-to-event submodel (without association terms)
#' @param dataEvent Event dataset
#' @param assoc Type of association (either 'RE' for random effects, 'LP' for linear predictor, 'CV' for current value).
#' @param a_K Number of association terms
#' @param basehaz Type of baseline hazard function (only 'bs' is implemented)
#' @param basehaz_aux List of parameters for the baseline hazard function.
#' @param qnodes Number of quadrature nodes.
#' @param prior_list List of prior arguments.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @examples
#' library(rstanarm)
#' f_logBili <- list(formula(logBili ~ sex + trt + year + (1 | id)),
#'                   formula(sigma ~ sex + trt + year + (1 | id)))
#' f_albumin <- list(formula(albumin ~ sex + trt + year + (1 | id)),
#'                   formula(sigma ~ sex + trt + year + (1 | id)))
#' f_Event <- survival::Surv(futimeYears, death) ~ sex + trt
#' prior_input <- list(b_prior_scale = rep(1, 4),
#'                     b_prior_df = rep(1, 4),
#'                     b_prior_regularization = 5)
#' fit_jmwiv <- jmwiv_stan(formulaLong1 = f_logBili,
#'                          formulaLong2 = f_albumin,
#'                          dataLong = pbcLong,
#'                          id_var = "id",
#'                          time_varLong = "year",
#'                          formulaEvent = f_Event,
#'                          dataEvent = pbcSurv,
#'                          assoc = "CV",
#'                          a_K = 4,
#'                          basehaz = "bs",
#'                          basehaz_aux = list(df = 6,
#'                                             knots = NULL,
#'                                             degree = 3),
#'                          qnodes = 15L,
#'                          prior_list = prior_input,
#'                          cores = 2,
#'                          chains = 2,
#'                          warmup = 1000,
#'                          iter = 2000,
#'                          control = list(max_treedepth = 15))
jmwiv_stan <- function(formulaLong1,
                       formulaLong2 = NULL,
                       dataLong,
                       id_var,
                       time_varLong,
                       formulaEvent,
                       dataEvent,
                       assoc,
                       a_K = 4,
                       basehaz = "bs",
                       basehaz_aux = list(df = 6,
                                          knots = NULL,
                                          degree = 3),
                       qnodes = 15L,
                       prior_list = NULL,
                       ...){
  standata <- make_standata_jmwiv(formulaLong1 = formulaLong1,
                                  formulaLong2 = formulaLong2,
                                  dataLong = dataLong,
                                  id_var = id_var,
                                  time_varLong = time_varLong,
                                  formulaEvent = formulaEvent,
                                  dataEvent = dataEvent,
                                  assoc = assoc,
                                  a_K = a_K,
                                  basehaz = basehaz,
                                  basehaz_aux = basehaz_aux,
                                  qnodes = qnodes,
                                  prior_list = prior_list)
  out <- rstan::sampling(stanmodels$jmwiv, data = standata, ...)
  return(out)
}
