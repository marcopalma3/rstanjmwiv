#' Print random effect correlation matrix
#'
#' @param fit An object of class `stanfit` returned by `rstan::sampling`
#'
#' @return A correlation matrix.
#' @export
ranef_corr_jmwiv <- function(fit){
  covmat <- matrix(0, ncol = 4, nrow = 4)
  covmat[lower.tri(covmat, diag = TRUE)] <- summary(fit, pars = "b_cov")$summary[,1]
  covmat[upper.tri(covmat, diag = FALSE)] <- t(covmat)[upper.tri(covmat, diag = FALSE)]
  cov2cor(covmat)
}



#' Posterior predictive checks for the longitudinal submodel of JM-WIV
#'
#' @param list_of_draws Output from `rstan::extract(fit)`.
#' @param data Longitudinal dataset.
#' @param x Time variable in the longitudinal dataset.
#' @param idvar ID variable in the longitudinal dataset.
#'
#' @return A dataframe with posterior predictive draws for each longitudinal biomarker.
#' @export
#'
#' @import dplyr
#'
#' @examples ppcheck(list_of_draws = rstan::extract(fit_jmwiv), data = pbcLong,
#'                   x = "year", idvar = "id")
ppcheck <- function(list_of_draws,
                    data,
                    x,
                    idvar){
  x <- select(data, all_of(x))
  # id <- select(data, all_of(id))

  ind <- sample(1:dim(list_of_draws[[1]]), size = 1)
  sim_ranef <- list_of_draws$b_mat[ind,,]


  y1mu_eta <- with(list_of_draws,
                   y1mu_Intercept[ind] +
                     y1mu_beta[ind]*x +
                     model.matrix(~ 0 + as.factor(get(idvar)),
                                  data= data) %*% sim_ranef[,1]) %>% .[,1]

  y1sigma_eta <- with(list_of_draws,
                      y1sigma_Intercept[ind] +
                        y1sigma_beta[ind]*x +
                        model.matrix(~ 0 + as.factor(get(idvar)),
                                     data= data) %*% sim_ranef[,2]) %>% .[,1]

  y2mu_eta <- with(list_of_draws,
                   y2mu_Intercept[ind] +
                     y2mu_beta[ind]*x +
                     model.matrix(~ 0 + as.factor(get(idvar)),
                                  data= data) %*% sim_ranef[,3]) %>% .[,1]

  y2sigma_eta <- with(list_of_draws,
                      y2sigma_Intercept[ind] +
                        y2sigma_beta[ind]*x +
                        model.matrix(~ 0 + as.factor(get(idvar)),
                                     data= data) %*% sim_ranef[,4]) %>% .[,1]

  y1 <- rnorm(length(with(data, get(idvar))), y1mu_eta, exp(y1sigma_eta))
  y2 <- rnorm(length(with(data, get(idvar))), y2mu_eta, exp(y2sigma_eta))

  return(data.frame(y1,y2))
}






#' Simulate data from a joint model from Cox specification with baseline hazard following exponential or Gompertz distribution and MELS specification for longitudinal submodel
#'
#' This function returns a list with longitudinal and survival data.
#' The baseline hazard function is assumed to be constant over time.
#' At the moment only Cox-Exponential for 'RE' is available, while for 'LP' you can have Cox-Exp or Cox-Gompertz.
#' @param seed seed used for simulation
#' @param nsub number of subjects
#' @param assoc either random effects (RE) or linear predictor (LP) - current value association not available
#' @param g_shape shape parameter for the Gompertz distribution (for g_shape=0, this reduces to Cox-exponential)
#' @param lambda baseline hazard function (constant over time)
#' @param e_beta regression coefficient in Cox model for a baseline, binary variable (e.g. treatment)
#' @param a_beta regression coefficients in Cox model for the random intercept of the mean and SD of y1 and mean and SD of y2, respectively
#' @param distObs time between consecutive visits
#' @param max_time censoring time
#' @param ranef_covmat covariance matrix for random effects
#' @param y1mu_fixed fixed effects for y1_mu
#' @param y1sigma_fixed fixed effects for y1_sigma
#' @param y2mu_fixed fixed effects for y2_mu
#' @param y2sigma_fixed fixed effects for y2_sigma
#'
#' @return A list with a simulated dataset for the longitudinal data and a dataset for the event data.
#' @export
#'
#' @import dplyr
simulate_jmwiv <- function(seed,
                           nsub,
                           assoc = "RE",
                           lambda = 1,
                           g_shape = 0,
                           e_beta = c(-0.5, 0),
                           a_beta = c(-0.5, -0.2, -0.4, 0.3),
                           distObs  = 0.5,
                           max_time = 10,
                           ranef_covmat = c(3,1,2,0.2) * diag(4),
                           y1mu_fixed = c(-0.8, 2),   #first value is intercept, second value is for visit_times, third for baseline covariate long_cov
                           y1sigma_fixed = c(0.2, -0.5),
                           y2mu_fixed = c(1.2, -0.3),
                           y2sigma_fixed = c(-0.6, 1.4)
){

  if(assoc != "LP" && assoc != "RE") stop("Association must be either 'RE' or 'LP'.")

  set.seed(seed)

  ### Generate IDs #################################

  ID <- seq.int(1, nsub) %>% str_pad(., 4, pad = "0") %>% as.factor() ###this covers up to nsub=9999
  binary_cov <- rbinom(n=nsub, size=1, prob=0.5)
  norm_cov <- rnorm(n = nsub, 0, 1)

  ### Generate random effects ######################

  ranef_dat <- MASS::mvrnorm(n = nlevels(ID),
                             mu = rep(0, nrow(ranef_covmat)),
                             Sigma = ranef_covmat)


  ### Generate survival times ######################

  if(assoc == "RE"){

    e_time <- - log(runif(nlevels(ID)))/(lambda * exp(e_beta[1] * binary_cov +
                                                      e_beta[2] * norm_cov +
                                                        ranef_dat %*% a_beta))
  }else{
    A <- a_beta[1]*(y1mu_fixed[1] + ranef_dat[, 1]) +
      a_beta[2]*(y1sigma_fixed[1] + ranef_dat[, 2]) +
      a_beta[3]*(y2mu_fixed[1] + ranef_dat[, 3]) +
      a_beta[4]*(y2sigma_fixed[1] + ranef_dat[, 4])

    B <- a_beta[1]*y1mu_fixed[2] +
      a_beta[2]*y1sigma_fixed[2] +
      a_beta[3]*y2mu_fixed[2] +
      a_beta[4]*y2sigma_fixed[2] +
      g_shape

  suppressWarnings(e_time <- log(1 - B*log(runif(nlevels(ID)))/(lambda*exp(A + e_beta[1] * binary_cov + e_beta[2] * norm_cov)))/B)
  }

  ### Generate censoring times ######################

  c_time <- runif(nlevels(ID), 0, max_time)
  e_status <- ifelse(c_time <= e_time| is.nan(e_time), 0, 1)
  e_time <- ifelse(e_status == 0, c_time, e_time)

  dataEvent <- data.frame(ID, e_time, e_status, binary_cov,
                          norm_cov,
                          ranef_dat,
                          row.names = NULL)

  ### Generate visit times at regular times from 0 to time ######################

  dataLong <- apply(dataEvent, MARGIN = 1,
                    function(x) data.frame(ID = x[1],
                                           visit_times = seq(0, x[2], by = distObs),
                                           row.names = NULL)) %>%
    rlist::list.stack()


  ### Generate y1 and y2 following the MELSM model ######################


  mu1 <- y1mu_fixed[1] +
    pull(y1mu_fixed[2] * dataLong["visit_times"]) +
    (model.matrix(~ 0 + ID, data= dataLong) %*% ranef_dat[,1]) %>%
    as.numeric()

  sigma1 <- exp(y1sigma_fixed[1] +
                  y1sigma_fixed[2]*dataLong["visit_times"] +
                  model.matrix(~ 0 + ID, data= dataLong) %*% ranef_dat[,2]) %>%
    pull()

  mu2 <- y2mu_fixed[1] +
    pull(y2mu_fixed[2] * dataLong["visit_times"]) +
    (model.matrix(~ 0 + ID, data= dataLong) %*% ranef_dat[,3])%>%
    as.numeric()

  sigma2 <- exp(y2sigma_fixed[1] +
                  y2sigma_fixed[2] * dataLong["visit_times"] +
                  model.matrix(~ 0 + ID, data= dataLong) %*% ranef_dat[,4]) %>%
    pull()


  dataLong$y1 <- rnorm(nrow(dataLong), mean = mu1, sd = sigma1)
  dataLong$y2 <- rnorm(nrow(dataLong), mean = mu2, sd = sigma2)
  dataLong <- cbind(dataLong,
                    (model.matrix(~ 0 + ID, data= dataLong) %*% ranef_dat))

  list(dataLong = dataLong, dataEvent = dataEvent)
}








#' Extract covariates for submodel
#'
#' Internal function
#' @param standata Stan data input
#' @param param Parameter
#' @param times Quadrature times
#' @param ind Index
#'
#' @import dplyr
#'
#' @export
#' @return A dataframe to be used within `predict_long_qpts`.
prepare_standata_prediction <- function(standata, param, times, ind){
  a1 <- standata[[paste0(param, "_X")]] %>%
    as.data.frame()
  if(ncol(a1)>1){
    a1 <- select(a1, -standata$time_varLong) %>%
      filter(row_number() %in% ind) %>%
      unique() %>%
      data.frame("time" = times, ., row.names = NULL) %>%
      rename(!!standata$time_varLong := "time")
  }else{
    data.frame("time" = times, row.names = NULL) %>%
      rename(!!standata$time_varLong := "time")}
}




#' Predict longitudinal biomarker value at quadrature points
#'
#' Internal function
#' @param standata Stan data input
#' @param times Quadrature times
#' @param ind Index for the individual
#' @param psb Dataset for individual
#' @param assoc_code Association code (0 for 'RE', 1 for 'LP', 2 for 'CV').
#'
#' @import dplyr
#' @export
#' @return Prediction of longitudinal biomarkers at quadrature points (to be used in `pred_plot_jmwiv`).
predict_long_qpts <- function(standata, stanfit, times = qpts, ind, id_selected, psb, assoc_code){


    newdat_y1mu <- prepare_standata_prediction(standata,
                                               param = "y1mu",
                                               times,
                                               ind)

    newdat_y1sigma <- prepare_standata_prediction(standata,
                                                  param = "y1sigma",
                                                  times,
                                                  ind)

    newdat_y2mu <- prepare_standata_prediction(standata,
                                               param = "y2mu",
                                               times,
                                               ind)

    newdat_y2sigma <- prepare_standata_prediction(standata,
                                                  param = "y2sigma",
                                                  times,
                                                  ind)

    psb <- as.matrix(stanfit)

    linpred <- array(NA, dim = c(nrow(psb), length(times), 4))
    dimnames(linpred)[[3]] <- c("mu1","sigma1","mu2","sigma2")

    linpred[,,"mu1"] <- c(psb[,"y1mu_Intercept"]) +
      tcrossprod(psb[,grep("^y1mu_beta", colnames(psb))], as.matrix(newdat_y1mu)) +
      as.numeric(psb[, c(paste0("b_mat[", id_selected,",1]"))])

    linpred[,,"sigma1"] <- c(psb[,"y1sigma_Intercept"]) +
      tcrossprod(psb[,grep("^y1sigma_beta", colnames(psb))], as.matrix(newdat_y1sigma)) +
      as.numeric(psb[, c(paste0("b_mat[", id_selected,",2]"))])

    linpred[,,"mu2"] <- c(psb[,"y2mu_Intercept"]) +
      tcrossprod(psb[,grep("^y2mu_beta", colnames(psb))], as.matrix(newdat_y2mu)) +
      as.numeric(psb[, c(paste0("b_mat[", id_selected,",3]"))])

    linpred[,,"sigma2"] <- c(psb[,"y2sigma_Intercept"]) +
      tcrossprod(psb[,grep("^y2sigma_beta", colnames(psb))], as.matrix(newdat_y2sigma)) +
      as.numeric(psb[, c(paste0("b_mat[", id_selected,",4]"))])

    if(assoc_code == 2){
      linpred[,,"sigma1"] <- exp(linpred[,,"sigma1"])
      linpred[,,"sigma2"] <- exp(linpred[,,"sigma2"])
    }

    linpred
}




#' Posterior prediction from JM-WIV for an individual
#'
#' Produce posterior predictions for one or more subjects (code is slow, please use only for up to 3 subjects).
#' @param stanfit Fit from `jmwiv_stan.R`
#' @param standata Stan data input
#' @param id_selected ID for one or more subject (if more than one value, you need to sort them)
#' @param xlim_input Range of x-axis
#' @param start_timezero Start prediction at time 0 or at last longitudinal observation?
#' @param var_labels Labels for outcomes in the plot
#' @param plot_jm Plot the predictions?
#'
#' @return Plot of prediction for one or more subjects.
#' @export
#'
#' @import dplyr ggplot2
#'
#' @examples pred_plot_jmwiv(stanfit = fit_jmwiv,
#' standata = standata_jmwiv_pbc,
#' var_labels = c("logBilirubin", "Albumin", "Survival probability"),
#' id_selected = c(1,5), plot_jm = FALSE) +
#'  labs(x = "Age")
pred_plot_jmwiv <- function(stanfit,
                            standata,
                            id_selected,
                            start_timezero = FALSE,
                            var_labels,
                            xlim_input = attr(standata$basehaz_X, "Boundary.knots"),
                            plot_jm = TRUE){

  ind <- which(standata$y1_Z_id %in% id_selected)

  psb <- as.matrix(stanfit)

  predicted_values_y1    <- matrix(NA, nrow = nrow(psb), ncol = length(ind))
  predicted_values_y2    <- matrix(NA, nrow = nrow(psb), ncol = length(ind))


  for(j in 1:nrow(psb)) {

    mu1 <- c(psb[j,"y1mu_Intercept"]) +
      as.matrix(sweep(matrix(standata$y1mu_X[ind,], nrow = length(ind)), 2, STATS = standata$y1mu_Xbar, FUN = "+")) %*% psb[j,grep("^y1mu_beta", colnames(psb))] +
      as.numeric(psb[j, c(paste0("b_mat[", standata$y1_Z_id[ind],",1]"))])

    sigma1 <- c(psb[j,"y1sigma_Intercept"]) +
      as.matrix(sweep(matrix(standata$y1sigma_X[ind,], nrow = length(ind)), 2, STATS = standata$y1sigma_Xbar, FUN = "+")) %*% psb[j,grep("^y1sigma_beta", colnames(psb))] +
      as.numeric(psb[j, c(paste0("b_mat[", standata$y1_Z_id[ind],",2]"))])

    mu2 <- c(psb[j,"y2mu_Intercept"]) +
      as.matrix(sweep(matrix(standata$y2mu_X[ind,], nrow = length(ind)), 2, STATS = standata$y2mu_Xbar, FUN = "+")) %*% psb[j,grep("^y2mu_beta", colnames(psb))] +
      as.numeric(psb[j, c(paste0("b_mat[", standata$y2_Z_id[ind],",3]"))])

    sigma2 <- c(psb[j,"y2sigma_Intercept"]) +
      as.matrix(sweep(matrix(standata$y2sigma_X[ind,], nrow = length(ind)), 2, STATS = standata$y2sigma_Xbar, FUN = "+")) %*% psb[j,grep("^y2sigma_beta", colnames(psb))] +
      as.numeric(psb[j, c(paste0("b_mat[", standata$y2_Z_id[ind],",4]"))])

    predicted_values_y1[j, ]   <- rnorm(length(ind), mu1, exp(sigma1))
    predicted_values_y2[j, ]   <- rnorm(length(ind), mu2, exp(sigma2))
  }

  summary_pred_y1 <-  apply(predicted_values_y1, 2, function(x) c(mean(x),quantile(x, probs = c(0.025,0.975)))) %>% t %>%
    data.frame(., "outcome" = 1)
  summary_pred_y2 <-  apply(predicted_values_y2, 2, function(x) c(mean(x),quantile(x, probs = c(0.025,0.975)))) %>% t %>%
    data.frame(., "outcome" = 2)


  colnames(summary_pred_y1) <- c("Estimate", "Q2.5", "Q97.5", "outcome")
  colnames(summary_pred_y2) <- c("Estimate", "Q2.5", "Q97.5", "outcome")


  data_pred <- data.frame("time" = c(with(standata, y1sigma_X[ind, time_varLong] + y1sigma_Xbar[time_varLong]),
                                     with(standata, y2sigma_X[ind, time_varLong] + y2sigma_Xbar[time_varLong])),
                          rbind(summary_pred_y1, summary_pred_y2),
                          "observed" = with(standata, c(y1[ind], y2[ind])),
                          "id" = with(standata, c(y1_Z_id[ind], y1_Z_id[ind])),
                          "time_vline" = NA
  )

  max_time <- attr(standata$basehaz_X, "Boundary.knots")[2]

  last_Longtime <- case_when(start_timezero == FALSE ~ data_pred %>%
                               group_by(id) %>%
                               slice_tail(n = 1) %>%
                               pull(time), TRUE ~ rep(0, length(id_selected)))

  posteriordraws_survcoef <- psb[, grep("^[ea]_beta\\[", colnames(psb))]

  summary_Survsub_list <- vector("list", length(id_selected))

  quad <- mvQuad::QuadRules[["GKr"]][[29]]

  for(idsel in 1:length(id_selected)){

    qwts <- (max_time - last_Longtime[idsel])*quad$w
    qpts <- last_Longtime[idsel] + (max_time - last_Longtime[idsel])*quad$n

    bs_basis <- splines::bs(qpts,
                            df = standata$basehaz_df,
                            knots = NULL,
                            degree = 3,
                            Boundary.knots = c(0, max_time),
                            intercept = TRUE)

    posteriordraws_logbasehaz <- tcrossprod(bs_basis, psb[, grep("^e_aux\\[", colnames(psb))]) %>%
      sweep(., 2, STATS = psb[, "e_Intercept"], FUN = "+")

    haz_subj <- posteriordraws_logbasehaz


    if(standata$assoc_code == 0){
      ## RE association


      newdat <- cbind(matrix(rep(standata$e_modelmat[id_selected[idsel], ], each = nrow(psb)),
                             nrow = nrow(psb), ncol = ncol(standata$e_modelmat)),
                      psb[, grep(paste0("b_mat\\[", id_selected[idsel]), colnames(psb))]) %>%
        {array(., dim = c(1, dim(.)))}

      for(j in 1:nrow(posteriordraws_survcoef)){
        haz_subj[, j] <- exp(posteriordraws_logbasehaz[,j] +
                               c(newdat[, j, ] %*% posteriordraws_survcoef[j, ])) ##sum(posteriordraws_survcoef[j, ]*newdat[j, ]))  ###can I rewrite this as newdat[1, j, ] %*% posteriordraws_survcoef[j, ]
      }

    }else if(standata$assoc_code == 1){

      ind_1 <- which(standata$y1_Z_id %in% id_selected[idsel])

      newdat <- predict_long_qpts(standata, stanfit, times = qpts, ind_1, id_selected[idsel], psb, assoc_code = 1) %>%
        plyr::aaply(., 2,
                    function(x){cbind(
                      matrix(rep(standata$e_modelmat[id_selected[idsel], ],
                                 each = nrow(psb)),
                             nrow = nrow(psb), ncol = ncol(standata$e_modelmat)), x)}
        )

      for(j in 1:nrow(posteriordraws_survcoef)){
        haz_subj[, j] <- exp(posteriordraws_logbasehaz[,j] + newdat[, j, ] %*% posteriordraws_survcoef[j, ])
      }


      # stop("Not yet implemented for LP association")
    }else if(standata$assoc_code == 2){
      ## CV association
      ind_1 <- which(standata$y1_Z_id %in% id_selected[idsel])

      newdat <- predict_long_qpts(standata, stanfit, times = qpts, ind_1, id_selected[idsel], psb, assoc_code = 2) %>%
        plyr::aaply(., 2,
                    function(x){cbind(
                      matrix(rep(standata$e_modelmat[id_selected[idsel], ],
                                 each = nrow(psb)),
                             nrow = nrow(psb), ncol = ncol(standata$e_modelmat)), x)}
        )

      for(j in 1:nrow(posteriordraws_survcoef)){
        haz_subj[, j] <- exp(posteriordraws_logbasehaz[,j] + newdat[, j, ] %*% posteriordraws_survcoef[j, ])
      }

      #stop("Not yet implemented for CV association")
    }

    summary_Survsub_list[[idsel]] <- apply(haz_subj, 2, function(x){-cumsum(qwts*x)}) %>%
      exp() %>%
      apply(., 1, function(x) c(mean(x),quantile(x, probs = c(0.025,0.975)))) %>%
      t() %>%
      data.frame("time" = qpts,
                 . ,
                 "outcome" = 3,
                 "observed" = NA,
                 "id" = id_selected[idsel],
                 "time_vline" = last_Longtime[idsel]) %>%
      rbind(c(last_Longtime[idsel], 1, 1, 1, 3, NA, id_selected[idsel], last_Longtime[idsel]),
            .)

  }

  summary_Survsub <- data.table::rbindlist(summary_Survsub_list) %>%
    rbind(list(last_Longtime[idsel], NA, 0, 0, 3, NA, id_selected[idsel], last_Longtime[idsel]))
  ##to show zero in y-axis of survival
  colnames(summary_Survsub) <- c("time", "Estimate", "Q2.5", "Q97.5", "outcome", "observed", "id", "time_vline")


  data_plotjm <- rbind(data_pred, summary_Survsub)
  outcome_labels <- if(is.null(var_labels)){
    c(as.character(formula.tools::lhs(standata$formulaLong1[[1]])),
      as.character(formula.tools::lhs(standata$formulaLong2[[1]])),
      "Event")
  }else{var_labels}
  names(outcome_labels) <- c("1", "2", "3")


  jm_plot <- data_plotjm %>%
    ggplot(data = ., aes(x = time, y = observed)) +
    facet_grid(outcome ~ id, scales="free_y", labeller = labeller(outcome = outcome_labels)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5), alpha = 0.2, fill = "red") +
    geom_vline(aes(xintercept = time_vline), linetype = 1, col = "black") +
    geom_line(aes(x = time, y = Estimate), linetype = 2, col = "red") +
    labs(x = standata$time_varLong, y = NULL) +
    lims(x = xlim_input)

  #suppressWarnings(print(jm_plot))
  if(plot_jm == TRUE){
    suppressWarnings(print(jm_plot))
  }else{
    return(jm_plot)
  }
}
