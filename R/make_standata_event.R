#' Create Stan data input for the event submodel
#'
#' @param formulaEvent Formula for the time-to-event submodel (without association terms)
#' @param dataEvent Event dataset
#' @param M Number of longitudinal biomarkers (at the moment, only 2)
#' @param time_varLong Character string for name of time variable in `dataLong`.
#' @param standata_long An object returned by `make_standata_long.R`
#' @param id_var Character string for name of ID variable in `dataEvent`.
#' @param a_K Number of association terms (use 4: mean and WIV for each biomarker)
#' @param assoc Type of association (either 'RE' for random effects, 'LP' for linear predictor, 'CV' for current value).
#' @param basehaz Type of baseline hazard function (only 'bs' is implemented)
#' @param basehaz_aux List of parameters for the baseline hazard function.
#' @param qnodes Number of quadrature nodes.
#'
#' @import dplyr
#'
#' @return A list to be passed to `make_standata_jmwiv.R`.
make_standata_event <- function(formulaEvent,
                                dataEvent,
                                id_var,
                                standata_long,
                                assoc,
                                a_K = 4,
                                basehaz = "bs",
                                basehaz_aux = list(df = 6,
                                                   knots = NULL,
                                                   degree = 3),
                                qnodes = 15L){



  if(assoc %in% c("RE", "LP", "CV")){
    assoc_code <- case_when(
      assoc == "RE" ~ 0,
      assoc == "LP" ~ 1,
      assoc == "CV" ~ 2,
    )
  }else stop("Association type is not valid. Association must be equal to either 'RE', 'LP' or 'CV'.")

  qnodes <- as.integer(qnodes)
  M <- standata_long$M
  time_varLong <- standata_long$time_varLong

  y1mu_Xbar <- standata_long$y1mu_Xbar
  y1sigma_Xbar <- standata_long$y1sigma_Xbar
  y2mu_Xbar <- standata_long$y2mu_Xbar
  y2sigma_Xbar <- standata_long$y2sigma_Xbar


  time_var <- as.character(terms(formulaEvent)[[2]])[2]
  event_indicator <- as.character(terms(formulaEvent)[[2]])[3]



  e_modelmat <- model.matrix(formulaEvent, data = dataEvent)[,-1]  %>%
    as.matrix()
  e_Xbar <- colMeans(e_modelmat)
  e_K <- ncol(e_modelmat)

  Npat <- pull(dataEvent, all_of(id_var)) %>%
    as.factor() %>%
    droplevels() %>%
    nlevels

  data_onlyevent <- dataEvent %>%
    filter(., get(event_indicator) != 0) %>%
    data.frame(row.names = pull(., all_of(id_var)))

  Nevents <- nrow(data_onlyevent)

  Npat_times_qnodes <- Npat*qnodes
  nrow_e_Xq <- Nevents + Npat_times_qnodes %>%
    as.integer

  quad <- mvQuad::QuadRules[["GKr"]][[qnodes]]
  qwts <- quad$w %x% pull(dataEvent, all_of(time_var)) %>%
    setNames(., rep(as.character(1:Npat), qnodes)) %>%
    as.array()
  qpts <- quad$n %x% pull(dataEvent, all_of(time_var)) %>%
    as.numeric()




  replicatedData <- dataEvent %>%
    slice(rep(row_number(), qnodes)) %>%  ###changed from 15
    bind_rows(data_onlyevent, .) %>%
    mutate("(Intercept)" = 1,
           !!time_varLong := c(pull(data_onlyevent, all_of(time_var)), qpts))


  e_times <- replicatedData %>%
    pull(!!time_varLong) %>%
    setNames(., pull(replicatedData, all_of(id_var))) %>%
    as.array()

  e_XqRAW <- model.matrix(formulaEvent, data = replicatedData)[,-1] %>%
    as.matrix()
  e_Xbar <- colMeans(e_XqRAW) %>% as.array()
  e_Xq <- e_XqRAW %>%
    scale(., scale = FALSE) %>%
    as.array()

  norm_const <- log(Nevents / sum(pull(data_onlyevent, all_of(time_var))))  ##not used in the model


  ###basehaz modified from https://github.com/stan-dev/rstanarm/blob/master/R/jm_data_block.R, line 1211
  if(basehaz == "bs"){
    basehaz_df <- as.integer(basehaz_aux[["df"]])
    basehaz_X <- splines::bs(pull(data_onlyevent, all_of(time_var)),
                            df = basehaz_aux$df,
                            knots = basehaz_aux$knots,
                            degree = basehaz_aux$degree,
                            Boundary.knots = c(0, max(pull(dataEvent, all_of(time_var)))),
                            intercept = TRUE) %>%
      predict(., e_times) %>%
      as.array()
  }




  nrow_y_Xq <- rep(nrow_e_Xq, M)

  y1mu_Xq <- if(length(names(y1mu_Xbar)) == 1 && names(y1mu_Xbar) == time_varLong){
    magrittr::subtract(e_times, standata_long$y1mu_Xbar[[time_varLong]]) %>%
      as.matrix()
  }else{
      with(standata_long, data.frame(y1_Z_id, y1mu_X)) %>%
      rename(!!id_var := 1) %>%
      select(all_of(id_var), names(y1mu_Xbar), -all_of(time_varLong)) %>%
      unique() %>%
      left_join(
        sweep(replicatedData[,c(id_var, time_varLong)],
              2,
              c(0, standata_long$y1mu_Xbar[[time_varLong]])
        ),
        ., by = id_var) %>%
      select(names(y1mu_Xbar)) %>%
      as.matrix()
  }

  y1sigma_Xq <- if(length(names(y1sigma_Xbar)) == 1 && names(y1sigma_Xbar) == time_varLong){
    magrittr::subtract(e_times, standata_long$y1sigma_Xbar[[time_varLong]]) %>%
      as.matrix()
  }else{
    with(standata_long, data.frame(y1_Z_id, y1sigma_X)) %>%
    rename(!!id_var := 1) %>%
    select(all_of(id_var), names(y1sigma_Xbar), -all_of(time_varLong)) %>%
    unique() %>% # it works only for baseline covariates in the longitudinal model
    left_join(
      sweep(replicatedData[,c(id_var, time_varLong)],
            2,
            c(0, standata_long$y1sigma_Xbar[[time_varLong]])
      ),
      ., by = id_var) %>%
    select(names(y1sigma_Xbar)) %>%
    as.matrix()   ##remove ID and reorder columns
  }



  y1mu_Zq <- replicatedData %>%
    select("(Intercept)") %>%
    as.matrix() %>%
    t()
  y1sigma_Zq <- replicatedData %>%
    select("(Intercept)") %>%
    as.matrix() %>%
    t()
  y1_Zq_id <- pull(replicatedData, id_var) %>% as.integer()

  if(M == 2){

    y2mu_Xq <- if(length(names(y2mu_Xbar)) == 1 && names(y2mu_Xbar) == time_varLong){
      magrittr::subtract(e_times, standata_long$y2mu_Xbar[[time_varLong]]) %>%
        as.matrix()
    }else{
      with(standata_long, data.frame(y2_Z_id, y2mu_X)) %>%
      rename(!!id_var := 1) %>%       #give ID name to first column
      select(all_of(id_var), names(y2mu_Xbar), -all_of(time_varLong)) %>% #you have already the quadrature times as timevar_Long
      unique() %>% # it works only for baseline covariates in the longitudinal model
      left_join(
        sweep(replicatedData[,c(id_var, time_varLong)],
              2,
              c(0, standata_long$y2mu_Xbar[[time_varLong]])
        ),
        ., by = id_var) %>%
      select(names(y2mu_Xbar)) %>%
      as.matrix()
    }



    y2sigma_Xq <- if(length(names(y2sigma_Xbar)) == 1 && names(y2sigma_Xbar) == time_varLong){
      magrittr::subtract(e_times, standata_long$y2sigma_Xbar[[time_varLong]]) %>%
        as.matrix()
    }else{
      with(standata_long, data.frame(y2_Z_id, y2sigma_X)) %>%
        rename(!!id_var := 1) %>%       #give ID name to first column
        select(all_of(id_var), names(y2sigma_Xbar), -all_of(time_varLong)) %>% #you have already the quadrature times as timevar_Long
        unique() %>% # it works only for baseline covariates in the longitudinal model
        left_join(
          sweep(replicatedData[,c(id_var, time_varLong)],
                2,
                c(0, standata_long$y2sigma_Xbar[[time_varLong]])
          ),
          ., by = id_var) %>%
        select(names(y2sigma_Xbar))%>%
        as.matrix()   ##remove ID and reorder columns
    }





    y2mu_Zq <- replicatedData %>%
      select("(Intercept)") %>%
      as.matrix() %>%
      t()
    y2sigma_Zq <- replicatedData %>%
      select("(Intercept)") %>%
      as.matrix() %>%
      t()
    y2_Zq_id <- pull(replicatedData, id_var) %>% as.integer()

    output2 <- rstan::nlist(y2mu_Xq, y2sigma_Xq, y2mu_Zq, y2sigma_Zq, y2_Zq_id)
  }



  output1 <- rstan::nlist(e_K, Npat, Nevents, qnodes, Npat_times_qnodes, nrow_e_Xq,
                   e_times, e_Xq, e_Xbar, basehaz_df, basehaz_X, qwts,
                   a_K, norm_const,
                   assoc_code, e_modelmat,
                   nrow_y_Xq, y1mu_Xq, y1sigma_Xq,
                   y1mu_Zq, y1sigma_Zq, y1_Zq_id, formulaEvent)


  if(M == 1) output1
  else c(output1, output2)
}
