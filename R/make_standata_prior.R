#' Create Stan data input for the prior distributions
#'
#' @param standata_long An output from `make_standata_long.R`.
#' @param standata_event An output from `make_standata_event.R`.
#' @param prior_input A list with prior parameters.
#'
#' @import dplyr rstan
#' @export
#' @return A list with all prior parameters to be passed to `make_standata_jmwiv.R`.


make_standata_prior <- function(standata_long,
                                standata_event,
                                prior_input = NULL
                                ){

  prior_list <- rstan::nlist(y1mu_prior_mean = NULL,
                      y2mu_prior_mean = NULL,
                      y1sigma_prior_mean= NULL,
                      y2sigma_prior_mean= NULL,
                      e_prior_mean = NULL,
                      a_prior_mean = NULL,
                      ymu_prior_mean_for_intercept = NULL,
                      ysigma_prior_mean_for_intercept = NULL,
                      e_prior_mean_for_aux = NULL,
                      y1mu_prior_scale = NULL,
                      y2mu_prior_scale = NULL,
                      y1sigma_prior_scale = NULL,
                      y2sigma_prior_scale = NULL,
                      e_prior_scale = NULL,
                      a_prior_scale = NULL,
                      ymu_prior_scale_for_intercept = NULL,
                      ysigma_prior_scale_for_intercept = NULL,
                      e_prior_scale_for_aux = NULL,
                      b_prior_scale = NULL,
                      b_prior_df = NULL,
                      b_prior_regularization = NULL)

  prior_input <- sapply(prior_input, as.array) ###need to use as.array to store dimensions as well

  if(!is.null(prior_input)){
    if(sum(names(prior_input) %in% names(prior_list)) != length(names(prior_input)))
      stop("One or more input names are not valid.\n
         Please provide prior values only for existing parameters.")

    for(i in 1:length(names(prior_input)))
      prior_list[[names(prior_input)[i]]] <- prior_input[[names(prior_input)[i]]]
  }





    list2env(prior_list, env = environment())



  if(!is.null(y1mu_prior_mean) & length(y1mu_prior_mean) != standata_long$ymu_K[1])
    stop("Please check the dimensions of y1mu_prior_mean!")
  if(!is.null(y1sigma_prior_mean) & length(y1sigma_prior_mean) != standata_long$ysigma_K[1])
    stop("Please check the dimensions of y1sigma_prior_mean!")
  if(!is.null(e_prior_mean) & length(e_prior_mean) != standata_event$e_K)
    stop("Please check the dimensions of e_prior_mean!")
  if(!is.null(a_prior_mean) & length(a_prior_mean) != standata_event$a_K)
    stop("Please check the dimensions of a_prior_mean!")
  if(!is.null(ymu_prior_mean_for_intercept) & length(ymu_prior_mean_for_intercept) != standata_long$M)
    stop("Please check the dimensions of ymu_prior_mean_for_intercept!")
  if(!is.null(ysigma_prior_mean_for_intercept) & length(ysigma_prior_mean_for_intercept) != standata_long$M)
    stop("Please check the dimensions of ysigma_prior_mean_for_intercept!")
  if(!is.null(e_prior_mean_for_aux) & length(e_prior_mean_for_aux) != standata_event$basehaz_df)
    stop("Please check the dimensions of e_prior_mean_for_aux!")


  if(!is.null(y1mu_prior_scale)){
    if(length(y1mu_prior_scale) != standata_long$ymu_K[1])
      stop("Please check the dimensions of y1mu_prior_scale!")
    else if(y1mu_prior_scale < 0)
      stop("Please provide a non-negative value for y1mu_prior_scale!")
  }

  if(!is.null(y1sigma_prior_scale)){
    if(length(y1sigma_prior_scale) != standata_long$ysigma_K[1])
      stop("Please check the dimensions of y1sigma_prior_scale!")
    else if(y1sigma_prior_scale< 0)
      stop("Please provide a non-negative value for y1sigma_prior_scale!")
  }

  if(!is.null(e_prior_scale)){
    if(length(e_prior_scale) != standata_event$e_K)
      stop("Please check the dimensions of e_prior_scale!")
    else if(sum(e_prior_scale< 0) != 0)
      stop("Please provide a non-negative value for e_prior_scale!")
  }



  if(!is.null(a_prior_scale)){
    if(length(a_prior_scale) != standata_event$a_K)
      stop("Please check the dimensions of a_prior_scale!")
    else if(sum(a_prior_scale< 0) != 0)
      stop("Please provide a non-negative value for a_prior_scale!")
  }


  if(!is.null(ymu_prior_scale_for_intercept)){
    if(length(ymu_prior_scale_for_intercept) != standata_long$M)
      stop("Please check the dimensions of ymu_prior_scale_for_intercept!")
    else if(sum(ymu_prior_scale_for_intercept< 0) != 0)
      stop("Please provide a non-negative value for ymu_prior_scale_for_intercept!")
  }

  if(!is.null(ysigma_prior_scale_for_intercept)){
    if(length(ysigma_prior_scale_for_intercept) != standata_long$M)
      stop("Please check the dimensions of ysigma_prior_scale_for_intercept!")
    else if(sum(ysigma_prior_scale_for_intercept< 0) != 0)
      stop("Please provide a non-negative value for ysigma_prior_scale_for_intercept!")
  }

  if(!is.null(e_prior_scale_for_aux)){
    if(length(e_prior_scale_for_aux) != standata_event$basehaz_df)
      stop("Please check the dimensions of e_prior_scale_for_aux!")
    else if(sum(e_prior_scale_for_aux< 0) != 0)
      stop("Please provide a non-negative value for e_prior_scale_for_aux!")
  }






  if(!is.null(y2mu_prior_mean) & length(y2mu_prior_mean) != standata_long$ymu_K[2])
    stop("Please check the dimensions of y2mu_prior_mean!")
  if(!is.null(y2sigma_prior_mean) & length(y2sigma_prior_mean) != standata_long$ysigma_K[2])
    stop("Please check the dimensions of y2sigma_prior_mean!")

  if(!is.null(y2mu_prior_scale)){
    if(length(y2mu_prior_scale) != standata_long$ymu_K[2])
      stop("Please check the dimensions of y2mu_prior_scale!")
    else if(y2mu_prior_scale < 0)
      stop("Please provide a non-negative value for y2mu_prior_scale!")
  }

  if(!is.null(y2sigma_prior_scale)){
    if(length(y2sigma_prior_scale) != standata_long$ysigma_K[2])
      stop("Please check the dimensions of y2sigma_prior_scale!")
    else if(y2sigma_prior_scale < 0)
      stop("Please provide a non-negative value for y2sigma_prior_scale!")
  }






  if(!is.null(b_prior_scale) & length(b_prior_scale) != standata_long$b_K)
    stop("Please check the dimensions of b_prior_scale!")
  if(!is.null(b_prior_df) & length(b_prior_df) != standata_long$b_K)
    stop("Please check the dimensions of b_prior_df!")
  if(!is.null(b_prior_regularization) & length(b_prior_regularization) != 1)
    stop("Please check the dimensions of b_prior_regularization!")







  if(is.null(y1mu_prior_mean)) y1mu_prior_mean <- rep(0, standata_long$ymu_K[1]) %>% as.array()
  if(is.null(y1sigma_prior_mean)) y1sigma_prior_mean <- rep(0, standata_long$ysigma_K[1]) %>% as.array()
  if(is.null(e_prior_mean)) e_prior_mean <- rep(0, standata_event$e_K) %>% as.array()
  if(is.null(a_prior_mean)) a_prior_mean <- rep(0, standata_event$a_K) %>% as.array()
  if(is.null(ymu_prior_mean_for_intercept)) ymu_prior_mean_for_intercept <- rep(0, standata_long$M) %>% as.array()
  if(is.null(ysigma_prior_mean_for_intercept)) ysigma_prior_mean_for_intercept <- rep(0, standata_long$M) %>% as.array()
  if(is.null(e_prior_mean_for_aux)) e_prior_mean_for_aux <- rep(0, standata_event$basehaz_df) %>% as.array()


  if(is.null(y1mu_prior_scale)) y1mu_prior_scale <- rep(1, standata_long$ymu_K[1]) %>% as.array()
  if(is.null(y1sigma_prior_scale)) y1sigma_prior_scale <- rep(1, standata_long$ysigma_K[1]) %>% as.array()
  if(is.null(e_prior_scale)) e_prior_scale <- rep(1, standata_event$e_K) %>% as.array()
  if(is.null(a_prior_scale)) a_prior_scale <- rep(1, standata_event$a_K) %>% as.array()
  if(is.null(ymu_prior_scale_for_intercept)) ymu_prior_scale_for_intercept <- rep(1, standata_long$M) %>% as.array()
  if(is.null(ysigma_prior_scale_for_intercept)) ysigma_prior_scale_for_intercept <- rep(1, standata_long$M) %>% as.array()
  if(is.null(e_prior_scale_for_aux)) e_prior_scale_for_aux <- rep(1, standata_event$basehaz_df) %>% as.array()

  if(standata_long$M == 2){
  if(is.null(y2mu_prior_mean)) y2mu_prior_mean <- rep(0, standata_long$ymu_K[2]) %>% as.array()
  if(is.null(y2sigma_prior_mean)) y2sigma_prior_mean <- rep(0, standata_long$ysigma_K[2]) %>% as.array()
  if(is.null(y2mu_prior_scale)) y2mu_prior_scale <- rep(1, standata_long$ymu_K[2]) %>% as.array()
  if(is.null(y2sigma_prior_scale)) y2sigma_prior_scale <- rep(1, standata_long$ysigma_K[2]) %>% as.array()
  }

  if(is.null(b_prior_scale)) b_prior_scale <- rep(1, standata_long$b_K) %>% as.array()
  if(is.null(b_prior_df)) b_prior_df <- rep(1, standata_long$b_K) %>% as.array()
  #if(is.null(b_prior_regularization)) b_prior_regularization <- 1L
  b_prior_regularization <- ifelse(is.null(b_prior_regularization), 1L, as.numeric(b_prior_regularization))


  rstan::nlist(y1mu_prior_mean,
                      y2mu_prior_mean,
                      y1sigma_prior_mean,
                      y2sigma_prior_mean,
                      e_prior_mean,
                      a_prior_mean,
                      ymu_prior_mean_for_intercept,
                      ysigma_prior_mean_for_intercept,
                      e_prior_mean_for_aux,
                      y1mu_prior_scale,
                      y2mu_prior_scale,
                      y1sigma_prior_scale,
                      y2sigma_prior_scale,
                      e_prior_scale,
                      a_prior_scale,
                      ymu_prior_scale_for_intercept,
                      ysigma_prior_scale_for_intercept,
                      e_prior_scale_for_aux,
                      b_prior_scale,
                      b_prior_df,
                      b_prior_regularization)

}

