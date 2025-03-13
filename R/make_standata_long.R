#' Create Stan data input for the longitudinal (MELS) submodel
#'
#' @param formulaLong1 Formula for the first longitudinal marker.
#' @param formulaLong2 Formula for the second longitudinal marker.
#' @param dataLong Longitudinal dataset.
#' @param id_var Character string for name of ID variable in `dataLong`.
#' @param time_varLong Character string for name of time variable in `dataLong`.
#'
#' @import dplyr
#' @export
#' @return A list to be passed to `make_standata_jmwiv.R`.
#'
#' @examples
#' library(rstanarm)
#' data(pbcLong)
#' f_logBili <- list(formula(logBili ~ sex + trt + year + (1 | id)),
#' formula(sigma ~ sex + trt + year + (1 | id)))
#' f_albumin <- list(formula(albumin ~ sex + trt + year + (1 | id)),
#' formula(sigma ~ sex + trt + year + (1 | id)))
#' u1 <- make_standata_long(formulaLong1 = f_logBili,
#'                          formulaLong2 = f_albumin,
#'                          dataLong = pbcLong,
#'                          id_var = "id",
#'                          time_varLong = "year")

make_standata_long <- function(formulaLong1,
                               formulaLong2 = NULL,
                               dataLong,
                               id_var,
                               time_varLong
){

  try(if(is.null(formulaLong1)) stop("You need to provide at least a formula!"))

  M <- ifelse(is.null(formulaLong2), 1L, 2L)
  ymu_K <- rep(NA, M)  %>% as.array()
  ysigma_K <- rep(NA, M)  %>% as.array()
  bmu_K <- rep(NA, M)  %>% as.array()
  bsigma_K <- rep(NA, M)  %>% as.array()
  y_N <- rep(nrow(dataLong), M)

  b_N <- dataLong %>%
    select(all_of(id_var)) %>%
    pull() %>%
    as.factor() %>%
    nlevels()

  f1mu_dat <- model.matrix(formulaLong1[[1]][-2] %>% lme4::nobars(),
                           data = dataLong)
  resp1 <-  as.character(terms(formulaLong1[[1]])[[2]])
  y1 <- dataLong %>%
    select(all_of(resp1)) %>%
    pull()
  y1mu_XRAW <- data.frame(f1mu_dat, check.names = FALSE) %>%
    select(., !contains(c("|", "(Intercept)", "s(")))
  y1mu_Xbar <- colMeans(y1mu_XRAW) %>% as.array()
  y1mu_X <- y1mu_XRAW %>%
    mutate_all(scale, scale = FALSE) %>%
    as.matrix()
  ymu_K[1] <- length(y1mu_Xbar)

  y1mu_ranef <- formulaLong1[[1]][[3]] %>%   #lme4::findbars(formulaLong1[[1]])[[1]]
    as.character() %>%
    stringr::str_subset(pattern = id_var) %>%
    stringr::str_match_all(., "(?<=\\().+?(?=\\))") %>%
    simplify2array()

  y1mu_Z <- sub("\\|.*", "", y1mu_ranef) %>%
    paste("~", .) %>%
    formula() %>%
    model.matrix(., data = dataLong) %>%
    t()

  y1_Z_id <- dataLong %>%
    select(all_of(id_var)) %>%
    pull()

  bmu_K[1] <- nrow(y1mu_Z)


  f1sigma_dat <- model.matrix(formulaLong1[[2]][-2] %>% lme4::nobars(),
                              data = dataLong)  ###you need to remove sigma which is not in the dataset
  y1sigma_XRAW <- data.frame(f1sigma_dat, check.names = FALSE) %>%
    select(., !contains(c("|", "(Intercept)", "s(")))
  y1sigma_Xbar <- colMeans(y1sigma_XRAW) %>% as.array()
  y1sigma_X <- y1sigma_XRAW %>%
    mutate_all(scale, scale = FALSE) %>%
    as.matrix()
  ysigma_K[1] <- length(y1sigma_Xbar)


  y1sigma_ranef <- formulaLong1[[2]][[3]] %>%
    as.character() %>%
    stringr::str_subset(pattern = id_var) %>%
    stringr::str_match_all(., "(?<=\\().+?(?=\\))") %>%
    simplify2array()

  y1sigma_Z <- sub("\\|.*", "", y1sigma_ranef) %>%
    paste("~", .) %>%
    formula() %>%
    model.matrix(., data = dataLong) %>%
    t()

  bsigma_K[1] <- nrow(y1sigma_Z)

  if(!is.null(formulaLong2)){
    f2_dat <- model.matrix(formulaLong2[[1]][-2] %>% lme4::nobars(),
                           data = dataLong)
    resp2 <-  as.character(terms(formulaLong2[[1]])[[2]])
    y2 <- dataLong %>%
      select(all_of(resp2)) %>%
      pull()
    y2mu_XRAW <- data.frame(f2_dat, check.names = FALSE) %>%
      select(., !contains(c("|", "(Intercept)", "s(")))
    y2mu_Xbar <- colMeans(y2mu_XRAW) %>% as.array()
    y2mu_X <- y2mu_XRAW %>%
      mutate_all(scale, scale = FALSE) %>%
      as.matrix()
    ymu_K[2] <- length(y2mu_Xbar)

    y2mu_Z <- formulaLong2[[1]][[3]] %>%
      as.character() %>%
      stringr::str_subset(pattern = id_var) %>%
      stringr::str_match_all(., "(?<=\\().+?(?=\\))") %>%
      simplify2array() %>%
      sub("\\|.*", "", .) %>%
      paste("~", .) %>%
      formula() %>%
      model.matrix(., data = dataLong) %>%
      t()

    y2_Z_id <- dataLong %>%
      select(all_of(id_var)) %>%
      pull()

    bmu_K[2] <- nrow(y2mu_Z)


    f2_datsigma <- model.matrix(formulaLong2[[2]][-2] %>% lme4::nobars(),
                                data = dataLong)  ###you need to remove sigma which is not in the dataset
    y2sigma_XRAW <- data.frame(f2_datsigma, check.names = FALSE) %>%
      select(., !contains(c("|", "(Intercept)", "s(")))
    y2sigma_Xbar <- colMeans(y2sigma_XRAW) %>% as.array()
    y2sigma_X <- y2sigma_XRAW %>%
      mutate_all(scale, scale = FALSE) %>%
      as.matrix()
    ysigma_K[2] <- length(y2sigma_Xbar)



    y2sigma_Z <- formulaLong2[[2]][[3]] %>%
      as.character() %>%
      stringr::str_subset(pattern = id_var) %>%
      stringr::str_match_all(., "(?<=\\().+?(?=\\))") %>%
      simplify2array() %>%
      sub("\\|.*", "", .) %>%
      paste("~", .) %>%
      formula() %>%
      model.matrix(., data = dataLong) %>%
      t()

    bsigma_K[2] <- nrow(y2sigma_Z)

    output2 <- nlist(y2, y2mu_X, y2sigma_X, y2mu_Xbar, y2sigma_Xbar,
                     y2mu_Z, y2sigma_Z, y2_Z_id, formulaLong2)
  }

  b_K <- sum(bmu_K, bsigma_K, na.rm = T)

  output1 <- nlist(M, time_varLong, id_var,
                   y_N, ymu_K, ysigma_K,
                   y1, y1mu_X, y1sigma_X, y1mu_Xbar, y1sigma_Xbar,
                   b_N, b_K, bmu_K, bsigma_K,
                   y1mu_Z, y1sigma_Z, y1_Z_id,
                   formulaLong1)


  if(M == 1) output1
  else c(output1, output2)
}
