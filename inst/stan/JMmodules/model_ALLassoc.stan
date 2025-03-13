#include pre.stan


//---- Log-lik for longitudinal submodels

{
  // declare linear predictors
  vector[y_N[1]] y1mu_eta; 
  vector[y_N[1]] y1sigma_eta;
  vector[y_N[2]] y2mu_eta;
  vector[y_N[2]] y2sigma_eta;
  
  // evaluate linear predictor (i.e. mu and log(sigma)) for each long. submodel - note that the shift is cumulative, so that you consider the correct part of b-mat for each one 
  y1mu_eta = evaluate_eta(y1mu_X, y1mu_Z, y1_Z_id, y1mu_gamma, y1mu_beta, b_mat, 0);
  y2mu_eta = evaluate_eta(y2mu_X, y2mu_Z, y2_Z_id, y2mu_gamma, y2mu_beta, b_mat, b_KM[1] + b_KM[2]);   // bmu_K[1] + bsigma_K[1]
  y1sigma_eta = evaluate_eta(y1sigma_X, y1sigma_Z, y1_Z_id, y1sigma_gamma, y1sigma_beta, b_mat, b_KM[1]);     // bmu_K[1]
  y2sigma_eta = evaluate_eta(y2sigma_X, y2sigma_Z, y2_Z_id, y2sigma_gamma, y2sigma_beta, b_mat, b_KM[1] + b_KM[2] + b_KM[3]);  // bmu_K[1] + bsigma_K[1] + bmu_K[2]
  
  // increment the target with the log-lik
  target += normal_lpdf(y1 | y1mu_eta, exp(y1sigma_eta));
  target += normal_lpdf(y2 | y2mu_eta, exp(y2sigma_eta));
}

//----- Log-lik for event submodel (Gauss-Kronrod quadrature)

{

  vector[nrow_e_Xq] e_eta_q; 
  vector[nrow_e_Xq] log_basehaz;  // log baseline hazard AT event time and quadrature points
  vector[nrow_e_Xq] log_haz_q;    // log hazard AT event time and quadrature points
  vector[Nevents] log_haz_etimes; // log hazard AT the event time only
  vector[Npat_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points
  
  // Event submodel: linear predictor at event time and quadrature points
  e_eta_q = e_Xq * e_beta;
  
  
  if (assoc_code == 0) //RE association
    e_eta_q += a_beta[1] * b_mat[y1_Zq_id, 1] + a_beta[2] * b_mat[y1_Zq_id, 2] + a_beta[3] * b_mat[y2_Zq_id, 3] + a_beta[4] * b_mat[y2_Zq_id, 4];   
  else { //LP or CV association
    
    // Long. submodel: linear predictor at event time and quadrature points
    vector[nrow_y_Xq[1]] y1mu_eta_q; 
    vector[nrow_y_Xq[2]] y2mu_eta_q;
    vector[nrow_y_Xq[1]] y1sigma_eta_q; 
    vector[nrow_y_Xq[2]] y2sigma_eta_q;
    
    y1mu_eta_q = evaluate_eta(y1mu_Xq, y1mu_Zq, y1_Zq_id, y1mu_gamma, y1mu_beta, b_mat, 0);
    y2mu_eta_q = evaluate_eta(y2mu_Xq, y2mu_Zq, y2_Zq_id, y2mu_gamma, y2mu_beta, b_mat, b_KM[1] + b_KM[2]);
    y1sigma_eta_q = evaluate_eta(y1sigma_Xq, y1sigma_Zq, y1_Zq_id, y1sigma_gamma, y1sigma_beta, b_mat, b_KM[1]);
    y2sigma_eta_q = evaluate_eta(y2sigma_Xq, y2sigma_Zq, y2_Zq_id, y2sigma_gamma, y2sigma_beta, b_mat, b_KM[1] + b_KM[2] + b_KM[3]);
    
    if (assoc_code == 2) { //CV association
      y1sigma_eta_q = exp(y1sigma_eta_q);
      y2sigma_eta_q = exp(y2sigma_eta_q);
    }
    
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    e_eta_q += a_beta[1] * y1mu_eta_q + a_beta[2] * y1sigma_eta_q + a_beta[3] * y2mu_eta_q + a_beta[4] * y2sigma_eta_q;
    
  }  
  
  
  
  // Log baseline hazard at event time and quadrature points
  log_basehaz = norm_const + basehaz_X * e_aux; //added with respect to Brilleman (it appears in rstanarm, https://github.com/stan-dev/rstanarm/blob/master/src/stan_files/model/event_lp.stan)
  //log_basehaz = basehaz_X * e_aux;
  
  // Log hazard at event time and quadrature points
  log_haz_q = log_basehaz + e_eta_q;
  
  // Log hazard at event times only
  log_haz_etimes = head(log_haz_q, Nevents);
  
  // Log hazard at quadrature points only
  log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);
  
  // Log likelihood for event submodel
  // NB The first term is the log hazard contribution to the log  
  // likelihood for the event submodel. The second term is the log  
  // survival contribution to the log likelihood for the event submodel.  
  // The latter is obtained by summing over the quadrature points to get 
  // the approximate integral (i.e. cumulative hazard). Note that the
  // 'qwts' vector already incorporates (b-a)/2 scaling such that the
  // integral is evaluated over limits (a,b) rather than (-1,+1), where
  // 'a' is baseline, i.e. time 0, and 'b' is the event or censoring
  // time for the individual.
  target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));  
}    

//----- Log-priors

// intercepts for long. submodels
target += normal_lpdf(y1mu_gamma | 
ymu_prior_mean_for_intercept[1], ymu_prior_scale_for_intercept[1]);    
target += normal_lpdf(y2mu_gamma | 
ymu_prior_mean_for_intercept[2], ymu_prior_scale_for_intercept[2]); 
target += normal_lpdf(y1sigma_gamma | 
ysigma_prior_mean_for_intercept[1], ysigma_prior_scale_for_intercept[1]); 
target += normal_lpdf(y2sigma_gamma | 
ysigma_prior_mean_for_intercept[2], ysigma_prior_scale_for_intercept[2]);    

// coefficients for long. submodels   
target += normal_lpdf(y1mu_z_beta | 0, 1);
target += normal_lpdf(y1sigma_z_beta | 0, 1);
target += normal_lpdf(y2mu_z_beta | 0, 1);
target += normal_lpdf(y2sigma_z_beta | 0, 1);

// coefficients for event submodel
target += normal_lpdf(e_z_beta | 0, 1);
target += normal_lpdf(a_z_beta | 0, 1);

// residual error SDs for long. submodels- not applicable to MELSM
//target += normal_lpdf(y1_aux_unscaled | 0, 1);
//target += normal_lpdf(y2_aux_unscaled | 0, 1);

// b-spline coefs for baseline hazard
target += normal_lpdf(e_aux_unscaled | 0, 1);

// group level terms
// sds
target += student_t_lpdf(b_sd | b_prior_df, 0, b_prior_scale); //half-student T because b_sd is defined as <lower = 0>
//target += gamma_lpdf(b_sd | 2, 0.5);  //parameters are shape and inverse scale (i.e., rate)
// primitive coefs
target += normal_lpdf(to_vector(z_b_mat) | 0, 1); 
// corr matrix
if (b_K > 1) 
  target += lkj_corr_cholesky_lpdf(b_cholesky | b_prior_regularization);

