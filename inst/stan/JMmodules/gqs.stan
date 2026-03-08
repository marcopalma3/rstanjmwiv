#include pre.stan


  

  
  real y1mu_Intercept; // transformed intercepts for long. submodels
  real y2mu_Intercept;
  real y1sigma_Intercept; // transformed intercepts for long. submodels
  real y2sigma_Intercept;
  real e_Intercept; // transformed intercept for event submodel 
  vector[size(b_cov_idx)] b_cov; // var-cov for REs
  
 
  
  vector[Npat] log_lik;
  
 
  vector[y_N[1]] y1mu_eta; 
  vector[y_N[1]] y1sigma_eta;
  vector[y_N[2]] y2mu_eta;
  vector[y_N[2]] y2sigma_eta;
  
  
  vector[nrow_e_Xq] e_eta_q; 
  vector[nrow_e_Xq] log_basehaz;
  vector[nrow_e_Xq] log_haz_q;
  vector[Nevents] log_haz_etimes;
  vector[Npat_times_qnodes] log_haz_qtimes;
 
  vector[nrow_y_Xq[1]] y1mu_eta_q; 
  vector[nrow_y_Xq[2]] y2mu_eta_q;
  vector[nrow_y_Xq[1]] y1sigma_eta_q; 
  vector[nrow_y_Xq[2]] y2sigma_eta_q;
  
  
// Transformed intercepts for long. submodels
y1mu_Intercept = y1mu_gamma - dot_product(y1mu_Xbar, y1mu_beta);
y2mu_Intercept = y2mu_gamma - dot_product(y2mu_Xbar, y2mu_beta);
y1sigma_Intercept = y1sigma_gamma - dot_product(y1sigma_Xbar, y1sigma_beta);
y2sigma_Intercept = y2sigma_gamma - dot_product(y2sigma_Xbar, y2sigma_beta);


// Transformed intercept for event submodel 
e_Intercept = 0.0 - dot_product(e_Xbar, e_beta);  //do I need to add norm_const?

// Transform variance-covariance matrix for REs
if (b_K == 1)
  b_cov[1] = b_sd[1] * b_sd[1];
else
  b_cov = to_vector(quad_form_diag(multiply_lower_tri_self_transpose(b_cholesky), 
  b_sd))[b_cov_idx];
  
  
 
  
  y1mu_eta = evaluate_eta(y1mu_X, y1mu_Z, y1_Z_id, y1mu_gamma, y1mu_beta, b_mat, 0);
  y2mu_eta = evaluate_eta(y2mu_X, y2mu_Z, y2_Z_id, y2mu_gamma, y2mu_beta, b_mat, b_KM[1] + b_KM[2]);
  y1sigma_eta = evaluate_eta(y1sigma_X, y1sigma_Z, y1_Z_id, y1sigma_gamma, y1sigma_beta, b_mat, b_KM[1]);
  y2sigma_eta = evaluate_eta(y2sigma_X, y2sigma_Z, y2_Z_id, y2sigma_gamma, y2sigma_beta, b_mat, b_KM[1] + b_KM[2] + b_KM[3]);
  
  
  
  

  e_eta_q = e_Xq * e_beta;
  

  y1mu_eta_q = evaluate_eta(y1mu_Xq, y1mu_Zq, y1_Zq_id, y1mu_gamma, y1mu_beta, b_mat, 0);
  y2mu_eta_q = evaluate_eta(y2mu_Xq, y2mu_Zq, y2_Zq_id, y2mu_gamma, y2mu_beta, b_mat, b_KM[1] + b_KM[2]);
  y1sigma_eta_q = evaluate_eta(y1sigma_Xq, y1sigma_Zq, y1_Zq_id, y1sigma_gamma, y1sigma_beta, b_mat, b_KM[1]);
  y2sigma_eta_q = evaluate_eta(y2sigma_Xq, y2sigma_Zq, y2_Zq_id, y2sigma_gamma, y2sigma_beta, b_mat, b_KM[1] + b_KM[2] + b_KM[3]);




// Event submodel contrib (LP / CV)
if (assoc_code == 2) {
    
    y1sigma_eta_q = exp(y1sigma_eta_q);
    y2sigma_eta_q = exp(y2sigma_eta_q);
}
e_eta_q += a_beta[1] * y1mu_eta_q + 
           a_beta[2] * y1sigma_eta_q + 
           a_beta[3] * y2mu_eta_q + 
           a_beta[4] * y2sigma_eta_q;
  
  
  log_basehaz = norm_const + basehaz_X * e_aux;
  
 
  log_haz_q = log_basehaz + e_eta_q;
  

  log_haz_etimes = head(log_haz_q, Nevents);
  log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);
  
  
  //log-likelihood per patient
  
  for (i in 1:Npat) {
    log_lik[i] = 0;
    
    // Long contrib of biomarker 1  : Find all the observations of patient i and add their log-likelihood ( could be optimized here instead of looping over all y_N[1] measurements) 
    for (j in 1:y_N[1]) {
      if (y1_Z_id[j] == i) {
        log_lik[i] += normal_lpdf(y1[j] | y1mu_eta[j], exp(y1sigma_eta[j]));
      }
    }
    
    // Long contrib of  biomarker 2
    for (j in 1:y_N[2]) {
      if (y2_Z_id[j] == i) {
        log_lik[i] += normal_lpdf(y2[j] | y2mu_eta[j], exp(y2sigma_eta[j]));
      }
    }
    
    // Event contrib (if event)
    int event_idx = 0;
    for (k in 1:Nevents) {
      if (y1_Zq_id[k] == i) {
        event_idx = k;
        break;
      }
    }
    if (event_idx > 0) {
      log_lik[i] += log_haz_etimes[event_idx];
    }
    
    // Cumulative hazard contrib
    real cumhaz = 0;
    for (q in 1:qnodes) {
      int quad_idx = (i-1)*qnodes + q;
      cumhaz += qwts[quad_idx] * exp(log_haz_qtimes[quad_idx]);
    }
    log_lik[i] -= cumhaz;
  }