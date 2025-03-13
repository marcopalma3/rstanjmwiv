#include pre.stan

vector[ymu_K[1]] y1mu_beta;     // population level params for long. submodels
vector[ysigma_K[1]] y1sigma_beta;     // population level params for long. submodels
vector[ymu_K[2]] y2mu_beta;     // population level params for long. submodels
vector[ysigma_K[2]] y2sigma_beta;
vector[e_K] e_beta;       // coefs in event submodel (log hazard ratios)
vector[a_K] a_beta;       // assoc params in event submodel (log hazard ratios) 
vector[basehaz_df] e_aux; // b-spline coefs for baseline hazard
matrix[b_N,b_K] b_mat;    // group level params

// coefs for long. submodels
y1mu_beta = y1mu_z_beta .* y1mu_prior_scale + y1mu_prior_mean;
y2mu_beta = y2mu_z_beta .* y2mu_prior_scale + y2mu_prior_mean;
y1sigma_beta = y1sigma_z_beta .* y1sigma_prior_scale + y1sigma_prior_mean;
y2sigma_beta = y2sigma_z_beta .* y2sigma_prior_scale + y2sigma_prior_mean;


// coefs for event submodel (incl. association parameters)
e_beta = e_z_beta .* e_prior_scale + e_prior_mean;
a_beta = a_z_beta .* a_prior_scale + a_prior_mean;

// residual error SDs for long. submodels - not applicable to MELSM
//y1_aux = y1_aux_unscaled * y_prior_scale_for_aux[1] + y_prior_mean_for_aux[1];
//y2_aux = y2_aux_unscaled * y_prior_scale_for_aux[2] + y_prior_mean_for_aux[2];

// b-spline coefs for baseline hazard
e_aux = e_aux_unscaled .* e_prior_scale_for_aux + e_prior_mean_for_aux;

// group level params
if (b_K == 1) 
  b_mat = (b_sd[1] * z_b_mat)'; 
else if (b_K > 1) 
  b_mat = (diag_pre_multiply(b_sd, b_cholesky) * z_b_mat)';

