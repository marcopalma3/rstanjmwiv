#include pre.stan


real y1mu_Intercept; // transformed intercepts for long. submodels
real y2mu_Intercept;
real y1sigma_Intercept; // transformed intercepts for long. submodels
real y2sigma_Intercept;
real e_Intercept; // transformed intercept for event submodel 
vector[size(b_cov_idx)] b_cov; // var-cov for REs

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

//vector[y_N[1]] y1_rep;
//vector[y_N[2]] y2_rep;
//for(i in 1:y_N[1]){
// y1_rep[i] = normal_rng(y1_etamu, exp(y1_etasigma));
//}
//for(i in 1:y_N[2]){
// y2_rep[i] = normal_rng(y2_etamu, exp(y2_etasigma));
//}
