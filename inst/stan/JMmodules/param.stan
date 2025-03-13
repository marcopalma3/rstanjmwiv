#include pre.stan

real y1mu_gamma; // primitive intercepts in long. submodels
real y1sigma_gamma; 
real y2mu_gamma; 
real y2sigma_gamma;
vector[ymu_K[1]] y1mu_z_beta; // primitive coefs in long. submodels
vector[ysigma_K[1]] y1sigma_z_beta; 
vector[ymu_K[2]] y2mu_z_beta; 
vector[ysigma_K[2]] y2sigma_z_beta;
// group level params   
vector<lower=0>[b_K] b_sd; // group level sds  
matrix[b_K,b_N] z_b_mat;   // unscaled group level params 
cholesky_factor_corr[b_K > 1 ? b_K : 0] b_cholesky; // cholesky factor of corr matrix



vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)
vector[a_K] a_z_beta; // primitive assoc params (log hazard ratios)
//real<lower=0> y1_aux_unscaled; // unscaled residual error SDs - in MELSM we are modelling them
//real<lower=0> y2_aux_unscaled; 
vector[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard      

