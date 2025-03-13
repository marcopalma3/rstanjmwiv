#include pre.stan

// modified from https://github.com/sambrilleman/2018-StanCon-Notebook/blob/master/Stan/jm.stan

//----- Hyperparameters for prior distributions

// means for priors
// coefficients
vector[ymu_K[1]]     y1mu_prior_mean;
vector[ymu_K[2]]     y2mu_prior_mean;
vector[ysigma_K[1]]  y1sigma_prior_mean;
vector[ysigma_K[2]]  y2sigma_prior_mean;
vector[e_K]          e_prior_mean;
vector[a_K]          a_prior_mean;
vector[M]            ymu_prior_mean_for_intercept;
vector[M]            ysigma_prior_mean_for_intercept;
vector[basehaz_df]   e_prior_mean_for_aux;

// scale for priors
vector<lower=0>[ymu_K[1]] y1mu_prior_scale;
vector<lower=0>[ymu_K[2]] y2mu_prior_scale;
vector<lower=0>[ysigma_K[1]] y1sigma_prior_scale;
vector<lower=0>[ysigma_K[2]] y2sigma_prior_scale;
vector<lower=0>[e_K]    e_prior_scale;
vector<lower=0>[a_K]    a_prior_scale;
vector<lower=0>[M]      ymu_prior_scale_for_intercept;
vector<lower=0>[M]      ysigma_prior_scale_for_intercept;
vector<lower=0>[basehaz_df] e_prior_scale_for_aux;

// lkj prior stuff
vector<lower=0>[b_K] b_prior_scale;
vector<lower=0>[b_K] b_prior_df;
real<lower=0> b_prior_regularization;   //shape parameter of LKJ prior, https://mc-stan.org/docs/2_21/functions-reference/lkj-correlation.html
