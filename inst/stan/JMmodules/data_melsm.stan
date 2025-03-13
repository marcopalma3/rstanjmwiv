#include pre.stan

// modified from https://github.com/sambrilleman/2018-StanCon-Notebook/blob/master/Stan/jm.stan

//----- Longitudinal submodels

// number of long. submodels
// NB this is fixed equal to 2 for this simplified jm.stan 
// file. See the jm.stan file in rstanarm for more general 
// code that allows for 1, 2 or 3 longitudinal submodels.
int<lower=1,upper=2> M; 

// population level dimensions
int<lower=0> y_N[M]; // num observations
int<lower=0> ymu_K[M]; // num predictors in the mean submodels
int<lower=0> ysigma_K[M]; // num predictors in the variability submodels

// population level data
// NB these design matrices are evaluated at the observation times
vector[y_N[1]] y1; // response vectors
vector[y_N[2]] y2;  
matrix[y_N[1],ymu_K[1]] y1mu_X; // fe (centered) design matrix for mean submodel
matrix[y_N[2],ymu_K[2]] y2mu_X; 
matrix[y_N[1],ysigma_K[1]] y1sigma_X; // fe (centered) design matrix for the variability submodel
matrix[y_N[2],ysigma_K[2]] y2sigma_X; 
vector[ymu_K[1]] y1mu_Xbar; // predictor means for mean submodel
vector[ymu_K[2]] y2mu_Xbar;
vector[ysigma_K[1]] y1sigma_Xbar; // predictor means for the variability submodel
vector[ysigma_K[2]] y2sigma_Xbar;



// group level dimensions
int<lower=0> b_N;     // num groups
int<lower=0> b_K; // total num params
int<lower=0> bmu_K[M]; // num params in each mean submodel
int<lower=0> bsigma_K[M]; // num params in each variability submodel



// group level data
vector[y_N[1]] y1mu_Z[bmu_K[1]]; // re design matrix for mean submodels
vector[y_N[2]] y2mu_Z[bmu_K[2]];
vector[y_N[1]] y1sigma_Z[bsigma_K[1]]; // re design matrix for variability submodels
vector[y_N[2]] y2sigma_Z[bsigma_K[2]];
int<lower=0> y1_Z_id[y_N[1]]; // group ids for y*_Z - they are equal for mu and sigma submodels of the same longitudinal marker
int<lower=0> y2_Z_id[y_N[2]];
