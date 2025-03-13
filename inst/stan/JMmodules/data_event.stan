#include pre.stan

// modified from https://github.com/sambrilleman/2018-StanCon-Notebook/blob/master/Stan/jm.stan

//----- Event submodel

// data for calculating event submodel linear predictor in GK quadrature
// NB these design matrices are evaluated AT the event time and
// the (unstandardised) quadrature points
int<lower=0> e_K;           // num. of predictors in event submodel
int<lower=0> Npat;          // num. individuals
int<lower=0> Nevents;       // num. events (ie. not censored)  
int<lower=0> qnodes;        // num. of nodes for GK quadrature 
int<lower=0> Npat_times_qnodes; 
int<lower=0> nrow_e_Xq;     // num. rows in event submodel predictor matrix
vector[nrow_e_Xq] e_times;  // event times and quadrature points
matrix[nrow_e_Xq,e_K] e_Xq; // centered predictor matrix (event submodel)
vector[e_K] e_Xbar;         // predictor means (event submodel)
int<lower=0> basehaz_df;    // df for B-splines baseline hazard
matrix[nrow_e_Xq,basehaz_df] basehaz_X; // design matrix (basis terms) for baseline hazard
vector[Npat_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling 
real norm_const;            // constant shift for log baseline hazard
real assoc_code;            // association type
int<lower=0> a_K;           // num. of association parameters

// data for calculating long. submodel linear predictor in GK quadrature
// NB these design matrices are evaluated AT the event time and
// the (unstandardised) quadrature points
int<lower=0> nrow_y_Xq[M]; // num. rows in long. predictor matrix at quadpoints
matrix[nrow_y_Xq[1],ymu_K[1]] y1mu_Xq; // fe design matrix at quadpoints (mean submodel)
matrix[nrow_y_Xq[2],ymu_K[2]] y2mu_Xq; 
matrix[nrow_y_Xq[1],ysigma_K[1]] y1sigma_Xq; // fe design matrix at quadpoints (variability submodel)
matrix[nrow_y_Xq[2],ysigma_K[2]] y2sigma_Xq; 
vector[nrow_y_Xq[1]] y1mu_Zq[bmu_K[1]]; // re design matrix at quadpoints (mean submodel)
vector[nrow_y_Xq[2]] y2mu_Zq[bmu_K[2]];
vector[nrow_y_Xq[1]] y1sigma_Zq[bsigma_K[1]]; // re design matrix at quadpoints (variability submodel)
vector[nrow_y_Xq[2]] y2sigma_Zq[bsigma_K[2]];
int<lower=0> y1_Zq_id[nrow_y_Xq[1]]; // group indexing for re design matrix
int<lower=0> y2_Zq_id[nrow_y_Xq[2]]; 
