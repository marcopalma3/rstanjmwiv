#include pre.stan

//transformed data {
  int<lower=0> b_KM[2*M] = {bmu_K[1], bsigma_K[1], bmu_K[2], bsigma_K[2]};//num params in each submodel (mu1, sigma1, mu2, sigma2)
  // bmu_K must be equal to elements 1 and 3 of b_KM
  // bsigma_K must be equal to elements 2 and 4 of b_KM


  // indexing used to extract lower tri of RE covariance matrix
  int b_cov_idx[b_K + choose(b_K, 2)];
  if (b_K > 0) 
    b_cov_idx = lower_tri_indices(b_K);
//}
