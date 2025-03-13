#include pre.stan

/** 
* Evaluate the linear predictor for the glmer submodel
*
* @param X Design matrix for fe
* @param Z Design matrix for re, for a single grouping factor
* @param Z_id Group indexing for Z
* @param gamma The intercept parameter
* @param beta Vector of population level parameters
* @param bMat Matrix of group level params
* @param shift Number of columns in bMat
*   that correpond to group level params from prior glmer submodels
* @return A vector containing the linear predictor for the glmer submodel
*/  
vector evaluate_eta(matrix X, vector[] Z, int[] Z_id, real gamma, 
vector beta, matrix bMat, int shift) {
  int N = rows(X);    // num rows in design matrix
  int K = rows(beta); // num predictors
  int p = size(Z);    // num group level params
  vector[N] eta;
  
  if (K > 0) eta = X * beta + gamma;    //ADDED: gamma was missing in Brilleman's version
  else eta = rep_vector(0.0, N) + gamma;
  
  for (k in 1:p)
  for (n in 1:N)
  eta[n] = eta[n] + (bMat[Z_id[n], k + shift]) * Z[k,n];
  
  return eta;
} 

/** 
* Get the indices corresponding to the lower tri of a square matrix
*
* @param dim The number of rows in the square matrix
* @return A vector of indices
*/ 
int[] lower_tri_indices(int dim) {
  int indices[dim + choose(dim, 2)];
  int mark = 1;
  for (r in 1:dim) {
    for (c in r:dim) {
      indices[mark] = (r - 1) * dim + c;
      mark = mark + 1;
    }
  }
  return indices;
}





/** function for random effects association VALID ONLY FOR RANDOM INTERCEPT!
*
* @return The matrix of random effects where each row is corresponding to Z_id[n] 
*/
matrix evaluate_bMat_q(matrix bMat, int[] Z_id, int nrow) {
  int bK = cols(bMat);
  matrix[nrow, bK] bMat_q; 
  
  for (i in 1:nrow) {
    bMat_q[i, ] = bMat[Z_id[i], ];
  }
  
  return bMat_q;
}


