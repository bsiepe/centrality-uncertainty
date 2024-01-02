////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of predictors
  int<lower=0> N; // number of total measurements
  
  array[I] int<lower=0> n_t; // number of time points per person

  array[N] vector[K] Y; // responses
  
  // Priors
  matrix[K,K] prior_Beta_loc; // locations for priors on Beta matrix
  matrix[K,K] prior_Beta_scale;  // scales for priors on Beta matrix
  matrix[K,K] prior_Rho_loc;  // locations for priors on partial correlations
  matrix[K,K] prior_Rho_scale;   // scales for priors on partial correlations
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  array[I] matrix[K,K] Beta_raw; // raw Beta matrix
  //real mu_Beta;
  //real<lower=0> sigma_Beta;
  
  // Contemporaneous
  array[I] cholesky_factor_corr[K] L_Theta;
  array[I] vector[K] sigma_theta;
} // end parameters
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  array[I] matrix[K,K] Beta;
  //matrix[K,K] Beta;
  
  // Covariance matrix from cholesky corr matrix and SDs
  array[I]matrix[K,K] Sigma; 
  // Partial correlation matrix
  array[I]matrix[K,K] Rho;
  //  Centrality
  array[I] vector[K] Beta_out_strength; 
  array[I] vector[K] Beta_in_strength;
  array[I] vector[K] Beta_density;
  array[I] vector[K] Rho_centrality;
  array[I] vector[K] Rho_density;

  for(i in 1:I) {
    
    Beta[i] = Beta_raw[i] .* prior_Beta_scale + prior_Beta_loc;
  //Beta[I] = Beta_raw * sigma_Beta + mu_Beta;
  
  // Covariance matrix from cholesky corr matrix and SDs
    Sigma[i] = diag_pre_multiply(exp(sigma_theta[i]), L_Theta[i]) * 
               diag_pre_multiply(exp(sigma_theta[i]), L_Theta[i])'; 
               
    // Precision matrix
    matrix[K,K] Theta = inverse_spd(Sigma[i]); 
    for(j in 1:K){
      for(k in 1:K){
        if(j != k){
          Rho[i,j,k] = -Theta[j,k] / sqrt(Theta[j,j] * Theta[k,k]);
        }else{
          Rho[i,j,k] = 0;
        } // end else
      } // end k
    } // end j
    
    //  Centrality, Density
    for(k in 1:K){
      Beta_out_strength[i,k] = sum(Beta[i, ,k]);
      Beta_out_strength[i,k] = sum(Beta[i,k ,]);
      Beta_density[i,k]      = sum(Beta[i, , ]);
      Rho_centrality[i, k]   = sum(Rho[i, ,k]);
      Rho_density[i, k]      = sum(Rho[i, , ]);
    } // end k
  } // end i
} // end transformed parameters
////////////////////////////////////////////////////////////////////////////////
model {
  int position = 1; // position counter
  for (i in 1:I) {
      // Priors
  target+=   std_normal_lpdf(to_vector(Beta_raw[i]));    // prior on Beta
  //target+= student_t_lpdf(mu_Beta | 3,0,2);
  //target+= student_t_lpdf(sigma_Beta | 3,0,2);
  target+=   lkj_corr_cholesky_lpdf(L_Theta[i] | 1);  // prior on Cholesky factor
  target+=   student_t_lpdf(sigma_theta[i] | 3,0,2);   // prior on sigma_theta
  // Priors on partial correlations
  for(j in 1:K){
    for(k in 1:K){
      if(j < k){
        // Scaled beta prior on partial correlations 
        // (Rho[j,k] / 2 + 0.5 is a scaling to the unit interval)
        target+= beta_proportion_lpdf(
          Rho[i,j,k] / 2 + 0.5 | prior_Rho_loc[j,k], prior_Rho_scale[j,k]);
        } // end if
      } // end k
    } // end j
    
    array[n_t[i]] vector[K] Y_temp = segment(Y, position, n_t[i]); // slice array
    
    
    position = position + n_t[i]; // increment position counter
    
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(exp(sigma_theta[i]), L_Theta[i]);
    for(t in 2:n_t[i]){
      // BS: What about intercept?
      vector[K] mu = Beta[i] * Y_temp[t-1,];
      target += multi_normal_cholesky_lpdf(Y_temp[t,] | mu, Sigma_chol);
    } // end t
  } // end i
} // end model
////////////////////////////////////////////////////////////////////////////////
// generated quantities{
//   vector[T-1] log_lik;
//   {
//     // Cholesky decomposition of the covariance matrix
//     matrix[K, K] Sigma_chol = diag_pre_multiply(exp(sigma_theta), L_Theta);
//     for(t in 2:T){
//       // BS: What about intercept?
//       vector[K] mu = Beta * Y[t-1,];
//       log_lik[t-1] = multi_normal_cholesky_lpdf(Y[t, ] | mu, Sigma_chol);
//     }
//   }
// }
