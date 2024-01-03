////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of predictors
  int<lower=0> N; // number of total measurements
  array[I] int<lower=0> n_t; // number of time points per person
  array[N] vector[K] Y; // responses
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  array[I] matrix<lower=-pi()/2, upper=pi()/2>[K,K] Beta_raw; // raw Beta matrix
  vector[I] mu_Beta;
  vector[I] sigma_Beta;
  
  // Contemporaneous
  array[I] cholesky_factor_corr[K] L_Theta;
  matrix[I,K] sigma_theta;
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
  array[I] real Beta_density;
  array[I] vector[K] Rho_centrality;
  array[I] real Rho_density;

  for(i in 1:I) {
    Beta[i] = mu_Beta[i] + exp(sigma_Beta[i] + log(.5)) .* tan(Beta_raw[i]); // Beta ~ cauchy(mu_Beta, sigma_Beta)
  
    // Covariance matrix from cholesky corr matrix and SDs
    Sigma[i] = diag_pre_multiply(exp(sigma_theta[i] + log(.5)), L_Theta[i]) * 
               diag_pre_multiply(exp(sigma_theta[i] + log(.5)), L_Theta[i])'; 
               
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
      Beta_in_strength[i,k]  = sum(Beta[i,k, ]);
      Rho_centrality[i, k]   = sum(Rho[i, ,k]);
    } // end k
    // Density
    Beta_density[i] = sum(Beta[i]);
    Rho_density[i]  = sum(Rho[i]);
  } // end i
} // end transformed parameters
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors  
  {
  vector[I*K*K] Beta_raw_vec;
  int start = 1;
  for (i in 1:I) {
    Beta_raw_vec[start:(start - 1 + K*K)] = to_vector(Beta_raw[i]);
    start = start + K*K;
  }
  target+= uniform_lpdf(Beta_raw_vec | -pi()/2, pi()/2); // prior on Beta
  }
  target+= std_normal_lpdf(mu_Beta);
  target+= std_normal_lpdf(sigma_Beta);
  target+= std_normal_lpdf(to_vector(sigma_theta));   // prior on sigma_theta
    
  int position = 1; // position counter
  for (i in 1:I) {
    // Priors
    // marginal beta: alpha = beta = eta -1 + K/2
    // cholesky prior: eta = alpha +1 -K/2
    target+= lkj_corr_cholesky_lpdf(L_Theta[i] | 10 + 1 - K/2.0);  
    
    // Data for one respondent
    array[n_t[i]] vector[K] Y_temp = segment(Y, position, n_t[i]); // slice array
    position = position + n_t[i]; // increment position counter
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(exp(sigma_theta[i]), L_Theta[i]);
    array[n_t[i]-1] vector[K] mu;
    for(t in 1:(n_t[i]-1)){
      mu[t] = Beta[i] * Y_temp[t,];
    } // end t
    target += multi_normal_cholesky_lpdf(Y_temp[2:n_t[i]] | mu, Sigma_chol);
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
