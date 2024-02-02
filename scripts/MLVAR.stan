////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of predictors
  int<lower=0> N_total; // number of total measurements
  array[I] int<lower=0> n_t; // number of time points per person
  array[N_total] vector[K] Y; // responses
  vector[I] outcome; // outcome
}
transformed data{
  int n_pc = K * (K-1) %/% 2; // number of off-diagonal elements
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  array[I] matrix<lower=-pi()/2, upper=pi()/2>[K,K] Beta_raw; // raw Beta matrix
  matrix[K,K] mu_Beta; // means of Betas
  matrix[K,K] sigma_Beta; // SDs of Betas
  
  // Contemporaneous: Partial Correlations
  array[I] cholesky_factor_corr[K] L_Theta; // cholesky factor of the correlation matrix of innovations

  // Contemporaneous: Variances
  matrix[I,K] theta_sd_raw; // SDs of the innovations
  row_vector[K] mu_theta_sd; // means of the innovation SDs
  row_vector[K] sigma_theta_sd; // SDs of the innovation SDs

  
} // end parameters
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  array[I] matrix[K,K] Beta;
  // Covariance matrix from cholesky corr matrix and SDs
  array[I]matrix[K,K] Sigma; 
  // Partial correlation matrix

  matrix[I,K] theta_sd;

  // loop over respondents
  for(i in 1:I) {
  
  // Temporal //////////////////////////////////////////////////////////  
    Beta[i] = mu_Beta + exp(sigma_Beta + log(.5)) .* tan(Beta_raw[i]); // Beta ~ cauchy(mu_Beta, sigma_Beta)
  
  // Contemporaneous ///////////////////////////////////////////////////
    
    // Covariance matrix from cholesky corr matrix and SDs
    theta_sd[i,] = exp(mu_theta_sd + exp(sigma_theta_sd + log(.5)) .* theta_sd_raw[i,]);
    Sigma[i] = diag_pre_multiply(exp(theta_sd[i] + log(.5)), L_Theta[i]) * 
               diag_pre_multiply(exp(theta_sd[i] + log(.5)), L_Theta[i])'; 
    // Precision matrix from covariance matrix
    matrix[K,K] Theta = inverse_spd(Sigma[i]); 
    
  } // end i
} // end transformed parameters
////////////////////////////////////////////////////////////////////////////////
model {
  // Priors Temporal
  vector[I*K*K] Beta_raw_vec; // vectorize Beta_raw
  int position_Beta = 1; // position counter for Beta
  for (i in 1:I) {
    Beta_raw_vec[position_Beta:(position_Beta - 1 + K*K)] = to_vector(Beta_raw[i]);
    position_Beta += K*K; // increment position counter  
  } // end i
  target+= uniform_lpdf(Beta_raw_vec | -pi()/2, pi()/2); // prior on Beta
  target+= std_normal_lpdf(to_vector(mu_Beta)); // prior on mu_Beta
  target+= std_normal_lpdf(to_vector(sigma_Beta)); // prior on sigma_Beta
  
  // Priors Contemporaneous
  // Theta
  target+= std_normal_lpdf(to_vector(theta_sd_raw)); // prior on sigma_theta
  target+= std_normal_lpdf(mu_theta_sd);             // prior on mu_theta_sd
  target+= std_normal_lpdf(sigma_theta_sd);          // prior on sigma_theta_sd


  int position_Y = 1; // position counter for the data
  for (i in 1:I) {
    // Precision Matrix
    target+= lkj_corr_cholesky_lpdf(L_Theta[i] | 1); 

    // Partition data for one respondent
    array[n_t[i]] vector[K] Y_temp = segment(Y, position_Y, n_t[i]); // slice array
    position_Y += n_t[i]; // increment position counter
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(exp(theta_sd[i]), L_Theta[i]);
    array[n_t[i]-1] vector[K] mu_network;
    // loop over time points
    for(t in 1:(n_t[i]-1)){
      mu_network[t] = Beta[i] * Y_temp[t,]; // predictions for the network
      } // end t
    
    // Network
    target += multi_normal_cholesky_lpdf(Y_temp[2:n_t[i]] | mu_network, Sigma_chol);
  } // end i
} // end model
////////////////////////////////////////////////////////////////////////////////
