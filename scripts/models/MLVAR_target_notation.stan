////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of respondents
  int<lower=0> N_total; // number of total response measurements
  int n_pc; // number of partial correlations
  array[n_pc] int idx_rho; // index for partial correlations
  array[I] int<lower=0> n_t; // number of time points per person
  array[N_total] vector[K] Y; // longitudinal responses for all persons
  array[I] real reg_covariate; // regression outcome
}
////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  array[I] matrix<lower=-pi()/2, upper=pi()/2>[K,K] Beta_raw; // raw Beta matrix
  matrix[K,K] mu_Beta; // means of Betas
  matrix[K,K] sigma_Beta; // SDs of Betas
  matrix[I,K] Intercepts_raw; // raw intercepts
  vector[K] mu_Intercepts; // means of the intercepts
  vector[K] sigma_Intercepts; // SDs of the intercepts
  
  // Contemporaneous: Partial Correlations
  array[I] cholesky_factor_corr[K] L_Theta; // cholesky factor of the correlation matrix of innovations
  matrix[I, n_pc] rho_loc_raw; // means of the partial correlations
  matrix[I, n_pc] rho_var_raw;// SDs of the partial correlations
  row_vector[n_pc] mu_rho_loc; // population mean of the partial correlations locations
  row_vector[n_pc] sigma_rho_loc; // population SD of the partial correlations locations
  row_vector[n_pc] mu_rho_var; // population mean of the partial correlations variances
  row_vector[n_pc] sigma_rho_var; // population SD of the partial correlations variances
  // Contemporaneous: Variances
  matrix[I,K] theta_sd_raw; // SDs of the innovations
  row_vector[K] mu_theta_sd; // means of the innovation SDs
  row_vector[K] sigma_theta_sd; // SDs of the innovation SDs
  // Regression
  real reg_intercept;   // Intercept of Regression
  real reg_slope_density; // Slope of Regression
  real reg_residual;  // Residual term of Regression
  
} // end parameters
////////////////////////////////////////////////////////////////////////////////
transformed parameters{
  // Non-centered parameterization for Beta matrix
  array[I] matrix[K,K] Beta;  // estimated Beta coefficients
  matrix[I,K] Intercepts; // raw intercepts
  
  // Covariance matrix from cholesky corr matrix and SDs
  array[I]matrix[K,K] Sigma; 
  // Partial correlation matrix
  array[I] matrix[K,K] Rho;
  matrix[I, n_pc] rho_loc;   // location of partial correlations
  matrix[I, n_pc] rho_var;   // scale of partial correlations
  matrix[I,K] theta_sd;      // 
  //  Centrality for each individual
  array[I] vector[K] Beta_out_strength; 
  array[I] vector[K] Beta_in_strength;
  array[I] vector[K] Rho_centrality;
  array[I] real Beta_density;
  array[I] real Rho_density;
  // Regression
  array[I] real mu_regression;

  
  // loop over respondents
  for(i in 1:I) {
  
  // Temporal //////////////////////////////////////////////////////////  
    // Implied: Beta ~ cauchy(mu_Beta, sigma_Beta)
    // exp(sigma_Beta + log(.5)) implies that the mode of the hyper-prior 
    // for sigma_Beta is at 0.5
    Beta[i] = mu_Beta + exp(sigma_Beta + log(.5)) .* tan(Beta_raw[i]); 
    Intercepts[i,] = mu_Intercepts' + exp(sigma_Intercepts' + log(.5)) .* Intercepts_raw[i, ];
  
  // Contemporaneous ///////////////////////////////////////////////////
    
    // Covariance matrix from cholesky corr matrix and SDs
    theta_sd[i,] = exp(mu_theta_sd + exp(sigma_theta_sd + log(.5)) .* theta_sd_raw[i,]);
    // Correlation Matrix
    Sigma[i] = multiply_lower_tri_self_transpose(L_Theta[i]); 
    // Precision matrix from covariance matrix
    matrix[K,K] Theta = inverse_spd(Sigma[i, , ]); 
    // Partial correlation matrix
    // computed from off-diagonal elements of precision matrix
    for(j in 1:K){
      for(k in 1:K){
        if(j != k){
          Rho[i,j,k] = -Theta[j,k] / sqrt(Theta[j,j] * Theta[k,k]);
        }else{
          Rho[i,j,k] = 0;
        } // end else
      } // end k
    } // end j
    
    // Rho priors non-centered parameterization: inplied cauchy priors
    rho_loc[i,] = Phi_approx(mu_rho_loc + exp(sigma_rho_loc + log(.5)) .* rho_loc_raw[i,]);
    rho_var[i,] = exp(mu_rho_var + exp(sigma_rho_var + log(.5)) .* rho_var_raw[i,]);
    
    // Aggregate Statistics /////////////////////////////////////////////
    //  In-Strength, Out-Strength, and Centrality
    for(k in 1:K){
      Beta_out_strength[i,k] = mean(abs(Beta[i, ,k]));
      Beta_in_strength[i,k]  = mean(abs(Beta[i,k, ]));
      Rho_centrality[i, k]   = mean(abs(Rho[i, ,k]));
    } // end k
    // Density
    Beta_density[i] = mean(abs(Beta[i]));
    Rho_density[i]  = mean(abs(Rho[i]));
    
    // Regression ////////////////////////////////////////////////////////
    mu_regression[i] = reg_intercept +  reg_slope_density * Beta_density[i];
   // mu_regression[i] = reg_intercept +  reg_slope_density * reg_covariate[i];
    
  } // end i
} // end transformed parameters
////////////////////////////////////////////////////////////////////////////////
model {
  
  // Priors Temporal //////////////////////////////////////////////////////////
  vector[I*K*K] Beta_raw_vec; // vectorize Beta_raw
  int position_Beta = 1; // position counter for Beta
  vector[I*n_pc] Rho_vec; // vectorize Rho
  int position_rho = 1; // position counter for Rho
  for (i in 1:I) {
    Beta_raw_vec[position_Beta:(position_Beta - 1 + K*K)] = to_vector(Beta_raw[i]);
    position_Beta += K*K; // increment position counter  
    Rho_vec[position_rho:(position_rho - 1 + n_pc)] = to_vector(Rho[i])[idx_rho];
    position_rho += n_pc; // increment position counter
  } // end i
  target+= uniform_lpdf(Beta_raw_vec | -pi()/2, pi()/2); // prior on Beta
  target+= std_normal_lpdf(to_vector(mu_Beta)); // prior on mu_Beta
  target+= std_normal_lpdf(to_vector(sigma_Beta)); // prior on sigma_Beta
  target += std_normal_lpdf(to_vector(Intercepts_raw)); // prior on Intercepts
  target += std_normal_lpdf(mu_Intercepts); // prior on mu_Intercepts
  target += std_normal_lpdf(sigma_Intercepts); // prior on sigma_Intercepts
  
  // Priors Contemporaneous ///////////////////////////////////////////////////
  // Theta
  target+= std_normal_lpdf(to_vector(theta_sd_raw)); // prior on sigma_theta
  target+= std_normal_lpdf(mu_theta_sd);             // prior on mu_theta_sd
  target+= std_normal_lpdf(sigma_theta_sd);          // prior on sigma_theta_sd
  
  // Partial correlations 
  target+= beta_proportion_lpdf(Rho_vec / 2 + 0.5 | to_vector(rho_loc'), to_vector(rho_var'));
  target+= std_normal_lpdf(to_vector(rho_loc_raw)); // prior on rho_loc_raw
  target+= std_normal_lpdf(to_vector(rho_var_raw)); // prior on rho_var_raw
  target+= std_normal_lpdf(mu_rho_loc); // prior on mu_rho_loc
  target+= std_normal_lpdf(sigma_rho_loc); // prior on sigma_rho_loc
  target+= std_normal_lpdf(mu_rho_var); // prior on mu_rho_var
  target+= std_normal_lpdf(sigma_rho_var); // prior on sigma_rho_var
  
  // Regression
  target+= student_t_lpdf(reg_intercept | 3, 0, 2); // prior on reg_intercept
  target+= student_t_lpdf(reg_slope_density | 3, 0, 2); // prior on reg_slope_density
  target+= student_t_lpdf(reg_residual | 3, 0, 1); // prior on reg_residual
  
  
  int position_Y = 1; // position counter for the data
  for (i in 1:I) {
    // Precision Matrix prior
    target+= lkj_corr_cholesky_lpdf(L_Theta[i] | 1); // prior on L_Theta
    
    //// Likelihood //////////////////////////////////////////////////////////
    
    // Partition data for one respondent
    array[n_t[i]] vector[K] Y_temp = segment(Y, position_Y, n_t[i]); // slice array
    position_Y += n_t[i]; // increment position counter
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(exp(theta_sd[i] + log(.5)), L_Theta[i]);
    array[n_t[i]-1] vector[K] mu_network;
    
    // network predictions: loop over time points
    for(t in 1:(n_t[i]-1)){
      mu_network[t] = Intercepts[i, ]' + Beta[i] * Y_temp[t,]; // predictions for the network
      } // end t
      
    // Network
    target += multi_normal_cholesky_lpdf(Y_temp[2:n_t[i]] | mu_network, Sigma_chol);
  } // end i
  
  // Regression
    target += normal_lpdf(reg_covariate | mu_regression, exp(reg_residual));
    // target += normal_lpdf(Beta_density | mu_regression, exp(reg_residual));
    
} // end model
////////////////////////////////////////////////////////////////////////////////
