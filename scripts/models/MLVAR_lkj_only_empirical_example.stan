////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of respondents
  int<lower=0> N_total; // number of total response measurements
  int<lower=0> P; // number of covariates
  int <lower=1, upper= 3> sparsity; // sparsity level affecting the hyper-priors for Beta
  int n_pc; // number of partial correlations
  int<lower = 0, upper = 1> mean_center; // mean centering for the predictors
  array[n_pc] int idx_rho; // index for partial correlations
  array[I] int<lower=0> n_t; // number of time points per person
  array[N_total] vector[K] Y; // longitudinal responses for all persons
  matrix[I,P] reg_covariate; // regression outcome
}
////////////////////////////////////////////////////////////////////////////////
transformed data{
  array[I*P] real reg_covariate_array=to_array_1d(to_vector(reg_covariate));
}

////////////////////////////////////////////////////////////////////////////////
parameters {
  // Temporal
  array[I] matrix[K,K] Beta_raw; // raw Beta matrix
  matrix[K,K] mu_Beta; // means of Betas
  matrix<lower=0>[K,K] sigma_Beta; // SDs of Betas
  matrix[I,K] Intercepts_raw; // raw intercepts
  vector[K] mu_Intercepts; // means of the intercepts
  vector<lower=0>[K] sigma_Intercepts; // SDs of the intercepts
  
  // Contemporaneous: Partial Correlations
  array[I] cholesky_factor_corr[K] L_Theta; // cholesky factor of the correlation matrix of innovations
  // Contemporaneous: Variances
  matrix<lower=0>[I,K] theta_sd; // SDs of the innovations, half-t
  // Regression
  vector[P] reg_intercept;   // Intercept of Regression
  vector[P] reg_slope_density; // Slope of Regression
  vector<lower=0>[P] sigma_residual;  // Residual term of Regression
  
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
  vector[I*n_pc] Rho_vec; // vectorize Rho

  //  Centrality for each individual
  vector[I]  Beta_density;
  vector[I]  Rho_density;
  matrix[I,K] Beta_in_strength;
  matrix[I,K] Beta_out_strength;
  matrix[I,K] Rho_strength;
  
  // Regression
  matrix[I,P] mu_regression;

   { 
  int position_rho = 1;
  // loop over respondents
  for(i in 1:I) {
  
  // Temporal //////////////////////////////////////////////////////////  
    // Implied: Beta ~ cauchy(mu_Beta, sigma_Beta)
    // 0.5*exp(sigma_Beta) implies that the mode of the hyper-prior 
    // for sigma_Beta is at 0.5

    if(sparsity == 1){
      // lowest sparsity level
      // wide normal prior on Beta means
      Beta[i] = mu_Beta + sigma_Beta .* Beta_raw[i]; 
    }
    if(sparsity == 2){
      // medium sparsity level
      // narrow, fat-tailed prior on Beta means
      Beta[i] = 0.05 * mu_Beta + sigma_Beta .* Beta_raw[i]; 
    }
    if(sparsity == 3){
      // high sparsity level
      // Beta means are fixed to zero
      // doubled variance of random effects to compensate fixed means
      Beta[i] = 0 + 2 .* sigma_Beta .* Beta_raw[i]; 
    }
    
    Intercepts[i,] = mu_Intercepts' + sigma_Intercepts' .* Intercepts_raw[i, ];
  
  // Contemporaneous ///////////////////////////////////////////////////
    
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
    Rho_vec[position_rho:(position_rho - 1 + n_pc)] = to_vector(Rho[i])[idx_rho];
    position_rho += n_pc; // increment position counter
    
    // Aggregate Statistics /////////////////////////////////////////////
    //  In-Strength, Out-Strength, and Centrality
    for(k in 1:K){
      
      // Ignore autoregressive effects in centrality estimation
      // by ignoring the diagonal elements
      Beta_in_strength[i,k] = 0; // Initialize to zero
      Beta_out_strength[i,k] = 0; // Initialize to zero
      Rho_strength[i,k] = 0; // Initialize to zero
      
      for(j in 1:K){
        if(j != k){
          Beta_in_strength[i,k]  += abs(Beta[i,k,j]); // colsums
          Beta_out_strength[i,k] += abs(Beta[i,j,k]); // rowsums
          Rho_strength[i,k]      += abs(Rho[i,j,k]); // sum of partial correlations
        }
      }
      // Divide by number of predictors to obtain the mean
      // Subtract 1 because we ignore the diagonal
      Beta_in_strength[i,k]  = Beta_in_strength[i,k] / (K -1);
      Beta_out_strength[i,k] = Beta_out_strength[i,k] / (K -1);
      Rho_strength[i,k]      = Rho_strength[i,k] / (K - 1);
      
      //Rho_centrality[i,k]    = mean(abs(Rho[i, ,k]));
    } // end k
    // Density
    Beta_density[i] = mean(abs(Beta[i]));
    Rho_density[i]  = mean(abs(Rho[i]));
  } // end i
  } // end block
  // Regression ////////////////////////////////////////////////////////
    mu_regression[,1] = reg_intercept[1] + reg_slope_density[1] * Beta_out_strength[,1];
    mu_regression[,2] = reg_intercept[2] + reg_slope_density[2] * Beta_density;

} // end transformed parameters
////////////////////////////////////////////////////////////////////////////////
model {
  
  // Priors Temporal //////////////////////////////////////////////////////////
  vector[I*K*K] Beta_raw_vec; // vectorize Beta_raw
  int position_Beta = 1; // position counter for Beta

  for (i in 1:I) {
    Beta_raw_vec[position_Beta:(position_Beta - 1 + K*K)] = to_vector(Beta_raw[i]);
    position_Beta += K*K; // increment position counter  
  } // end i
  Beta_raw_vec              ~ std_normal(); // prior on Beta
  if(sparsity == 1){
    to_vector(mu_Beta)      ~ normal(0, 1); // prior on mu_Beta
  }
  if(sparsity == 2){
    to_vector(mu_Beta)      ~ cauchy(0, 1); // prior on mu_Beta
  }
  to_vector(sigma_Beta)     ~ student_t(3, 0, .2); // prior on sigma_Beta, half-t
  to_vector(Intercepts_raw) ~ std_normal(); // prior on Intercepts
  mu_Intercepts             ~ student_t(3, 0, 2); // prior on mu_Intercepts
  sigma_Intercepts          ~ student_t(3, 0, .2); // prior on sigma_Intercepts, half-t
  
  // Priors Contemporaneous ///////////////////////////////////////////////////
  // Theta
  to_vector(theta_sd) ~ student_t(3, 0, 1); // prior on residual variances, half-t

  // Regression
  reg_intercept     ~ student_t(3, 0, 2); // prior on reg_intercept
  reg_slope_density ~ student_t(3, 0, 2); // prior on reg_slope_density
  sigma_residual    ~ student_t(3, 0, 1); // prior on sigma_residual
  
  
  int position_Y = 1; // position counter for the data
  for (i in 1:I) {
    // Precision Matrix prior
    L_Theta[i] ~ lkj_corr_cholesky(3); // prior on L_Theta
    
    //// Likelihood //////////////////////////////////////////////////////////
    
    // Partition data for one respondent
    array[n_t[i]] vector[K] Y_temp = segment(Y, position_Y, n_t[i]); // slice array
    position_Y += n_t[i]; // increment position counter
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(theta_sd[i], L_Theta[i]);
    array[n_t[i]-1] vector[K] mu_network;
    
    // network predictions: loop over time points
    // depends on setting by user
    if(mean_center == 1){
      for(t in 1:(n_t[i]-1)){
        mu_network[t] = 
          to_vector(Intercepts[i, ]) + Beta[i] * (Y_temp[t,] - to_vector(Intercepts[i, ])); // predictions for the network
        } // end t
    }
    if(mean_center == 0){
      for(t in 1:(n_t[i]-1)){
        mu_network[t] = 
          to_vector(Intercepts[i, ]) + Beta[i] * Y_temp[t,]; // predictions for the network
        } // end t
    }
      
    // Network
    Y_temp[2:n_t[i]] ~ multi_normal_cholesky(mu_network, Sigma_chol);
  } // end i
  
  
  // Regression
    reg_covariate_array ~ normal(
      to_array_1d(to_vector(mu_regression)), to_array_1d(to_vector(rep_matrix(sigma_residual', I)))
      );

    
} // end model
////////////////////////////////////////////////////////////////////////////////
generated quantities{
  vector[P] reg_slope_density_z;
  vector[P] reg_intercept_z;
  
  // Standardize Regression Coefficients
  // beta * sd(x) = beta_std;
  reg_slope_density_z[1] = reg_slope_density[1] * sd(Beta_out_strength[,1]);
  reg_slope_density_z[2] = reg_slope_density[2] * sd(Beta_density);
  
  // Standardize Regression Intercept
  // alpha_std = alpha + (beta_std * mean(x)) / sd(x);
  reg_intercept_z[1] = reg_intercept[1] + (reg_slope_density_z[1] * mean(Beta_out_strength[,1])) / sd(Beta_out_strength[,1]);
  reg_intercept_z[2] = reg_intercept[2] + (reg_slope_density_z[2] * mean(Beta_density)) / sd(Beta_density);

  // Posterior predictive checks
  // Posterior predictive checks
  int position_Y = 1; // position counter 
  array[N_total] vector[K] Y_rep; // posterior predictive time series
  
  for (i in 1:I) {
    // Partition observed data
    array[n_t[i]] vector[K] Y_temp = segment(Y, position_Y, n_t[i]); // slice array
    position_Y += n_t[i]; 

    // Cholesky decomposition
    matrix[K, K] Sigma_chol = diag_pre_multiply(theta_sd[i], L_Theta[i]);
    array[n_t[i]-1] vector[K] mu_network;

    // Compute network predictions for each time point
    if (mean_center == 1) {
      for (t in 1:(n_t[i]-1)) {
        mu_network[t] = 
          to_vector(Intercepts[i, ]) + Beta[i] * (Y_temp[t,] - to_vector(Intercepts[i, ]));
      }
    }
    if (mean_center == 0) {
      for (t in 1:(n_t[i]-1)) {
        mu_network[t] = 
          to_vector(Intercepts[i, ]) + Beta[i] * Y_temp[t,];
      }
    }

    // Simulate replicated data
    for (t in 2:n_t[i]) {
      Y_rep[position_Y - n_t[i] + t - 1] = 
        multi_normal_cholesky_rng(mu_network[t-1], Sigma_chol);
    }
  }

  // Replicated regression outcomes
  matrix[I,P] reg_covariate_rep;
  for (i in 1:I) {
    for (p in 1:P) {
      reg_covariate_rep[i, p] = normal_rng(mu_regression[i, p], sigma_residual[p]);
    }
  }
  
} // end generated quantities


  