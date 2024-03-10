////////////////////////////////////////////////////////////////////////////////
// VAR-Model with Custom Priors
////////////////////////////////////////////////////////////////////////////////
data {
  int<lower=0> K; // number of predictors
  int<lower=0> I; // number of respondents
  int<lower=0> N_total; // number of total response measurements
  int<lower=0> P; // number of covariates
  int n_pc; // number of partial correlations
  array[n_pc] int idx_rho; // index for partial correlations
  array[I] int<lower=0> n_t; // number of time points per person
  array[N_total] vector[K] Y; // longitudinal responses for all persons
  matrix[I,P] reg_covariate; // regression outcome
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
  // Contemporaneous: Variances
  matrix[I,K] theta_sd_raw; // SDs of the innovations
  row_vector[K] mu_theta_sd; // means of the innovation SDs
  row_vector[K] sigma_theta_sd; // SDs of the innovation SDs
  // Regression
  vector[P] reg_intercept;   // Intercept of Regression
  vector[P] reg_slope_density; // Slope of Regression
  vector[P] reg_residual;  // Residual term of Regression
  
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
    
  matrix[I,K] theta_sd;      // 
  //  Centrality for each individual
  vector[I]  Beta_density;
  vector[I]  Rho_density;
  matrix[I,K] Beta_out_strength;
  matrix[I,K] Beta_in_strength;
  matrix[I,K] Rho_centrality;
  
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
    Beta[i] = mu_Beta + 0.5*exp(sigma_Beta) .* tan(Beta_raw[i]); 
    Intercepts[i,] = mu_Intercepts' + 0.5*exp(sigma_Intercepts') .* Intercepts_raw[i, ];
  
  // Contemporaneous ///////////////////////////////////////////////////
    
    // Covariance matrix from cholesky corr matrix and SDs
    theta_sd[i,] = exp(mu_theta_sd + 0.5*exp(sigma_theta_sd) .* theta_sd_raw[i,]);
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
      Beta_out_strength[i,k] = mean(abs(Beta[i, ,k]));
      Beta_in_strength[i,k]  = mean(abs(Beta[i,k, ]));
      Rho_centrality[i,k]    = mean(abs(Rho[i, ,k]));
    } // end k
    // Density
    Beta_density[i] = mean(abs(Beta[i]));
    Rho_density[i]  = mean(abs(Rho[i]));
  } // end i
   }
  // Regression ////////////////////////////////////////////////////////
    mu_regression[,1] = reg_intercept[1] + reg_slope_density[1] * Beta_density;
    mu_regression[,2] = reg_intercept[2] + reg_slope_density[2] * Beta_density;
    mu_regression[,3] = reg_intercept[3] + reg_slope_density[3] * Beta_density;
    
    
    mu_regression[,4] = reg_intercept[4] + reg_slope_density[4] * Rho_density;
    mu_regression[,5] = reg_intercept[5] + reg_slope_density[5] * Rho_density;
    mu_regression[,6] = reg_intercept[6] + reg_slope_density[6] * Rho_density;
    
    mu_regression[,7] = reg_intercept[7] + reg_slope_density[7] * Beta_out_strength[,1];
    mu_regression[,8] = reg_intercept[8] + reg_slope_density[8] * Beta_out_strength[,1];
    mu_regression[,9] = reg_intercept[9] + reg_slope_density[9] * Beta_out_strength[,1];

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
  Beta_raw_vec ~ uniform(-pi()/2, pi()/2); // prior on Beta
  to_vector(mu_Beta) ~ std_normal(); // prior on mu_Beta
  to_vector(sigma_Beta) ~ std_normal(); // prior on sigma_Beta
  to_vector(Intercepts_raw) ~ std_normal(); // prior on Intercepts
  mu_Intercepts ~ student_t(3, 0, 2); // prior on mu_Intercepts
  sigma_Intercepts ~ std_normal(); // prior on sigma_Intercepts
  
  // Priors Contemporaneous ///////////////////////////////////////////////////
  // Theta
  to_vector(theta_sd_raw) ~ std_normal(); // prior on sigma_theta
  mu_theta_sd             ~ std_normal(); // prior on mu_theta_sd
  sigma_theta_sd          ~ std_normal(); // prior on sigma_theta_sd

  // Regression
  reg_intercept     ~ student_t(3, 0, 2); // prior on reg_intercept
  reg_slope_density ~ student_t(3, 0, 2); // prior on reg_slope_density
  reg_residual      ~ student_t(3, 0, 1); // prior on reg_residual
  
  
  int position_Y = 1; // position counter for the data
  for (i in 1:I) {
    // Precision Matrix prior
    L_Theta[i] ~ lkj_corr_cholesky(1); // prior on L_Theta
    
    //// Likelihood //////////////////////////////////////////////////////////
    
    // Partition data for one respondent
    array[n_t[i]] vector[K] Y_temp = segment(Y, position_Y, n_t[i]); // slice array
    position_Y += n_t[i]; // increment position counter
    
    // Cholesky decomposition of the covariance matrix
    matrix[K, K] Sigma_chol = diag_pre_multiply(0.5*exp(theta_sd[i]), L_Theta[i]);
    array[n_t[i]-1] vector[K] mu_network;
    
    // network predictions: loop over time points
    for(t in 1:(n_t[i]-1)){
      mu_network[t] = 
        to_vector(Intercepts[i, ]) + Beta[i] * (Y_temp[t,] - to_vector(Intercepts[i, ])); // predictions for the network
      } // end t
      
    // Network
    Y_temp[2:n_t[i]] ~ multi_normal_cholesky(mu_network, Sigma_chol);
  } // end i
  
  // Regression
  for(p in 1:P){
    reg_covariate[,p] ~ normal(mu_regression[,p], exp(reg_residual[p]));
  } // end p
    
} // end model
////////////////////////////////////////////////////////////////////////////////
generated quantities{
  vector[P]  reg_slope_density_z;
  vector[P] reg_intercept_z;
  
  // Standardize Regression Coefficients
  // beta * sd(x) = beta_std;
  reg_slope_density_z[1] = reg_slope_density[1] * sd(Beta_density);
  reg_slope_density_z[2] = reg_slope_density[2] * sd(Beta_density);
  reg_slope_density_z[3] = reg_slope_density[3] * sd(Beta_density);
  
  reg_slope_density_z[4] = reg_slope_density[4] * sd(Rho_density);
  reg_slope_density_z[5] = reg_slope_density[5] * sd(Rho_density);
  reg_slope_density_z[6] = reg_slope_density[6] * sd(Rho_density);
  
  reg_slope_density_z[7] = reg_slope_density[7] * sd(Beta_out_strength[,1]);
  reg_slope_density_z[8] = reg_slope_density[8] * sd(Beta_out_strength[,1]);
  reg_slope_density_z[9] = reg_slope_density[9] * sd(Beta_out_strength[,1]);
  
  // Standardize Regression Intercept
  // alpha_std = alpha + (beta_std * mean(x)) / sd(x);
  reg_intercept_z[1] = reg_intercept[1] + (reg_slope_density_z[1] * mean(Beta_density)) / sd(Beta_density);
  reg_intercept_z[2] = reg_intercept[2] + (reg_slope_density_z[2] * mean(Beta_density)) / sd(Beta_density);
  reg_intercept_z[3] = reg_intercept[3] + (reg_slope_density_z[3] * mean(Beta_density)) / sd(Beta_density);
  
  reg_intercept_z[4] = reg_intercept[4] + (reg_slope_density_z[4] * mean(Rho_density)) / sd(Rho_density);
  reg_intercept_z[5] = reg_intercept[5] + (reg_slope_density_z[5] * mean(Rho_density)) / sd(Rho_density);
  reg_intercept_z[6] = reg_intercept[6] + (reg_slope_density_z[6] * mean(Rho_density)) / sd(Rho_density);
  
  reg_intercept_z[7] = reg_intercept[7] + (reg_slope_density_z[7] * mean(Beta_out_strength[,1])) / sd(Beta_out_strength[,1]);
  reg_intercept_z[8] = reg_intercept[8] + (reg_slope_density_z[8] * mean(Beta_out_strength[,1])) / sd(Beta_out_strength[,1]);
  reg_intercept_z[9] = reg_intercept[9] + (reg_slope_density_z[9] * mean(Beta_out_strength[,1])) / sd(Beta_out_strength[,1]);
  
} // end generated quantities


  