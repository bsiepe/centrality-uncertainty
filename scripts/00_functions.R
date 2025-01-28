# DGP Functions -----------------------------------------------------------

# Write to LaTeX ----------------------------------------------------------
# (using a function adapted from https://www.r-bloggers.com/2020/08/matrix-to-latex/):
arr_to_latex <- function(arr){
  rows <- apply(arr, MARGIN=1, paste, collapse = " & ")
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(cat(paste("\\begin{bmatrix}", matrix_string, "\\end{bmatrix}")))
}


# Create VAR and Innov Matrix ---------------------------------------------
create_var_matrix <- function(p, 
                              diag_val = 0.35, 
                              off_diag_val = 0.075, 
                              boost_factor = 1.5, 
                              sparse = FALSE, 
                              sparsity_proportion = 0.15) {
  mat <- matrix(off_diag_val, p, p)
  diag(mat) <- diag_val
  
  if (sparse) {
    # Get only off-diagonal indices of the matrix, exclude first row and column
    # if we include it, the boost factor can remove zeroes
    off_diag_indices <- which(row(mat) != col(mat), arr.ind = TRUE)
    off_diag_indices <- off_diag_indices[-which(off_diag_indices[, 1] == 1), ]
    off_diag_indices <- off_diag_indices[-which(off_diag_indices[, 2] == 1), ]
    
    
    # number of effects to be set to zero
    num_to_zero <- round(sparsity_proportion * p * p)
    
    # sample indices
    zero_indices <- off_diag_indices[sample(1:nrow(off_diag_indices), num_to_zero), ]
    
    for (i in 1:nrow(zero_indices)) {
      mat[zero_indices[i, 1], zero_indices[i, 2]] <- 0
    }
  }
  
  first_col_sum <- sum(mat[, 1])
  boost_amount <- boost_factor * first_col_sum - first_col_sum
  
  mat[1, ] <- mat[1, ] + boost_amount / p
  mat[, 1] <- mat[, 1] + boost_amount / p
  
  mat <- round(mat, digits = 3)
  
  ev_b <- eigen(mat)$values
  if (isFALSE(all(Re(ev_b) ^ 2 + Im(ev_b) ^ 2 < 1))) {
    stop("VAR matrix is not stationary")
  }
  
  
  return(mat)
}


create_cov_matrix <- function(p, 
                              diag_val = 1, 
                              off_diag_val = 0.1, 
                              boost_factor = 1.5, 
                              sparse = FALSE, 
                              sparsity_proportion = 0.15) {
  mat <- matrix(off_diag_val, p, p)
  diag(mat) <- diag_val
  
  if (sparse) {
    upper_tri_indices <- which(row(mat) < col(mat), arr.ind = TRUE)
    # exclude first row (or else the boost factor will remove zeroes)
    upper_tri_indices <- upper_tri_indices[-which(upper_tri_indices[, 1] == 1), ]
    
    num_to_zero <- round(sparsity_proportion * length(upper_tri_indices[, 1]))
    
    zero_indices <- upper_tri_indices[sample(1:nrow(upper_tri_indices), num_to_zero), ]
    
    # set off-diagonal elements and their symmetric counterparts to zero
    for (i in 1:nrow(zero_indices)) {
      mat[zero_indices[i, 1], zero_indices[i, 2]] <- 0
      mat[zero_indices[i, 2], zero_indices[i, 1]] <- 0 # ensure symmetry
    }
  }
  
  # Exclude diagonal values from the boost calculation
  first_col_sum <- sum(mat[-1, 1])
  boost_amount <- boost_factor * first_col_sum - first_col_sum
  
  mat[1, -1] <- mat[1, -1] + boost_amount / (p - 1)
  mat[-1, 1] <- mat[-1, 1] + boost_amount / (p - 1)
  
  mat <- round(mat, digits = 3)
  
  ev_b <- eigen(mat)$values
  if (isFALSE(all(ev_b > 0))) {
    stop("Covariance matrix is not positive definite")
  }
  
  return(mat)
}



# GVAR based ML sim function ----------------------------------------------
# argument sparse_sim: if TRUE, no random noise/random effects will be added
# to zero fixed effects

# innov_var_fixed_sigma: if TRUE, the innovation variance is not sampled randomly

sim_gvar_loop <- function(graph,
                          beta_sd = 1,
                          kappa_sd = 1,
                          sigma_sd = 1,
                          n_person,
                          n_time,
                          n_node,
                          # which package is used to simulate the data
                          # either graphicalVAR or mlVAR
                          sim_pkg = "graphicalVAR",
                          max_try = 1000,
                          sparse_sim = FALSE,
                          failsafe = FALSE,
                          # if listify is TRUE, convert 3D-Array to list
                          listify = FALSE,
                          # should there be a minimum difference between the most central
                          # and the second most central node in the temporal network
                          most_cent_diff_temp = FALSE,
                          # if most_cent_diff_temp is TRUE, this is the minimum difference
                          # as a proportion of the centrality of the most central node
                          most_cent_diff_temp_min = 0.1,
                          # same for the contemporaneous network
                          most_cent_diff_cont = FALSE,
                          most_cent_diff_cont_min = 0.1, 
                          # do we induce differences in the density of the temporal network
                          change_density = FALSE,
                          # if change_density is TRUE, these are the minimum and 
                          # maximum scaling factors
                          change_density_factors = c(0.75, 1.25),
                          verbose = FALSE,
                          # do we fix the innovation variance
                          # i.e., do we add interindividual differences here?
                          innov_var_fixed_sigma = FALSE
                          ) {

  # Create a list to store the data for each person
  data <- vector("list", n_person)
  beta <- array(NA, c(n_node, n_node, n_person))   # VAR
  kappa <- array(NA, c(n_node, n_node, n_person))  # precision
  sigma <- array(NA, c(n_node, n_node, n_person))  # covariance
  pcor <- array(NA, c(n_node, n_node, n_person))   # pcors
  
  if(isTRUE(sparse_sim)){
    # helper function to identify zero elements by index
    find_zero_indices <- function(mat) {
      return(which(mat == 0, arr.ind = TRUE))
    }
    zeros_beta <- find_zero_indices(graph$beta)
    zeros_kappa <- find_zero_indices(graph$kappa)
    if(!is.null(graph$sigma)){
      zeros_sigma <- find_zero_indices(graph$sigma)
    }
  }
  
  # if change_density is true, generate uniform scaling factors
  # across the n persons
  if(isTRUE(change_density)){
    density_factors <- runif(n_person, 
                             min = change_density_factors[1], 
                             max = change_density_factors[2])
  }
  
  
  # Simulate data for each person
  for (i in seq_len(n_person)) {
    counter <- 0
    beta_counter <- 0
    kappa_counter <- 0
    sigma_counter <- 0
    repeat {
      beta_counter <- beta_counter + 1
      if (beta_counter > max_try)
        stop("Exceeded maximum number of attempts to generate stable beta matrix.")
      
      # Generate beta matrix using graph$beta as mean
      beta[, , i] <- graph$beta + matrix(rnorm(n_node * n_node,
                                               mean = 0,
                                               sd = beta_sd),
                                         nrow = n_node,
                                         ncol = n_node)
      # if sparse matrix should be generated, set true zero effects to zero
      if(isTRUE(sparse_sim)){
        beta[,,i][zeros_beta] <- 0
      }
      
      # if change_density is TRUE, scale the beta matrix
      if(isTRUE(change_density)){
        beta[,,i] <- beta[,,i] * density_factors[i]
      }
      
      
      # if minimum centrality difference is required
      if(isTRUE(most_cent_diff_temp)){
        cent_check <- check_centrality_diff(matrix = beta[, , i],
                                            n_node = n_node,
                                            min_diff = most_cent_diff_temp_min)
        # if not, try again
        if(!isTRUE(cent_check)){
          if(isTRUE(verbose)){
            print(paste0("Centrality difference too small, trying again. Difference: ", 
                         max(centralities) - sort(centralities, decreasing = TRUE)[2],
                         "Individual:", i))
          }
          next
        }
        
      }
      
      
      
      # Check if beta matrix is stable
      # choose value slightly below 1 to avoid almost instable estimates
      # cf. Haslbeck, Epskamp & Waldorp (2023)
      ev_b <- eigen(beta[, , i])$values
      if (all(Re(ev_b) ^ 2 + Im(ev_b) ^ 2 < .95)){
        break
      }
      else {
        if(isTRUE(verbose)){
          print(paste0("Beta matrix not stable, trying again. Individual: ", i))
        }
      }
    } # end repeat statement
    
    # use kappa if no covariance matrix is provided
    if(is.null(graph$sigma)){
    repeat {
      kappa_counter <- counter + 1
      if (kappa_counter > max_try)
        stop(
          "Exceeded maximum number of attempts to generate semi-positive definite kappa matrix."
        )
      
      # Generate kappa matrix using graph$kappa as mean
      kappa[, , i] <-
        as.matrix(Matrix::forceSymmetric(graph$kappa + matrix(
          rnorm(n_node * n_node,
                mean = 0,
                sd = kappa_sd),
          nrow = n_node,
          ncol = n_node
        )))
      
      # if sparse matrix should be generated, set true zero effects to zero
      if(isTRUE(sparse_sim)){
        if(is.null(graph$sigma)){
          kappa[ , , i][zeros_kappa] <- 0
        } else {
          kappa[ , , i][zeros_sigma] <- 0
        }
      }
    
      
      # Check if kappa matrix is semi-positive definite
      ev_k <- eigen(kappa[, , i])$values
      
      # add small tolerance 
      if (all(ev_k >= 0 - 1e-6)){
        pcor[, , i] <- -stats::cov2cor(kappa[, , i])
        diag(pcor[, , i]) <- 0
        # same pcor effects are zero as in kappa/sigma
        if(isTRUE(sparse_sim)){
          if(is.null(graph$sigma)){
            pcor[ , , i][zeros_kappa] <- 0
          } else{
            pcor[ , , i][zeros_sigma] <- 0
          }
          
        }
        if(!any(is.na(pcor[,,i]))){
          break
        }
        
      } 
      else {
        if(isTRUE(verbose)){
          print("Kappa matrix not semi-positive definite, trying again.")
        }
      }
    # if minimum centrality difference is required
      if(isTRUE(most_cent_diff_cont)){
        cent_check <- check_centrality_diff(matrix = pcor[, , i],
                                            n_node = n_node,
                                            min_diff = most_cent_diff_cont_min)
        # if not, try again
        if(!isTRUE(cent_check)){
          if(isTRUE(verbose)){
            print(paste0("Centrality difference too small, trying again. Difference: ", 
                         max(centralities) - sort(centralities, decreasing = TRUE)[2],
                         "Individual:", i))
          }
          next
        }
        
      }
      
      
    } # end repeat statement
  } # end if(is.null(graph$sigma))  
    
    # if covariance matrix is provided
    if(!is.null(graph$sigma)){
      repeat{
        sigma_counter <- counter + 1
        if (sigma_counter > max_try)
          stop(
            "Exceeded maximum number of attempts to generate semi-positive definite sigma matrix."
          )
        
        if(isTRUE(innov_var_fixed_sigma )){
         noise_mat <- matrix(rnorm(n_node * n_node,
                                   mean = 0,
                                   sd = sigma_sd),
                             nrow = n_node,
                             ncol = n_node)
         noise_mat <- as.matrix(Matrix::forceSymmetric(noise_mat))
         diag(noise_mat) <- 0
         sigma[ , , i] <- as.matrix(Matrix::forceSymmetric(graph$sigma + noise_mat))
        }
        else{
          sigma[ , , i] <- as.matrix(Matrix::forceSymmetric(graph$sigma + matrix(
            rnorm(n_node * n_node,
                  mean = 0,
                  sd = sigma_sd),
            nrow = n_node,
            ncol = n_node
          )))
        }

        # if sparse matrix should be generated, set true zero effects to zero
        if(isTRUE(sparse_sim)){
          sigma[ , , i][zeros_sigma] <- 0
        }
        

        # Check if sigma matrix is semi-positive definite
        ev_s <- eigen(sigma[, , i])$values
        
        # add small tolerance
        if (all(ev_s >= 0 - 1e-6)){
          kappa[, , i] <- solve(sigma[, , i])
          pcor[, , i] <- -stats::cov2cor(kappa[, , i])
          diag(pcor[, , i]) <- 0
          if(!any(is.na(pcor[,,i]))){
            break
          }
          
        }
        else {
          if(isTRUE(verbose)){
            print("Sigma matrix not semi-positive definite, trying again.")
          }
        }
      # if minimum centrality difference is required
        if(isTRUE(most_cent_diff_cont)){
          cent_check <- check_centrality_diff(matrix = pcor[, , i],
                                              n_node = n_node,
                                              min_diff = most_cent_diff_cont_min)
          # if not, try again
          if(!isTRUE(cent_check)){
            if(isTRUE(verbose)){
              print(paste0("Centrality difference too small, trying again. Difference: ", 
                           max(centralities) - sort(centralities, decreasing = TRUE)[2],
                           "Individual:", i))
            }
            next
          }
          
        }  

      }  # end repeat statement
      

    } # end sigma generation
    
    
    
    if(sim_pkg == "graphicalVAR"){
      if (failsafe) {
        data[[i]] <-
          tryCatch({
            graphicalVAR::graphicalVARsim(nTime = n_time,
                                          beta = beta[, , i],
                                          kappa = kappa[, , i],
                                          warmup = 250)
          }, error = function(e)
            NA)
      } else {
        repeat{
          counter <- counter + 1
          data[[i]] <- graphicalVAR::graphicalVARsim(nTime = n_time,
                                                     beta = beta[, , i],
                                                     kappa = kappa[, , i],
                                                     warmup = 250)
          if(!is.null(data[[i]])) break
          
          
        }
      }
    }
    if(sim_pkg == "mlVAR"){
      if (failsafe) {
        data[[i]] <-
          tryCatch({
            if(!is.null(graph$sigma)){
              resid_cov = graph$sigma
            } else {
              resid_cov = solve(kappa[,,i])
            }
            
            mlVAR::simulateVAR(Nt = n_time,
                               pars = beta[, , i],
                               # means = rnorm(n = n_node, mean = 0, sd = 0.3),
                               residuals = resid_cov,
                               burnin = 250)
          }, error = function(e)
            NA)
      } else {
        repeat{
          counter <- counter + 1
          if(!is.null(graph$sigma)){
            resid_cov = graph$sigma
          } else {
            resid_cov = solve(kappa[,,i])
          }
          data[[i]] <- mlVAR::simulateVAR(Nt = n_time,
                                          pars = beta[, , i],
                                          # means = (n = n_node, mean = 0, sd = 0.3),
                                          # kappa = kappa[, , i],
                                          residuals = resid_cov,
                                          burnin = 250)
          if(!is.null(data[[i]])) break
          
          
        }
      }
    }

  }
  ret <- list()
  if(isTRUE(listify)) {
  # Convert 3D arrays beta, kappa, pcor to list
  ret$beta_l <- lapply(1:n_person, function(i) beta[, , i])
  ret$kappa_l <- lapply(1:n_person, function(i) kappa[, , i])
  ret$pcor_l <- lapply(1:n_person, function(i) pcor[, , i])
  if(!is.null(graph$sigma)){
    ret$sigma_l <- lapply(1:n_person, function(i) sigma[, , i])
  }
  }
  
  # Return the list of simulated data
  # Also return the parameters used in the simulation
  ret$data <- data
  ret$beta <- beta
  ret$kappa <- kappa
  ret$pcor <- pcor
  if(!is.null(graph$sigma)){
    ret$sigma <- sigma
  }
  return(ret)
}


# Helper function that checks for minimum required centrality difference
# between most central and second most central node
check_centrality_diff <- function(matrix, n_node, min_diff) {
  matrix_tmp <- matrix
  diag(matrix_tmp) <- 0
  centralities <- colSums(abs(matrix_tmp)) / (n_node - 1)
  diff_check <- max(centralities) - sort(centralities, decreasing = TRUE)[2] > min_diff * max(centralities)
  return(diff_check)
}






# -------------------------------------------------------------------------
# graphicalVAR helper functions -------------------------------------------
# -------------------------------------------------------------------------
# Function to extract centralities from graphicalVAR fit object
# can ignore autoregressive effects for centrality estimation
centrality_gvar <- function(fit,
                            ignore_ar = TRUE){  # should AR effects be ignored?
  
  #--- Prepare matrices
  n_var <- ncol(fit$kappa)
  n_var_temp <- n_var
  
  beta <- abs(fit$beta[,-1])
  beta_cent <- beta
  
  # if autoregressive effects should be ignored in centrality estimation
  if(isTRUE(ignore_ar)){
    diag(beta_cent) <- 0
    n_var_temp <- n_var - 1
  }
  
  pcor <- abs(-cov2cor(fit$kappa))
  diag(pcor) <- 0
  # pcor[lower.tri(pcor)] <- 0L
  
  #--- Outstrength 
  outstrength <- colSums(beta_cent)/n_var_temp
  
  #--- Instrength
  instrength <- rowSums(beta_cent)/n_var_temp
  
  #--- Strength
  # subtract 1 because diagonal is always zero
  strength <- colSums(pcor)/ (n_var - 1)
  
  #--- Temporal Density
  dens_temp <- sum(beta)/(n_var^2)
  
  #--- Contemporaneous Density
  # This should be equivalent to strength
  pcor[lower.tri(pcor)] <- 0L
  dens_cont <- sum(pcor)/(n_var * (n_var-1)/2)
  
  #--- Return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
}



# -------------------------------------------------------------------------
# GIMME helper functions --------------------------------------------------
# -------------------------------------------------------------------------
# Function to extract correlation matrix from GIMME fit object
gimme_cor_mat <- function(gimme_res, 
                              id, 
                              n_vars,
                              pcor = FALSE) {
  

  # obtain variable names
  var_names <- rownames(gimme_res$path_est_mats[[1]])  # Ensure this gets the correct names
  
  # filter input dataset
  df_id <- gimme_res$path_se_est |>
    filter(op == "~~") |>
    filter(pval < 0.05) |>
    mutate(file = str_remove(file, "subj")) |>
    filter(file == id)
  
  # initialize a zero correlation matrix with correct variable names
  corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  rownames(corr_matrix) <- var_names
  colnames(corr_matrix) <- var_names
  
  # early return if no significant correlations exist for this subject
  if (nrow(df_id) == 0) {
    return(corr_matrix)
  }
  
  # reshape `df_id` into a wide format matrix with `lhs` and `rhs` as row/column names
  df_id <- df_id |>
    select(lhs, rhs, beta.std) |>
    pivot_wider(names_from = lhs, values_from = beta.std) |>
    column_to_rownames(var = "rhs") |>
    mutate(across(everything(), ~ replace_na(., 0))) |>
    as.matrix()
  
  # all variables in `df_id` aligned with `corr_matrix`
  all_vars <- intersect(var_names, union(rownames(df_id), colnames(df_id)))
  for (i in all_vars) {
    for (j in all_vars) {
      # if a value exists in `df_id`, add it; otherwise leave it as 0
      if (i %in% rownames(df_id) && j %in% colnames(df_id)) {
        corr_matrix[i, j] <- df_id[i, j]
        corr_matrix[j, i] <- df_id[i, j]
      }
    }
  }
  
  diag(corr_matrix) <- 0
  pcor_matrix <- corr_matrix
  # If partial correlations are requested
  if (isTRUE(pcor)) {
    pcor_matrix <- try(corpcor::cor2pcor(corr_matrix), silent = TRUE)
    if (inherits(corr_matrix, "try-error") | any(is.na(pcor_matrix))) {
      warning("Partial correlation estimation failed. Returning correlation matrix.")
    }
  }
  
  
  
  return(pcor_matrix)
}




# Function to extract centralities from GIMME fit object
centrality_gimme <- function(fit, 
                             var_only = FALSE,
                             ignore_ar = TRUE){ # should AR effects be ignored?
  
  #--- Prepare 
  n_var <- fit$n_vars_total
  n_var_temp <- n_var
  n_id <- length(fit$data)
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
  
  
  #--- Density
  dens_temp <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      sum(abs(x[, temp_ind]))/(n_var^2)
    }
  })
  
  if(isTRUE(var_only)){
    dens_cont <- lapply(fit$contemp_mat, function(x){
      if(!is.double(x)){
        NA
      }
      else{
        diag(x) <- 0
        x[lower.tri(x)] <- 0L
        sum(abs(x))/ (n_var * (n_var-1)/2)
      }
    })
  }


  #--- Centrality
  # if autoregressive effects should be ignored in centrality estimation
  # diagonal of beta matrix is set to zero
  outstrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      beta_tmp <- x[, temp_ind]
      if(isTRUE(ignore_ar)){
        diag(beta_tmp) <- 0
        n_var_temp <- n_var - 1
      }
      colSums(abs(beta_tmp))/n_var_temp
    }
  })
  instrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      beta_tmp <- x[, temp_ind]
      if(isTRUE(ignore_ar)){
        diag(beta_tmp) <- 0
        n_var_temp <- n_var - 1
      }
      rowSums(abs(x[, temp_ind]))/n_var_temp
    }
  })
  if(isFALSE(var_only)){
    strength <- lapply(fit$path_est_mats, function(x){
      if(!is.double(x)){
        NA
      }
      else{
        # because we ignore the diagonal
        colSums(abs(x))/(n_var - 1)
      }
    })
  }
  if(isTRUE(var_only)){
    strength <- lapply(fit$contemp_mat, function(x){
      if(!is.double(x)){
        NA
      }
      else{
        # because we ignore the diagonal
        colSums(abs(x))/(n_var - 1)
      }
    })
  }
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
  
}

# -------------------------------------------------------------------------
# mlVAR helper functions --------------------------------------------------
# -------------------------------------------------------------------------
# Obtain centrality of fitted mlVAR object
centrality_mlvar <- function(fit,
                             ignore_ar = TRUE){
  
  #--- Prepare
  n_var <- length(fit$fit$var)
  n_id <- length(fit$IDs)
  n_var_temp <- n_var

  #--- Obtain networks
  l_beta <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "temporal",
                  nonsig = "show")
  })
  l_pcor <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit,
                  subject = i,
                  type = "contemporaneous",
                  nonsig = "show")
  })
  
  # Obtain overall adjacency matrix
  adj_beta <- ifelse(
    mlVAR::getNet(fit,
                  type = "temporal",
                  nonsig = "hide") != 0, 1, 0)
  adj_pcor <- ifelse(
    mlVAR::getNet(fit,
                  type = "contemporaneous",
                  nonsig = "hide") != 0, 1, 0) 
  
  # Set zero fixed effects to zero
  l_beta <- lapply(l_beta, function(x){
    x * adj_beta
  })
  l_pcor <- lapply(l_pcor, function(x){
    x * adj_pcor
  })
  
  
  #--- Density
  dens_temp <- lapply(l_beta, function(x){
    sum(abs(x))/(n_var^2)
  })
  dens_cont <- lapply(l_pcor, function(x){
    x <- x
    diag(x) <- 0
    x[lower.tri(x)] <- 0L
    sum(abs(x))/(n_var * (n_var-1)/2)
  })
  
  
  #--- Centrality
  # Important: in mlVAR, the lagged vars are rows, not columns
  # but only for the getNet function!
  outstrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    rowSums(abs(x))/n_var_temp
  })
  instrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    colSums(abs(x))/n_var_temp
  })
  strength <- lapply(l_pcor, function(x){
    # because diagonal is always zero, we subtract 1
    colSums(abs(x))/(n_var - 1)
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
}


# Obtain centrality of simulated mlVAR object
# works with mlvar_sim or sim_gvar_loop
centrality_mlvar_sim <- function(simobj,
                                 sim_fn = "sim_gvar_loop",
                                 ignore_ar = TRUE){
  
  if(sim_fn == "mlvar_sim"){
    #--- Prepare
    n_id <- length(simobj$model$mu$subject)
    n_var <- length(simobj$vars)
    n_var_temp <- n_var
    
    #--- Obtain networks
    l_beta <- lapply(1:n_id, function(i){
      # transpose to get the same format as graphicalVAR
      t(simobj$model$Beta$subject[[i]][,,1])
    })
    
    l_pcor <- lapply(1:n_id, function(i){
      x <- simobj$model$Theta$pcor$subject[[i]]
      diag(x) <- 0
      x
    })
    
  }
  
  if(sim_fn == "sim_gvar_loop"){
    #--- Prepare
    n_id <- dim(simobj$beta)[3]
    n_var <- dim(simobj$beta)[2]
    
    #--- Obtain networks
    l_beta <- lapply(1:n_id, function(i){
      simobj$beta[,,i]
    })
    
    l_pcor <- lapply(1:n_id, function(i){
      x <- simobj$pcor[,,i]
      diag(x) <- 0
      x
    })
    
  }

  #--- Density
  dens_temp <- lapply(l_beta, function(x){
    sum(abs(x))/(n_var^2)
  })
  dens_cont <- lapply(l_pcor, function(x){
    x <- x
    x[lower.tri(x)] <- 0L
    sum(abs(x))/(n_var * (n_var-1)/2)
  })
  
  #--- Centrality
  # Important: in mlVAR, the lagged vars are rows, not columns
  # this is only true for objects obtained with mlVAR::getNet
  outstrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    colSums(abs(x))/n_var_temp
  })
  instrength <- lapply(l_beta, function(x){
    if(isTRUE(ignore_ar)){
      diag(x) <- 0
      n_var_temp <- n_var - 1
    }
    rowSums(abs(x))/n_var_temp
  })
  strength <- lapply(l_pcor, function(x){
    # because the diagonal is zero, we subtract 1
    colSums(abs(x))/(n_var - 1)
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
  
}


#------------------------------------------------------------------------------>
# Simulation Helper Functions
#------------------------------------------------------------------------------>
abs_mean <- function(x){
  mean(abs(x), na.rm = TRUE)
}
abs_med <- function(x){
  stats::median(abs(x), na.rm = TRUE)
}

abs_sum <- function(x){
  sum(abs(x), na.rm = TRUE)
}

# Summary functions
# Compute mean pairwise RMSE for list of matrices
rmse_mean_list <- function(x, y){
  Map(function(x, y){
    mean(sqrt((x-y)^2), na.rm = TRUE)
  }, x, y)
}

mse_mean_list <- function(x, y){
  Map(function(x, y){
    mean((x-y)^2, na.rm = TRUE)
  }, x, y)
}


# Compute mean pairwise absolute bias for list of matrices
abs_mean_bias_list <- function(x, y){
  Map(function(x, y){
    abs_mean(x - y)
  }, x, y)
}

# Compute mean pairwise bias for list of matrices
bias_mean_list <- function(x, y){
  Map(function(x, y){
    mean(x - y, na.rm = TRUE)
  }, x, y)
}


# Find most central node and compare to true centrality
# for list of centralities
most_cent_ident <- function(x, y){
  Map(function(x, y){
    which.max(x) == which.max(y)
  }, x, y)
}

# Calculate MCSE for generic performance measures
mcse_generic <- function(x, n){
  sd(x, na.rm = TRUE)/sqrt(n)
}




#------------------------------------------------------------------------------>
# Helper functions to transform between different draws formats
#------------------------------------------------------------------------------>

# Helper function to transform draws matrix into a list of matrices for Beta or Sigma
draws_matrix2list <- function(draws_matrix) {
  iterations_list <-
    lapply(
      X = 1:nrow(draws_matrix),
      FUN = function(X) {
        matrix(draws_matrix[X,], ncol = sqrt(ncol(draws_matrix)), byrow = FALSE)
      }
    )
  return(iterations_list)
}

# Helper function to transform point estimates array of matrices into a list of matrices for Beta or Sigma
estimates_array2list <- function(draws_array) {
  est_list <-
    lapply(
      X = 1:dim(draws_array)[3],
      FUN = function(X) {
        draws_array[,,X]
      }
    )
  return(est_list)
}

# Helper function to transform draws array into a draws matrix
draws_matrix2array <- function(draws_matrix) {
  array <-
    array(t(draws_matrix),
          dim = c(sqrt(ncol(draws_matrix)),
                  sqrt(ncol(draws_matrix)),
                  nrow(draws_matrix)))
  
  return(array)
}

# Helper function to transform draws array into a draws matrix
# additionally deletes warmup samples specified manually
# or identified from BGGM object
draws_array2matrix <- function(array_3d,
                               warmup = 50) { # set to zero to keep everything
  iterations_list <-
    lapply(
      X = (warmup+1):dim(array_3d)[3],
      FUN = function(X) {
        as.vector(array_3d[, , X])
      }
    )
  matrix <- do.call(rbind, iterations_list)
  return(matrix)
}





# Compare fit to DGP ------------------------------------------------------
array_compare_dgp <- function(post_samples,
                              dgp = NULL,
                              plot = TRUE,
                              samples_pcor_name = "rho",
                              dgp_pcor_name = "pcor") {
  # Compute mean for each array element across the third dimension
  # of post_samples
  post_samples_mean <- lapply(post_samples, function(x) {
    apply(x, c(1, 2), mean)
  })
  post_samples_median <- lapply(post_samples, function(x) {
    apply(x, c(1, 2), median)
  })
  
  # Compare median of beta to DGP
  beta_diff <- post_samples_median$beta - dgp$beta
  rho_diff <-
    post_samples_median[[samples_pcor_name]] - dgp[[dgp_pcor_name]]
  
  
  result <- list(beta_diff = beta_diff, rho_diff = rho_diff)
  
  if (isTRUE(plot)) {
    # Plot both difference matrixes using cor.plot
    par(mfrow = c(1, 2))
    psych::cor.plot(beta_diff, main = "Beta difference")
    psych::cor.plot(rho_diff, main = "Rho difference", upper = FALSE)
  }
  
  return(result)
}




# ggplot theme ------------------------------------------------------------
# Create common theme for all plots
theme_centrality <- function(){
  # add google font
  sysfonts::font_add_google("News Cycle", "news")
  # use showtext
  showtext::showtext_auto()
  
  # theme
  ggplot2::theme_minimal(base_family = "news") +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(face = "plain",
                                         size = 22,
                                         hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 16,
                                            hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "plain", size = 18),
      axis.title.x = ggplot2::element_text(face = "plain", size = 18),
      axis.text.y = ggplot2::element_text(face = "plain", size = 18),
      axis.title.y = ggplot2::element_text(face = "plain", size = 18),
      axis.line = element_line(colour = "#6d6d6e"),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain",
                                         size = 20,
                                         hjust = 0.5),
      strip.text.x.top = ggplot2::element_text(face = "plain", 
                                               size = 20,
                                               hjust = 0.5),
      # strip.text.y = element_blank(),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.x = ggplot2::unit(1.25, "lines"),
      panel.spacing.y = ggplot2::unit(1.25, "lines"),
      # Remove vertical grid lines
      panel.grid.major.x = ggplot2::element_blank()
      
    )
}


okabe_fill_enh <- ggokabeito::palette_okabe_ito(order = c(5,1,3,4,2,6,7,8,9))
okabe_color_enh <- ggokabeito::palette_okabe_ito(order = c(5,1,3,4,2,6,7,8,9))


# ------------------------------------------------------------------------------
# Extract MCMC Draws from rstan fit  -------------------------------------------
# ------------------------------------------------------------------------------

# Extract estimates for one parameter
extract_estimates <- function(fit, par) {
  # extract draws
  draws <-
    rstan::extract(object = fit, pars = par, permute = FALSE) |>
    posterior::as_draws_matrix()
  # calculate estimates
  median <-
    try(bayestestR::point_estimate(draws, centrality = "median")[, "Median"])
  # use quantile-based due to potential problems with HDI
  ci_95_l <- try(bayestestR::ci(draws, ci = .95)[, "CI_low"])
  ci_95_u <- try(bayestestR::ci(draws, ci = .95)[, "CI_high"])
  ci_oneside <- try(bayestestR::ci(draws, ci = .90)[, "CI_low"])
  
  # split into matrices and return list of lists
  ret_list <- list(
    median = median,
    ci_95_l = ci_95_l,
    ci_95_u = ci_95_u,
    ci_oneside = ci_oneside
  )
  # return list
  return(ret_list)
}


# Extract estimates for all parameters
# if transpose_beta is TRUE, the beta estimates are transposed
# which means that they then have the lagged variables as columns
# if pcor_matrix is TRUE, pcor estimates will be returned as a matrix
# instead of as a vector
extract_all_estimates <- function(fit, 
                                  n_id, 
                                  n_var, 
                                  transpose_beta = FALSE,
                                  pcor_matrix = FALSE) {
  
  # compute estimates
  beta_est <- extract_estimates(fit, "Beta") |>
    map(function(x) {
      est_vector2matrix(x, n_id, n_var)
    })
  if(isTRUE(transpose_beta)){
    beta_est <- map(beta_est, ~map(.x, t))
  }
  
  pcor_est <-  extract_estimates(fit, "Rho_vec") |>
    map(function(x) {
      est_vector2vector(x, n_id, (n_var * (n_var - 1)) / 2)
    })
  
  # convert to matrix format
  if(isTRUE(pcor_matrix)){
    pcor_est <- lapply(pcor_est, function(x){
      # loop over medians and CIs
      lapply(x, function(y){
        pcor_mat <- matrix(data = 0, nrow = n_var, ncol = n_var)
        pcor_mat[upper.tri(pcor_mat)] <- as.numeric(y)
        pcor_mat <- pcor_mat + t(pcor_mat)
        pcor_mat
      })
    })
  }
  
  pcor_centrality_est <-
    extract_estimates(fit, "Rho_strength") |>
    map(function(x) {
      est_vector2vector(x, n_id, n_var)
    })
  
  contdens_est <-  extract_estimates(fit, "Rho_density")
  tempdens_est <-  extract_estimates(fit, "Beta_density")
  outstrength_est <-  extract_estimates(fit, "Beta_out_strength") |>
    map(function(x) {
      est_vector2vector(x, n_id, n_var)
    })
  instrength_est <-  extract_estimates(fit, "Beta_in_strength") |>
    map(function(x) {
      est_vector2vector(x, n_id, n_var)
    })
  regression_slope_est <-
    extract_estimates(fit, "reg_slope_density")
  regression_intercept_est <-
    extract_estimates(fit, "reg_intercept")
  regression_slope_est_z <-
    extract_estimates(fit, "reg_slope_density_z")
  regression_intercept_est_z <-
    extract_estimates(fit, "reg_intercept_z")
  
  # return list of lists
  ret_list <-
    list(
      beta_est = beta_est,
      pcor_est = pcor_est,
      pcor_centrality_est = pcor_centrality_est,
      contdens_est = contdens_est,
      tempdens_est = tempdens_est,
      outstrength_est = outstrength_est,
      instrength_est = instrength_est,
      regression_slope_est = regression_slope_est,
      regression_intercept_est = regression_intercept_est
    )
  
  return(ret_list)
}


### Helper functions to transform draws to a suitable format -------------------

# split vector into list of matrices (for beta)
est_vector2matrix <- function(est_vector, n_id, n_var) {
  # initialize list of matrices for each individual
  est_list <- lapply(1:n_id, function(x) matrix(NA, n_var, n_var))
  
  for (i in 1:n_id) {
    # create indicator for individual draws
    idx <- seq(i, length(est_vector), by = n_id)
    est_list[[i]] <-
      est_vector[idx] |> matrix(n_var, n_var, byrow = FALSE)
  }
  return(est_list)
}

# split vector into list of vectors (for pcor and pcor centrality)
est_vector2vector <- function(est_vector, n_id, n_var) {
  # initialize list of matrices for each individual
  est_list <- lapply(1:n_id, function(x) rep(NA, n_var))
  
  for (i in 1:n_id) {
    idx <- seq(i, length(est_vector), by = n_id)
    est_list[[i]] <-
      est_vector[idx]
  }
  return(est_list)
}