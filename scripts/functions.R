# GVAR based ML sim function ----------------------------------------------
# argument sparse_sim: if TRUE, no random noise/random effects will be added
# to zero fixed effects
sim_gvar_loop <- function(graph,
                          beta_sd,
                          kappa_sd,
                          n_person,
                          n_time,
                          n_node,
                          max_try = 1000,
                          sparse_sim = FALSE,
                          failsafe = FALSE,
                          listify = FALSE) {
  # browser()
  
  # Create a list to store the data for each person
  data <- vector("list", n_person)
  beta <- array(NA, c(n_node, n_node, n_person))
  kappa <- array(NA, c(n_node, n_node, n_person))
  pcor <- array(NA, c(n_node, n_node, n_person))
  
  if(isTRUE(sparse_sim)){
    # helper function to identify zero elements by index
    find_zero_indices <- function(mat) {
      return(which(mat == 0, arr.ind = TRUE))
    }
    zeros_beta <- find_zero_indices(graph$beta)
    zeros_kappa <- find_zero_indices(graph$kappa)
  }
  
  
  # Simulate data for each person
  for (i in seq_len(n_person)) {
    counter <- 0
    beta_counter <- 0
    kappa_counter <- 0
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
      
      # Check if beta matrix is stable
      ev_b <- eigen(beta[, , i])$values
      if (all(Re(ev_b) ^ 2 + Im(ev_b) ^ 2 < 1))
        break
    }
    
    
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
        kappa[ , , i][zeros_kappa] <- 0
      }
      
      # Check if kappa matrix is semi-positive definite
      ev_k <- eigen(kappa[, , i])$values
      
      # add small tolerance 
      if (all(ev_k >= 0 - 1e-6)){
        pcor[, , i] <- -stats::cov2cor(kappa[, , i])
        diag(pcor[, , i]) <- 0
        if(!any(is.na(pcor[,,i]))){
          break
        }
        
        }
    }
    
    if (failsafe) {
      data[[i]] <-
        tryCatch({
          graphicalVAR::graphicalVARsim(nTime = n_time,
                                        beta = beta[, , i],
                                        kappa = kappa[, , i])
        }, error = function(e)
          NA)
    } else {
      repeat{
      counter <- counter + 1
      data[[i]] <- graphicalVAR::graphicalVARsim(nTime = n_time,
                                                 beta = beta[, , i],
                                                 kappa = kappa[, , i])
      if(!is.null(data[[i]])) break
      
      
    }
    }
  }
  ret <- list()
  if(isTRUE(listify)) {
  # Convert 3D arrays beta, kappa, pcor to list
  ret$beta_l <- lapply(1:n_person, function(i) beta[, , i])
  ret$kappa_l <- lapply(1:n_person, function(i) kappa[, , i])
  ret$pcor_l <- lapply(1:n_person, function(i) pcor[, , i])
  }
  
  # Return the list of simulated data
  # Also return the parameters used in the simulation
  ret$data <- data
  ret$beta <- beta
  ret$kappa <- kappa
  ret$pcor <- pcor
  return(ret)
}

# -------------------------------------------------------------------------
# mlVAR simulation function -----------------------------------------------
# -------------------------------------------------------------------------
# Copied from: https://github.com/SachaEpskamp/mlVAR/blob/master/R/mlVARmodel.R
# and adapted
# 1. fixed effects (means) for every parameter and a variance-covariance matrix.
# 2. Generate Beta's 
# 3. Shrink all beta's untill all eigens are in unit circle.
mlVARsim_mod <- function(
    # Simulation setup:
  nPerson = 10, # # of persons.
  nNode = 5, # # of nodes 
  nTime = 100, # Or vector with time points per person. 
  lag = 1,
  thetaVar = rep(1,nNode),
  DF_theta = nNode*2,
  mu_SD = c(1,1),
  init_beta_SD = c(0.1,1),
  fixedMuSD = 1,
  shrink_fixed = 0.9,
  shrink_deviation = 0.9,
  pcor_sparsity = 0.5,
  beta_sparsity = 0.5
){
  contemporaneous <- "wishart"
  GGMsparsity = pcor_sparsity
  
  if (length(nTime)==1){
    nTime <- rep(nTime,nPerson)
  }
  # Number of temporal effects:
  nTemporal <- nNode^2 * lag
  
  # 1. Generate structures:
  # Generate Omega:
  
  # Simulate mu means:
  Omega_mu <- cov2cor(solve(diag(nNode)-mlVAR:::simGraph(nNode)))
  # Omega_mu <- genPositiveDefMat(nNode, "onion", rangeVar = c(1,1))$Sigma
  # Simulate temporal:
  Omega_Beta <- cov2cor(solve(diag(nTemporal)-mlVAR:::simGraph(nTemporal)))
  
  #   mat <- Omega_mu
  #   diag(mat) <- 0
  #   rowSums(mat)
  #   
  Omega <- rbind(
    cbind(Omega_mu, matrix(0,nNode,nTemporal)),
    cbind(matrix(0,nTemporal,nNode), Omega_Beta)
  )
  
  # Omega <- genPositiveDefMat(nNode + nTemporal, "onion", rangeVar = c(1,1))$Sigma
  
  # Generate SD and scale:
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1],nNode),rep(init_beta_SD[1],nNode)), c(rep(mu_SD[2],nNode),rep(init_beta_SD[2],nNode)))
  Omega <- diag(SD) %*%Omega %*% diag(SD)
  
  # Generate fixed contemporaneous:
  if (contemporaneous=="wishart"){
    # Theta_fixed <- genPositiveDefMat(nNode, "onion", rangeVar = c(1,1))$Sigma
    Theta_fixed <- cov2cor(solve(diag(nNode)-mlVAR:::simGraph(nNode)))
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    
    # 2. Generate residual covariance matrices:
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else {
    
    if (contemporaneous == "randomGGM"){
      Theta <- lapply(1:nPerson,function(x){
        net <- mlVAR:::simGraph(nNode,GGMsparsity)
        cov2cor(solve(diag(nNode) - net))
      })
      Theta <- do.call(abind,c(Theta,along=3))
      
      Theta_fixed <- apply(Theta,1:2,mean)
    } else {
      net <- mlVAR:::simGraph(nNode,GGMsparsity)
      Theta_fixed <- cov2cor(solve(diag(nNode) - net))
      Theta <- lapply(1:nPerson,function(x)Theta_fixed)
      Theta <- do.call(abind,c(Theta,along=3))
    }
  }
  
  
  
  # Generate fixed means:
  mu_fixed <- rnorm(nNode,0,fixedMuSD)
  # Generate fixed betas:
  beta_fixed <- rnorm(nTemporal,0)
  # set weakest beta_sparsity% to zero:
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal * beta_sparsity)]] <- 0
  
  # Include auto-regressions:
  mat <- matrix(0,nNode,nNode*lag)
  diag(mat) <- 1
  beta_fixed[c(mat)==1] <- runif(sum(c(mat)==1),0,1)
  
  # 3. Generate random parameter sets:
  if (lag > 0){
    repeat{
      Pars <- rmvnorm(nPerson, c(mu_fixed,beta_fixed), sigma = Omega)
      Mus <- Pars[,1:nNode]
      
      Betas <- array(c(t(Pars[,-(1:nNode)])), c(nNode,nNode*lag,nPerson))
      
      # 4. Construct the matrices:
      if (lag>1){
        under <- cbind(diag(nNode*(lag-1)),matrix(0,nNode*(lag-1),nNode))
        
        ev <- sapply(seq_len(nPerson),function(i){
          mat <- rbind(Betas[,,i],under)
          eigen(mat)$values
        })
        
      } else {
        
        ev <- sapply(seq_len(nPerson),function(i){
          eigen(Betas[,,i])$values
        })
        
      }
      # 5. Store results:
      allEV <- c(ev)
      
      # 6. Break if all Re(ev)^2 + Im(ev)^2 < 1
      if (all(Re(ev)^2 + Im(ev)^2 < 1)){
        
        # simulate VAR for every person:
        DataList <- lapply(1:nPerson,function(p){
          
          pars <- lapply(seq_len(lag),function(l)array(c(Betas[,,p]),c(nNode,nNode,lag))[,,l])
          # If lag > 0 simulate VAR:
          if (lag > 0){
            res <- mlVAR::simulateVAR(pars, 
                                      means = Mus[p,], 
                                      lags = seq_len(lag), 
                                      Nt = nTime[p],
                                      init = Mus[p,],
                                      burnin = 100,
                                      residuals = Theta[,,p]) 
          } else {
            res <- rmvnorm(nTime[p],Mus[p,],Theta[,,p])
          }
          colnames(res) <- paste0("V",1:nNode)
          res$ID <- p
          res
        })
        
        # Rbind data:
        Data <- do.call(rbind,DataList)
        
        # 10. If any absolute > 100, go to 6a
        if (!any(abs(Data[,1:nNode]) > 100)){
          break
        }
      } 
      
      # Else shrink:
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode),-(1:nNode)] <- shrink_deviation * D[-(1:nNode),-(1:nNode)] 
      Omega <- D %*% cov2cor(Omega) %*% D
    }
    
  } else {
    Pars <- rmvnorm(nPerson, mu_fixed, sigma = Omega)
    Mus <- Pars[,1:nNode]
    Betas <- array(dim = c(0,0,nPerson))
    
    # simulate VAR for every person:
    DataList <- lapply(1:nPerson,function(p){
      res <- as.data.frame(rmvnorm(nTime[p],Mus[p,],Theta[,,p]))
      colnames(res) <- paste0("V",1:nNode)
      res$ID <- p
      res
    })
    
    # Rbind data:
    Data <- do.call(rbind,DataList)
    
  }
  
  
  # Create the list:
  model <- list(
    mu = mlVAR:::modelArray(mean = mu_fixed, 
                            SD = mu_SD, 
                            subject = lapply(1:nrow(Mus), function(i)Mus[i,])),
    Beta = mlVAR:::modelArray(mean = array(beta_fixed,c(nNode,nNode,lag)), 
                              SD = array(sqrt(diag(Omega[-(1:nNode),-(1:nNode)])),c(nNode,nNode,lag)), 
                              subject = lapply(1:nPerson, 
                                               function(p)
                                                 array(Betas[,,p],c(nNode,nNode,lag)))),
    Omega_mu = mlVAR:::modelCov(
      cov = mlVAR:::modelArray(mean = Omega[1:nNode,1:nNode])
    ),
    Theta = mlVAR:::modelCov(
      cov = mlVAR:::modelArray(mean = Theta_fixed, 
                               subject = lapply(1:nPerson,function(p)Theta[,,p]))
    ),
    Omega = mlVAR:::modelCov(
      cov = mlVAR:::modelArray(mean = Omega)
    )
    
  )
  
  
  # Data generated! Now return in sensible list:
  Results <- list(
    Data = Data,
    #     beta_fixed = array(beta_fixed,c(nNode,nNode,lag)),
    #     beta_SD = array(sqrt(diag(Omega[-(1:nNode),-(1:nNode)])),c(nNode,nNode,lag)),
    #     mu_fixed = mu_fixed,
    #     mu_SD = sqrt(diag(Omega[1:nNode,1:nNode])),
    #     Theta_fixed = Theta_fixed,
    #     DF_theta = DF_theta,
    #     Omega = Omega,
    #     Omega_mu = Omega[1:nNode,1:nNode],
    #     Omega_beta = Omega[-(1:nNode),-(1:nNode)],
    #     Theta = Theta,
    #     Mus = Mus,
    #     Betas = Betas,
    vars = paste0("V",1:nNode),
    idvar = "ID",
    lag=lag,
    model=model
  )
  
  class(Results) <- "mlVARsim"
  
  return(Results)
}






# -------------------------------------------------------------------------
# graphicalVAR helper functions -------------------------------------------
# -------------------------------------------------------------------------
# Function to extract centralities from graphicalVAR fit object
centrality_gvar <- function(fit){
  
  #--- Prepare matrices
  n_var <- ncol(fit$kappa)
  
  
  beta <- abs(fit$beta[,-1])
  pcor <- abs(-cov2cor(fit$kappa))
  diag(pcor) <- 0
  # pcor[lower.tri(pcor)] <- 0L
  
  #--- Outstrength 
  outstrength <- colSums(beta)/n_var
  
  #--- Instrength
  instrength <- rowSums(beta)/n_var
  
  #--- Strength
  strength <- colSums(pcor)/n_var
  
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
# Function to extract correlations from GIMME fit object with VAR = TRUE
gimme_cor_mat <- function(gimme_res, id, n_vars) {
  df_id <- gimme_res %>%
    filter(op == "~~") |>
    filter(pval < 0.05) |>
    mutate(file = str_remove(file, "subj")) %>%
    filter(file == id)
  
  corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  rownames(corr_matrix) <- paste0("V", 1:n_vars)
  colnames(corr_matrix) <- paste0("V", 1:n_vars)
  
  
  df_id <- df_id %>%
    select(lhs, rhs, beta.std) %>%
    spread(lhs, beta.std) %>%
    column_to_rownames(var = "rhs") %>%
    # replace all NAs with zero
    replace(is.na(.), 0) %>%
    as.matrix()
  
  # if there is no correlation
  if (nrow(df_id) == 0) {
    return(corr_matrix)
  }
  else{
    for (i in rownames(df_id)) {
      for (j in colnames(df_id)) {
        corr_matrix[i, j] <- df_id[i, j]
        corr_matrix[j, i] <- df_id[i, j]
      }
    }
    
    diag(corr_matrix) <- 0
    return(corr_matrix)
  }
  
}



# Function to extract centralities from GIMME fit object
centrality_gimme <- function(fit, var_only = FALSE){
  
  #--- Prepare 
  n_var <- fit$n_vars_total
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
  # if var_only = TRUE, first need to extract partial correlations
  # gimme somehow does not return them 
  

  #--- Centrality
  outstrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      colSums(abs(x[, temp_ind]))/n_var
    }
  })
  instrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      rowSums(abs(x[, temp_ind]))/n_var
    }
  })
  if(isFALSE(var_only)){
    strength <- lapply(fit$path_est_mats, function(x){
      if(!is.double(x)){
        NA
      }
      else{
        colSums(abs(x))/n_var
      }
    })
  }
  if(isTRUE(var_only)){
    strength <- lapply(fit$contemp_mat, function(x){
      if(!is.double(x)){
        NA
      }
      else{
        colSums(abs(x))/n_var
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
centrality_mlvar <- function(fit){
  
  #--- Prepare
  n_var <- length(fit$fit$var)
  n_id <- length(fit$IDs)

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
  outstrength <- lapply(l_beta, function(x){
    rowSums(abs(x))/n_var
  })
  instrength <- lapply(l_beta, function(x){
    colSums(abs(x))/n_var
  })
  strength <- lapply(l_pcor, function(x){
    colSums(abs(x))/n_var
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
                                 sim_fn = "sim_gvar_loop"){
  
  if(sim_fn == "mlvar_sim"){
    #--- Prepare
    n_id <- length(simobj$model$mu$subject)
    n_var <- length(simobj$vars)
    
    #--- Obtain networks
    l_beta <- lapply(1:n_id, function(i){
      simobj$model$Beta$subject[[i]][,,1]
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
  outstrength <- lapply(l_beta, function(x){
    rowSums(abs(x))/n_var
  })
  instrength <- lapply(l_beta, function(x){
    colSums(abs(x))/n_var  
  })
  strength <- lapply(l_pcor, function(x){
    colSums(abs(x))/n_var
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
  
}


# Double check with an ml_sim object
# beta_arr <- array(0, dim = c(6, 6,50))
# for(i in 1:50){
#   beta_arr[,,i]<- ml_sim$model$Beta$subject[[i]][,,1]
# }
# beta_arr_mean <- apply(beta_arr, c(1,2),mean)
# beta_mean <- ml_sim$model$Beta$mean[,,1]
# beta_arr_mean - beta_mean

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
theme_compare <- function(){
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
      plot.title = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.2), hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.15)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1), hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.y = ggplot2::unit(1.5, "lines")
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
    rstan::extract(object = fit, pars = par, permute = FALSE) %>%
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
extract_all_estimates <- function(fit, n_id, n_var) {
  # compute estimates
  beta_est <- extract_estimates(fit, "Beta") %>%
    map(function(x){
      est_vector2matrix(x, n_id, n_var)})
  pcor_est <-  extract_estimates(fit, "Rho_vec") %>%
    map(function(x){
      est_vector2vector(x, n_id, (n_var*(n_var-1))/2)})
  pcor_centrality_est <-
    extract_estimates(fit, "Rho_centrality") %>%
    map(function(x){
      est_vector2vector(x, n_id, n_var)})
  contdens_est <-  extract_estimates(fit, "Rho_density")
  tempdens_est <-  extract_estimates(fit, "Beta_density")
  outstrength_est <-  extract_estimates(fit, "Beta_out_strength") %>% 
    map(function(x){
      est_vector2vector(x, n_id, n_var)})
  regression_slope_est <-
    extract_estimates(fit, "reg_slope_density_z")
  regression_intercept_est <-
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
      est_vector[idx] %>% matrix(n_var, n_var, byrow = FALSE)
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