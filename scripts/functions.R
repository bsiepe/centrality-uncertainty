# -------------------------------------------------------------------------
# graphicalVAR helper functions -------------------------------------------
# -------------------------------------------------------------------------
# Function to extract centralities from graphicalVAR fit object
centrality_graphicalvar <- function(fit){
  
  #--- Prepare matrices
  beta <- abs(fit$beta[,-1])
  pcor <- abs(-cov2cor(fit$kappa))
  diag(pcor) <- 0
  pcor[lower.tri(pcor)] <- 0L
  
  #--- Outstrength 
  outstrength <- colSums(beta)
  
  #--- Instrength
  instrength <- rowSums(beta)
  
  #--- Strength
  strength <- colSums(pcor)
  
  #--- Temporal Density
  dens_temp <- sum(beta)
  
  #--- Contemporaneous Density
  # This should be equivalent to strength
  dens_cont <- sum(pcor)
  
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
# Function to extract centralities from GIMME fit object
centrality_gimme <- function(fit){
  
  
  #--- Density
  dens_temp <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      sum(abs(x[, temp_ind]))
    }
  })
  
  dens_cont <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      diag(x) <- 0
      x[lower.tri(x)] <- 0L
      sum(abs(x))
    }
  })

  #--- Centrality
  outstrength <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      colSums(abs(x))
    }
  })
  instrength <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      rowSums(abs(x))
    }
  })
  strength <- lapply(ref_path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      diag(x) <- 0
      x[lower.tri(x)] <- 0L
      colSum(abs(x))
    }
  })
  
  
  
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

#------------------------------------------------------------------------------>
# Function to extract the log_lik for the gVAR
#------------------------------------------------------------------------------>
log_lik_gVAR <- function(Y, draws_beta, draws_sigma, n_cores = 1) {
  # # save chain ids
  # chain_ids <- draws_beta %>%
  #   as_draws_df() %>%
  #   select(.chain) %>% unlist()
  
  # prepare matrices from draws
  Beta <- draws_matrix2list(draws_beta)
  Sigma <- draws_matrix2list(draws_sigma)
  # number of iterations
  n_iter <- length(Beta)
  n_t <- nrow(Y)
  # parallelization
  if (n_cores > 1) {
    future::plan("multisession", workers = n_cores)
  } else{
    future::plan("sequential")
  }
  future::plan("multisession", workers = n_cores)
  # progress bar
  p <- progressr::progressor(along = 1:n_iter)
  # loop over iteraions
  progressr::with_progress(log_lik_list <- furrr::future_map(
    .x = 1:n_iter,
    .f = function(n) {
      # loop over time points
      log_lik_row <-
        lapply(2:n_t, function(t) {
          return(mvtnorm::dmvnorm(
            x = Y[t,],
            mean = Beta[[n]] %*% Y[t - 1,],
            sigma = Sigma[[n]],
            log = TRUE
          ))
        }) # end timepoint
      log_lik_row <- do.call(cbind, log_lik_row)
      return(log_lik_row)
      p()
    }
    # end iterations
    # end progress bar
  ))
  # end paralleliztion
  
  if (n_cores > 1) {
    future::plan("sequential")
  }
  log_lik_list <- do.call(rbind, log_lik_list)
  
  log_lik_mat <- posterior::as_draws_matrix(log_lik_list)
  return(log_lik_mat)
}

#------------------------------------------------------------------------------>
# Function to fit the gVAR Model in Stan
#------------------------------------------------------------------------------>
fit_gVAR_stan <-
  function(data,
           # a vector of beeps with length of nrow(data)
           beep = NULL,
           priors = NULL,
           backend = "rstan",
           method = "sampling",
           cov_prior = "LKJ",
           # c("LKJ", "IW")
           rmv_overnight = FALSE,
           iter_sampling = 500,
           iter_warmup = 500,
           n_chains = 4,
           n_cores = 4,
           server = FALSE,  # temporary option to run code on Linux server
           center_only = FALSE,   # only center (not scale)
           ...) {
    if(isTRUE(center_only)){
      Y <- data %>% apply(., 2, scale, center = TRUE, scale = FALSE)
    } else{
      Y <- data %>% apply(., 2, scale, center = TRUE, scale = TRUE)
    }
    
    K <- ncol(data)
    n_t <- nrow(data)
    
    # Specify Priors
    if (is.null(priors)) {
      prior_Rho_loc <- matrix(.5, nrow = K, ncol = K)
      prior_Rho_scale <- matrix(.4, nrow = K, ncol = K)
      prior_Beta_loc <- matrix(0, nrow = K, ncol = K)
      prior_Beta_scale <- matrix(.5, nrow = K, ncol = K)
    } else{
      prior_Rho_loc <- priors[["prior_Rho_loc"]]
      prior_Rho_scale <- priors[["prior_Rho_scale"]]
      prior_Beta_loc <- priors[["prior_Beta_loc"]]
      prior_Beta_scale <- priors[["prior_Beta_scale"]]
    }
    
    # Stan Data
    stan_data <- list(
      K = K,
      "T" = n_t,
      Y = as.matrix(Y),
      beep = beep,
      prior_Rho_loc = prior_Rho_loc,
      prior_Rho_scale = prior_Rho_scale,
      prior_Beta_loc = prior_Beta_loc,
      prior_Beta_scale = prior_Beta_scale
    )
    
    # Choose model to fit
    if (cov_prior == "LKJ") {
      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_LKJ_beep"
      } else{
        # standard model
        model_name <- "VAR_LKJ"
      }
    }
    if (cov_prior == "IW") {
      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_wishart_beep"
      } else{
        # standard model
        model_name <- "VAR_wishart"
      }
    }
    
    
    if (backend == "rstan") {
      # Compile model
      if(isFALSE(server)){
        stan_model <-
          rstan::stan_model(file = here::here("scripts", paste0(model_name, ".stan")))        
      } else {
        stan_model <-
          rstan::stan_model(file = paste0("~/stan-gvar/scripts/", model_name, ".stan"))
      }
      
      
      if (method == "sampling") {
        # Run sampler
        stan_fit <- rstan::sampling(
          object = stan_model,
          data = stan_data,
          #pars = c("Beta_raw"),
          #include = FALSE,
          chains = n_chains,
          cores = n_cores,
          iter = iter_sampling + iter_warmup,
          warmup = iter_warmup,
          refresh = 500,
          thin = 1,
          init = .1,
          control = list(adapt_delta = .8),
          ...
        )
      }
      if (method == "variational") {
        stan_fit <- rstan::vb(
          object = stan_model,
          data = stan_data,
          #pars = c("Beta_raw"),
          #include = FALSE,
          init = .1,
          tol_rel_obj = .001,
          output_samples = iter_sampling * n_chains,
          ...
        )
      }
    } else{
      # Compile model
      if(isFALSE(server)){
      stan_model <-
        cmdstanr::cmdstan_model(stan_file = here::here("scripts", paste0(model_name, ".stan")),
                                pedantic = TRUE)
      } else {
        stan_model <-
          cmdstanr::cmdstan_model(file = paste0("~/stan-gvar/scripts/", model_name, ".stan"),
                                  pedantic = TRUE)        
      }
      if (method == "sampling") {
        # Run sampler
        stan_fit <- stan_model$sample(
          data = stan_data,
          chains = n_chains,
          parallel_chains = n_cores,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          refresh = 500,
          thin = 1,
          adapt_delta = .8,
          init = .1,
          ...
        )
      }
      if (method == "variational") {
        stan_fit <- stan_model$variational(
          data = stan_data,
          tol_rel_obj = .001,
          init = .1,
          output_samples = iter_sampling * n_chains,
          ...
        )
      }
    }
    return(stan_fit)
  }

#------------------------------------------------------------------------------>
# Function to Compute the LOO-CV for Independent Stan Model and Data
#------------------------------------------------------------------------------>
loo_gVAR <- function(stan_fit, data, n_cores = 1) {
  c <- class(stan_fit)
  if (attr(c, "package") == "rstan") {
    log_lik <-
      log_lik_gVAR(
        Y = data %>% apply(., 2, scale),
        draws_beta = posterior::as_draws_matrix(rstan::extract(
          stan_fit, pars = "Beta", permuted = FALSE
        )),
        draws_sigma = posterior::as_draws_matrix(rstan::extract(
          stan_fit, pars = "Sigma", permuted = FALSE
        )),
        n_cores = n_cores
      )
    chain_ids <-
      rstan::extract(stan_fit, pars = "Beta", permuted = FALSE) %>%
      posterior::as_draws_df() %>%
      dplyr::select(.chain) %>%
      unlist()
  } else{
    log_lik <-
      log_lik_gVAR(
        Y = data %>% apply(., 2, scale),
        draws_beta = posterior::as_draws_matrix(stan_fit$draws("Beta")),
        draws_sigma = posterior::as_draws_matrix(stan_fit$draws("Sigma")),
        n_cores = n_cores
      )
    chain_ids <- stan_fit$draws("Beta") %>%
      posterior::as_draws_df() %>%
      dplyr::select(.chain) %>%
      unlist()
  }
  loo <- loo::loo(log_lik, r_eff = relative_eff(log_lik, chain_ids))
  return(loo)
}

#------------------------------------------------------------------------------>
# Function to Compute the Bayes Factor for a Posterior Difference Matrix
#------------------------------------------------------------------------------>
# Helpers
fisher_z <- function(r) {
  return(0.5 * log((1 + r) / (1 - r)))
}
fisher_z_inv <- function(z) {
  return((exp(2 * z) - 1) / (exp(2 * z) + 1))
}








# -------------------------------------------------------------------------
# Helper functions for model evaluation -----------------------------------
# -------------------------------------------------------------------------
# Convert Stan fit to array -----------------------------------------------
# TODO should maybe be flexible to incorporate something else besides Sigma?
# i.e. also theta (precision matrix?)

stan_fit_convert <-
  function(stan_fit) {
    # check fitting backend
    c <- class(stan_fit)
    
    if (attr(c, "package") == "rstan") {
      draws_beta <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Beta", permuted = FALSE))
      draws_sigma <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Sigma", permuted = FALSE))
      draws_rho <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Rho", permuted = FALSE))
    }
    else{
      draws_beta <- posterior::as_draws_matrix(stan_fit$draws("Beta"))
      draws_sigma <-
        posterior::as_draws_matrix(stan_fit$draws("Sigma"))
      draws_rho <- posterior::as_draws_matrix(stan_fit$draws("Rho"))
    }
    # Convert to array of p x p matrices
    nvar <- sqrt(ncol(draws_beta))
    
    # Beta
    split_beta <- split(draws_beta, seq(nrow(draws_beta)))
    beta_l <- lapply(split_beta, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    beta_array <-
      array(unlist(beta_l), dim = c(nvar, nvar, nrow(draws_beta)))
    
    # Sigma
    split_sigma <- split(draws_sigma, seq(nrow(draws_sigma)))
    sigma_l <- lapply(split_sigma, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    sigma_array <-
      array(unlist(sigma_l), dim = c(nvar, nvar, nrow(draws_sigma)))
    
    # Rho
    split_rho <- split(draws_rho, seq(nrow(draws_rho)))
    rho_l <- lapply(split_rho, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    rho_array <-
      array(unlist(rho_l), dim = c(nvar, nvar, nrow(draws_rho)))
    
    # Return
    return(list(
      beta = beta_array,
      sigma = sigma_array,
      rho = rho_array
    ))
    
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
