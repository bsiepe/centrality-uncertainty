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



compare_matrices <-
  function(mat1,
           mat2,
           # number of bootstrap samples
           bootstrap_samples = 0,
           # c("Beta", "Rho")
           parameter_type = "Beta",
           # must be a list containing 2 matrices
           H1_custom_priors = NULL,
           # normal SD
           H1_prior_scale = 0.5,
           # "normal" = SD, "uniform" = range / 2
           H0_prior_scale = NULL,
           # c("normal", "uniform", "posterior_uncertainty")
           H0_distribution = "uniform",
           # plot densities for null, posterior, and prior
           plot = TRUE,
           # upper limit of the x-axis
           plot_xlim = NULL,
           # upper limit of the y-axis
           plot_ylim = NULL) {
    
    # fisher z-transform partial correlations
    if (parameter_type == "Rho") {
      mat1 <- draws_matrix2list(mat1)
      mat2 <- draws_matrix2list(mat2)
      
      mat1 <-
        purrr::map(mat1, function(x) {
          fisher_z(x) %>% .[lower.tri(.)] %>% as.vector()
        })
      mat2 <-
        purrr::map(mat2, function(x) {
          fisher_z(x) %>% .[lower.tri(.)] %>% as.vector()
        })
      
      mat1 <- do.call(rbind, mat1)
      mat2 <- do.call(rbind, mat2)
    }
    
    #### with bootstrapping ###
    if (bootstrap_samples > 0) {
      # sample rows from posterior iterations
      rows1 <-
        sample(1:nrow(mat1),
               size = bootstrap_samples,
               replace = TRUE)
      rows2 <-
        sample(1:nrow(mat2),
               size = bootstrap_samples,
               replace = TRUE)
      # assemble new posterior draws_matrices
      mat1 <-
        purrr::map(
          .x = rows1,
          .f = function(.x) {
            mat1[.x,]
          }
        ) %>%
        do.call(rbind, .)
      mat2 <-
        purrr::map(
          .x = rows2,
          .f = function(.x) {
            mat2[.x, ]
          }
        ) %>%
        do.call(rbind, .)
      # compute differences between matrices
      diff_post_mat <- abs(mat1 - mat2)
      diff_post <- diff_post_mat %>%
        apply(., 1, sum)
      
      
      ### compute prior difference matrices
      
      # use custom priors if exist
      if (typeof(H1_custom_priors) == "list" &
          length(H1_custom_priors) == 2) {
        loc <- H1_custom_priors[[1]] %>% as.vector()
        scale <- H1_custom_priors[[2]] %>% as.vector()
        diff_prior <- abs(
          lapply(1:length(loc), function(n) {
            rnorm(n = bootstrap_samples,
                  mean = loc[n],
                  sd = scale[n])
          }) %>%
            do.call(cbind, .) -
            lapply(1:length(loc), function(n) {
              rnorm(n = bootstrap_samples,
                    mean = loc[n],
                    sd = scale[n])
            }) %>%
            do.call(cbind, .)
        ) %>%
          matrix(., ncol = ncol(mat1)) %>%
          apply(., 1, sum)
        
        # TODO need potential warning message if priors are not
        # NULL but do not satisfy correct format
        
      } else{
        # use default priors
        diff_prior_mat <-
          abs(
            rnorm(
              n = bootstrap_samples * ncol(mat1),
              mean = 0,
              sd = H1_prior_scale
            ) -
              rnorm(
                n = bootstrap_samples * ncol(mat1),
                mean = 0,
                sd = H1_prior_scale
              )
          ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
      }
    } else{
      ### NO bootstrappping ###
      
      diff_post_mat <- mat1 - mat2
      attr(diff_post_mat, which = "class") <- "matrix"
      
      diff_post_mat <- abs(diff_post_mat)
      
      diff_post <- diff_post_mat %>%
        apply(., 1, sum)
      
      # prior difference matrices
      if (typeof(H1_custom_priors) == "list" &
          length(H1_custom_priors) == 2) {
        loc <- H1_custom_priors[[1]] %>% as.vector()
        scale <- H1_custom_priors[[2]] %>% as.vector()
        
        diff_prior_mat <- abs(
          lapply(1:length(loc), function(n) {
            rnorm(n = max(4e4, dim(mat1)[1]),
                  mean = loc[n],
                  sd = scale[n])
          }) %>%
            do.call(cbind, .) -
            lapply(1:length(loc), function(n) {
              rnorm(n = max(4e4, dim(mat1)[1]),
                    mean = loc[n],
                    sd = scale[n])
            }) %>%
            do.call(cbind, .)
        ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
        # use default priors
      } else{
        diff_prior_mat <- abs(
          rnorm(
            n = max(4e4, dim(mat1)[1]) * ncol(mat1),
            mean = 0,
            sd = H1_prior_scale
          ) -
            rnorm(
              n = max(4e4, dim(mat1)[1]) * ncol(mat1),
              mean = 0,
              sd = H1_prior_scale
            )
        ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
      }
    }
    
    # distribution for null rope range
    if (H0_distribution == "normal") {
      if (is.null(H0_prior_scale)) {
        H0_prior_scale <- .005
      }
      # diff_null_mat <- abs(
      #   rnorm(
      #     n = 4e4 * ncol(mat1),
      #     mean = 0,
      #     sd = H0_prior_scale
      #   ) -
      #     rnorm(
      #       n = 4e4 * ncol(mat1),
      #       mean = 0,
      #       sd = H0_prior_scale
      #     )
      # ) %>%
      #   matrix(., ncol = ncol(mat1))
      # diff_null <-  diff_null_mat %>%
      #   apply(., 1, sum)
      
      # TODO What does this do?
      null_lim <- H0_prior_scale * 4 * ncol(mat1)
    }
    if (H0_distribution == "uniform") {
      if (is.null(H0_prior_scale)) {
        H0_prior_scale <- .01
      }
      
      # diff_null_mat <- abs(
      #   runif(
      #     n = 4e4 * ncol(mat1),
      #     min = -H0_prior_scale,
      #     max = H0_prior_scale
      #   ) -
      #     runif(
      #       n = 4e4 * ncol(mat1),
      #       min = -H0_prior_scale,
      #       max = H0_prior_scale
      #     )
      # ) %>%
      #   matrix(., ncol = ncol(mat1))
      # diff_null <-  diff_null_mat %>%
      #   apply(., 1, sum)
      
      # TODO why twice actually? it is a one-sided test
      null_lim <- H0_prior_scale * 2 * ncol(mat1)
    }
    if (H0_distribution == "posterior_uncertainty") {
      # TODO Combine uncertainty from both networks?
      # sample rows from posterior iterations
      rows1_null <-
        sample(1:nrow(mat1),
               size = nrow(mat1),
               replace = FALSE)
      
      # assemble new posterior draws_matrices
      mat1_null <-
        purrr::map(
          .x = rows1_null,
          .f = function(.x) {
            mat1[.x, ]
          }
        ) %>%
        do.call(rbind, .)
      # compute differences between matrices
      diff_null <- abs(mat1 - mat1_null) %>%
        apply(., 1, sum)
      # HDI for Null distribution, CI_high will be used as upper ROPE limit
      null_lim <- bayestestR::hdi(diff_null, ci = .99)$CI_high
    }
    
    # TODO remove percentage framing here, then adjust that in simulation
    post_in_rope <-
      diff_post[diff_post < null_lim] %>%
      length() / length(diff_post) * 100
    
    prior_in_rope <-
      diff_post[diff_prior < null_lim] %>%
      length() / length(diff_prior) * 100
    
    ### Compute BF for H0
    prec <- 100 # precision for mpfr (floating point)
    
    mean_post_u <- Rmpfr::mpfr(mean(diff_post),precBits = prec)
    sd_post_u <- Rmpfr::mpfr(sd(diff_post),precBits = prec)
    mean_prior_u <- Rmpfr::mpfr(mean(diff_prior),precBits = prec)
    sd_prior_u <- Rmpfr::mpfr(sd(diff_prior),precBits = prec)
    
    # TODO needs to be described
    q <- Rmpfr::mpfr(null_lim,precBits = prec)
    .mpfr_erange_set(value = (1-2^-52)*.mpfr_erange(c("min.emin","max.emax")))
    
    log_BF_01_mpfr <- 
      Rmpfr::pnorm(
        q,
        mean_post_u,
        sd_post_u,
        log.p = TRUE,
        lower.tail = TRUE
        ) -
      Rmpfr::pnorm(
        q,
        mean_prior_u,
        sd_prior_u,
        log.p = TRUE,
        lower.tail = TRUE
      )
    log_BF_01 <- Rmpfr::asNumeric(log_BF_01_mpfr)
    BF_01 <- Rmpfr::asNumeric(exp(log_BF_01_mpfr))
    
    
    # plot prior vs. posterior
    if (isTRUE(plot)) {
      # combine diff_post and diff_prior in a single column of a dataframe
      # with another column as an indicator if it is a posterior or prior
      # distribution
      if (is.null(plot_xlim)) {
        plot_xlim <- max(median(diff_prior),
                         median(diff_post))#,
                         #median(diff_null))
        
      }
      if (is.null(plot_ylim)) {
        plot_ylim <- max(max(density.default(diff_post)[["y"]]),
                         max(density.default(diff_prior)[["y"]])) * 1.5
        
      }
      df_samples <- data.frame(posterior = diff_post,
                               prior = diff_prior) %>%
        tidyr::pivot_longer(
          cols = c(posterior, prior),
          names_to = "distribution",
          values_to = "value"
        ) #%>%
      # rbind(data.frame(distribution = "null",
      #                  value = diff_null))
      
      plot <- df_samples %>%
        ggplot2::ggplot() +
        ggplot2::geom_density(aes(value, fill = distribution),
                              alpha = .5) +
        ggplot2::annotate(
          'rect',
          xmin = 0,
          xmax = null_lim,
          ymin = 0,
          ymax = plot_ylim,
          alpha = .4,
          fill = 'grey'
        ) +
        ggplot2::geom_vline(xintercept = null_lim) +
        ggplot2::scale_x_continuous(limits = c(0, plot_xlim),
                                    expand = expansion()) +
        ggplot2::scale_y_continuous(limits = c(0, plot_ylim),
                                    expand = expansion()) +
        ggplot2::labs(x = "Difference", y = "Density") +
        ggplot2::theme_light() +
        ggplot2::theme(panel.grid = element_blank())
      
      print(plot)
    }
    
    # df with log_BF and overlap coef
    df_results <- data.frame(
      round(BF_01, 4),
      round(log_BF_01, 2),
      round(post_in_rope, 2),
      round(prior_in_rope, 2),
      row.names = NULL
    )
    names(df_results) <-
      c("BF_01",
        "log BF_01",
        "% posterior < null CI_high",
        "% prior < null CI_high")
    
    print(df_results)
    
    # return list with results silently
    return(invisible(
      list(
        BF_01 = BF_01,
        log_BF_01 = log_BF_01,
        post_below_null_ub = post_in_rope,
        prior_below_null_ub = prior_in_rope
      )
    ))
    
    # clean up
    rm(mat1,
       mat2,
       diff_post,
       diff_post_mat,
       diff_prior,
       diff_prior_mat,
       diff_null)
  }



# -------------------------------------------------------------------------
# Full Revamp compare_matrices function -----------------------------------
# -------------------------------------------------------------------------
compare_matrices_new <-
  function(mat1,
           mat2,
           # number of bootstrap samples
           bootstrap_samples = 0,
           # c("Beta", "Rho")
           parameter_type = "Beta",
           # must be a list containing 2 matrices
           H1_custom_priors = NULL,
           # normal SD
           H1_prior_scale = 0.5,
           # "normal" = SD, "uniform" = range / 2
           H0_prior_scale = NULL,
           # c("normal", "uniform", "posterior_uncertainty")
           H0_distribution = "uniform",
           # plot densities for null, posterior, and prior
           plot = TRUE,
           # upper limit of the x-axis
           plot_xlim = NULL,
           # upper limit of the y-axis
           plot_ylim = NULL,
           
           #- Settings for distances
           # Should the posterior distances be returned? (TRUE/FALSE)
           return_dist = FALSE,
           
           #- Settings for Posterior Uncertainty Null
           # Type of null limit: c("sd", "ci")
           null_lim_type = NULL,
           # Scale for SD limit - on the parameters itself
           null_lim_sd_scale = 1,
           # Scale for hdi limit - on the distances
           null_lim_hdi_scale = .99
           ) {
    
    require(Rmpfr)
    
    # browser()
    
    # fisher z-transform partial correlations
    if (parameter_type == "Rho") {
      mat1 <- draws_matrix2list(mat1)
      mat2 <- draws_matrix2list(mat2)
      
      mat1 <-
        purrr::map(mat1, function(x) {
          fisher_z(x) %>% .[lower.tri(.)] %>% as.vector()
        })
      mat2 <-
        purrr::map(mat2, function(x) {
          fisher_z(x) %>% .[lower.tri(.)] %>% as.vector()
        })
      
      mat1 <- do.call(rbind, mat1)
      mat2 <- do.call(rbind, mat2)
    }
    
    #### with bootstrapping ###
    if (bootstrap_samples > 0) {
      # sample rows from posterior iterations
      rows1 <-
        sample(1:nrow(mat1),
               size = bootstrap_samples,
               replace = TRUE)
      rows2 <-
        sample(1:nrow(mat2),
               size = bootstrap_samples,
               replace = TRUE)
      # assemble new posterior draws_matrices
      mat1 <-
        purrr::map(
          .x = rows1,
          .f = function(.x) {
            mat1[.x,]
          }
        ) %>%
        do.call(rbind, .)
      mat2 <-
        purrr::map(
          .x = rows2,
          .f = function(.x) {
            mat2[.x, ]
          }
        ) %>%
        do.call(rbind, .)
      # compute differences between matrices
      diff_post_mat <- abs(mat1 - mat2)
      diff_post <- diff_post_mat %>%
        apply(., 1, sum)
      
      
      ### compute prior difference matrices
      
      # use custom priors if exist
      if (typeof(H1_custom_priors) == "list" &
          length(H1_custom_priors) == 2) {
        loc <- H1_custom_priors[[1]] %>% as.vector()
        scale <- H1_custom_priors[[2]] %>% as.vector()
        diff_prior <- abs(
          lapply(1:length(loc), function(n) {
            rnorm(n = bootstrap_samples,
                  mean = loc[n],
                  sd = scale[n])
          }) %>%
            do.call(cbind, .) -
            lapply(1:length(loc), function(n) {
              rnorm(n = bootstrap_samples,
                    mean = loc[n],
                    sd = scale[n])
            }) %>%
            do.call(cbind, .)
        ) %>%
          matrix(., ncol = ncol(mat1)) %>%
          apply(., 1, sum)
        
        # TODO need potential warning message if priors are not
        # NULL but do not satisfy correct format
        
      } else{
        # use default priors
        diff_prior_mat <-
          abs(
            rnorm(
              n = bootstrap_samples * ncol(mat1),
              mean = 0,
              sd = H1_prior_scale
            ) -
              rnorm(
                n = bootstrap_samples * ncol(mat1),
                mean = 0,
                sd = H1_prior_scale
              )
          ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
      }
    } else{
      ### NO bootstrappping ###
      
      diff_post_mat <- mat1 - mat2
      attr(diff_post_mat, which = "class") <- "matrix"
      
      diff_post_mat <- abs(diff_post_mat)
      
      diff_post <- diff_post_mat %>%
        apply(., 1, sum)
      
      # prior difference matrices
      if (typeof(H1_custom_priors) == "list" &
          length(H1_custom_priors) == 2) {
        loc <- H1_custom_priors[[1]] %>% as.vector()
        scale <- H1_custom_priors[[2]] %>% as.vector()
        
        diff_prior_mat <- abs(
          lapply(1:length(loc), function(n) {
            rnorm(n = max(4e4, dim(mat1)[1]),
                  mean = loc[n],
                  sd = scale[n])
          }) %>%
            do.call(cbind, .) -
            lapply(1:length(loc), function(n) {
              rnorm(n = max(4e4, dim(mat1)[1]),
                    mean = loc[n],
                    sd = scale[n])
            }) %>%
            do.call(cbind, .)
        ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
        # use default priors
        # TODO needs to be adjusted depending on beta vs. rho
      } else{
        diff_prior_mat <- abs(
          rnorm(
            n = max(4e4, dim(mat1)[1]) * ncol(mat1),
            mean = 0,
            sd = H1_prior_scale
          ) -
            rnorm(
              n = max(4e4, dim(mat1)[1]) * ncol(mat1),
              mean = 0,
              sd = H1_prior_scale
            )
        ) %>%
          matrix(., ncol = ncol(mat1))
        diff_prior <- diff_prior_mat %>%
          apply(., 1, sum)
      }
    }
    
    # distribution for null rope range
    if (H0_distribution == "normal") {
      if (is.null(H0_prior_scale)) {
        H0_prior_scale <- .005
      }
      null_lim <- H0_prior_scale * 4 * ncol(mat1)
    }
    if (H0_distribution == "uniform") {
      if (is.null(H0_prior_scale)) {
        H0_prior_scale <- .01
      }

      # Take maximum distance
      null_lim <- H0_prior_scale *  ncol(mat1)
    }
    if (H0_distribution == "posterior_uncertainty") {
      # Combine uncertainty from both networks?
      # sample rows from posterior iterations
      # TODO potentially use tsnet function for that
      rows1_null <-
        sample(1:nrow(mat1),
               size = nrow(mat1),
               replace = FALSE)
      rows2_null <-
        sample(1:nrow(mat2),
               size = nrow(mat2),
               replace = FALSE)
      
      # assemble new posterior draws_matrices
      mat1_null <-
        purrr::map(
          .x = rows1_null,
          .f = function(.x) {
            mat1[.x, ]
          }
        ) %>%
        do.call(rbind, .)
      mat2_null <-
        purrr::map(
          .x = rows2_null,
          .f = function(.x) {
            mat2[.x, ]
          }
        ) %>%
        do.call(rbind, .)
      # compute differences between matrices
      diff_null1 <- abs(mat1 - mat1_null) %>%
        apply(., 1, sum)
      diff_null2 <- abs(mat2 - mat2_null) %>%
        apply(., 1, sum)
      
      # combine differences
      diff_null <- c(diff_null1, diff_null2)
      
      # HDI for Null distribution, CI_high will be used as upper ROPE limit
      if(null_lim_type == "hdi"){
        null_lim <- bayestestR::hdi(diff_null, ci = null_lim_hdi_scale )$CI_high
      } else if(null_lim_type == "sd"){
        sd1 <- sd(mat1_null)
        sd2 <- sd(mat2_null)
        null_lim <- (sd1 + sd2)/2 * null_lim_sd_scale
      } else{
        stop("null_lim_type must be either 'hdi' or 'sd'")
      
      }
    }

    # TODO removed percentage framing here, then adjust that in simulation
    post_in_rope <-
      diff_post[diff_post < null_lim] %>%
      length() / length(diff_post)
    
    prior_in_rope <-
      diff_post[diff_prior < null_lim] %>%
      length() / length(diff_prior)
    
    ### Compute BF for H0
    prec <- 100 # precision for mpfr (floating point)
    
    mean_post_u <- Rmpfr::mpfr(mean(diff_post),precBits = prec)
    sd_post_u <- Rmpfr::mpfr(sd(diff_post),precBits = prec)
    mean_prior_u <- Rmpfr::mpfr(mean(diff_prior),precBits = prec)
    sd_prior_u <- Rmpfr::mpfr(sd(diff_prior),precBits = prec)
    
    # TODO needs to be described
    q <- Rmpfr::mpfr(null_lim,precBits = prec)
    .mpfr_erange_set(value = (1-2^-52)*.mpfr_erange(c("min.emin","max.emax")))
    
    log_BF_01_mpfr <- 
      Rmpfr::pnorm(
        q,
        mean_post_u,
        sd_post_u,
        log.p = TRUE,
        lower.tail = TRUE
      ) -
      Rmpfr::pnorm(
        q,
        mean_prior_u,
        sd_prior_u,
        log.p = TRUE,
        lower.tail = TRUE
      )
    log_BF_01 <- Rmpfr::asNumeric(log_BF_01_mpfr)
    BF_01 <- Rmpfr::asNumeric(exp(log_BF_01_mpfr))
    
    
    # plot prior vs. posterior
    if (isTRUE(plot)) {
      # combine diff_post and diff_prior in a single column of a dataframe
      # with another column as an indicator if it is a posterior or prior
      # distribution
      if (is.null(plot_xlim)) {
        plot_xlim <- max(median(diff_prior),
                         median(diff_post))#,
        #median(diff_null))
        
      }
      if (is.null(plot_ylim)) {
        plot_ylim <- max(max(density.default(diff_post)[["y"]]),
                         max(density.default(diff_prior)[["y"]])) * 1.5
        
      }
      df_samples <- data.frame(posterior = diff_post,
                               prior = diff_prior) %>%
        tidyr::pivot_longer(
          cols = c(posterior, prior),
          names_to = "distribution",
          values_to = "value"
        ) #%>%
      # rbind(data.frame(distribution = "null",
      #                  value = diff_null))
      
      plot <- df_samples %>%
        ggplot2::ggplot() +
        ggplot2::geom_density(aes(value, fill = distribution),
                              alpha = .5) +
        ggplot2::annotate(
          'rect',
          xmin = 0,
          xmax = null_lim,
          ymin = 0,
          ymax = plot_ylim,
          alpha = .4,
          fill = 'grey'
        ) +
        ggplot2::geom_vline(xintercept = null_lim) +
        ggplot2::scale_x_continuous(limits = c(0, plot_xlim),
                                    expand = expansion()) +
        ggplot2::scale_y_continuous(limits = c(0, plot_ylim),
                                    expand = expansion()) +
        ggplot2::labs(x = "Difference", y = "Density") +
        ggplot2::theme_light() +
        ggplot2::theme(panel.grid = element_blank())
      
      print(plot)
    }
    
    # df with log_BF and overlap coef
    df_results <- data.frame(
      round(BF_01, 4),
      round(log_BF_01, 2),
      round(post_in_rope, 2),
      round(prior_in_rope, 2),
      row.names = NULL
    )
    names(df_results) <-
      c("BF_01",
        "log BF_01",
        "posterior < null CI_high",
        "prior < null CI_high")
    
    print(df_results)
    
    # return list with results silently
    if(isFALSE(return_dist)){
      return(invisible(
        list(
          BF_01 = BF_01,
          log_BF_01 = log_BF_01,
          post_below_null_ub = post_in_rope,
          prior_below_null_ub = prior_in_rope,
          null_lim = null_lim
        )
      ))
      
    } else
      return(invisible(
        list(
          BF_01 = BF_01,
          log_BF_01 = log_BF_01,
          post_below_null_ub = post_in_rope,
          prior_below_null_ub = prior_in_rope,
          null_lim = null_lim,
          diff_post = diff_post,
          diff_prior = diff_prior
        )
      ))

    
    # clean up
    rm(mat1,
       mat2,
       diff_post,
       diff_post_mat,
       diff_prior,
       diff_prior_mat,
       diff_null)
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






# -------------------------------------------------------------------------
# Mplus -------------------------------------------------------------------
# -------------------------------------------------------------------------

# Create VAR Syntax in Mplus ----------------------------------------------
mplus_var_syntax <- function(data,
                             varnames = colnames(data),
                             laggedvars = NULL,
                             cores = 1,
                             iterations = 4000,
                             model = NULL,
                             datafile = "mplus_data.dat",
                             syntaxfile = "mplus_model.inp",
                             samplesfile = "mplus_samples.dat") {
  # browser()
  
  # Default variable names (V1, V2, ...)
  if (is.null(varnames)) {
    varnames <- paste0("V", 1:ncol(data))
  }
  
  # Default lagged variables (V1(1), V2(1), ...)
  if (is.null(laggedvars)) {
    laggedvars <- paste0(varnames, "(1)")
  }
  
  # Default model if not specified
  if (is.null(model)) {
    model <-
      paste(paste(
        varnames,
        " ON ",
        paste(varnames, "&1", sep = "", collapse = " "),
        ";",
        sep = ""
      ), collapse = "\n")
  }
  
  # # Write data to Mplus input format using MplusAutomation
  suppressMessages(MplusAutomation::prepareMplusData(data,
                                                     filename = datafile,
                                                     quiet = TRUE))
  
  # Fill in the template
  mplus_template <- "
DATA:	    FILE = @@DATAFILE@@;    ! data file (should be in same folder)

VARIABLE:	NAMES = @@VARNAMES@@;          ! providing names to the variables
          USEVARIABLES = @@USEVARIABLES@@;   ! select variables for the analysis
	        LAGGED = @@LAGGEDVARS@@;  ! creating first-order
                                    ! lagged observed variables
            MISSING = *;            ! missing value code

ANALYSIS:	ESTIMATOR = BAYES;      ! set estimator (must be Bayes for DSEM)
	        PROCESSORS = @@CORES@@;         ! using 2 processors
	        BITERATIONS = @@ITERATIONS@@;   ! choose number of iterations;
                                    ! minimum is now 2000; will be more if
                                    ! the convergence criterion indicates
                                    ! convergence was not reached

MODEL:	    @@MODEL@@

OUTPUT:	    TECH1 TECH8;            ! asking additional output
SAVEDATA: BPARAMETERS IS @@SAMPLESFILE@@; ! saving posterior samples
"
  mplus_syntax <- mplus_template
  
  # Substitute placeholders with actual values
  mplus_syntax <- gsub("@@DATAFILE@@", datafile, mplus_syntax)
  mplus_syntax <-
    gsub("@@VARNAMES@@", paste(varnames, collapse = " "), mplus_syntax)
  mplus_syntax <-
    gsub("@@USEVARIABLES@@",
         paste(varnames, collapse = " "),
         mplus_syntax)
  mplus_syntax <-
    gsub("@@LAGGEDVARS@@",
         paste(laggedvars, collapse = " "),
         mplus_syntax)
  mplus_syntax <- gsub("@@CORES@@", cores, mplus_syntax)
  mplus_syntax <- gsub("@@ITERATIONS@@", iterations, mplus_syntax)
  mplus_syntax <- gsub("@@MODEL@@", model, mplus_syntax)
  mplus_syntax <- gsub("@@SAMPLESFILE@@", samplesfile, mplus_syntax)
  
  # Write syntax to file
  invisible(writeLines(mplus_syntax, syntaxfile))
  
  return(list(syntax = mplus_syntax))
}


# Convert Mplus samples to array ------------------------------------------
convert_mplus_samples <- function(mplus_output_file) {
  # browser()
  # Extract posterior samples
  mplus_res <- MplusAutomation::readModels(mplus_output_file)
  mplus_samples <-
    as.data.frame(do.call(rbind, mplus_res$bparameters$valid_draw))
  
  # Extract columns
  samples_beta <- mplus_samples %>%
    dplyr::select(contains(".ON"))
  samples_sigma <- mplus_samples %>%
    dplyr::select(contains("with"))
  selected_columns <- names(mplus_samples) %>%
    str_detect("\\d+_V\\d$")
  samples_diag <- mplus_samples %>%
    dplyr::select(which(selected_columns))
  nvar <- sum(selected_columns)  # Number of variables
  
  # Convert posterior samples to matrix format
  split_beta <- split(samples_beta, 1:nrow(samples_beta))
  beta_samples <- lapply(split_beta, function(x) {
    matrix(x,
           nrow = nvar,
           ncol = nvar,
           byrow = TRUE)
  })
  # convert this list to an array
  beta_samples <- array(unlist(beta_samples),
                        dim = c(nvar, nvar, length(beta_samples)))
  split_sigma <- split(samples_sigma, 1:nrow(samples_sigma))
  
  
  # Function to extract covariance values and create covariance matrix
  create_cov_mat <- function(dat, diag_values) {
    cov_matrix <- matrix(NA, nrow = nvar, ncol = nvar)
    diag(cov_matrix) <- 0
    for (i in 1:ncol(dat)) {
      col_name <- names(dat)[i]
      values <- as.vector(strsplit(sub(".*_", "", col_name), "\\."))
      
      # Extract row and column indices from the column name
      row_index <- as.numeric(gsub("V", "", values[[1]][1]))
      col_index <- as.numeric(gsub("V", "", values[[1]][3]))
      
      # Assign the value to the corresponding position in the covariance matrix
      cov_matrix[row_index, col_index] <- as.numeric(dat[, i])
      cov_matrix[col_index, row_index] <-
        as.numeric(dat[, i])  # Covariance matrix is symmetric
    }
    diag(cov_matrix) <- as.numeric(diag_values)
    return(cov_matrix)
  }
  
  # this is a very ugly mix of a loop and lapply, sorry...
  sigma_samples <-
    simplify2array(lapply(seq_along(split_sigma), function(i) {
      create_cov_mat(split_sigma[[i]], samples_diag[i, ])
    }))
  
  
  # Convert to partial correlations
  pcor_samples <- array(NA, dim = dim(sigma_samples))
  for (i in 1:dim(sigma_samples)[3]) {
    # take inverse of covariance matrix (i.e. precision)
    # and convert to partial correlation matrix
    tmp <- -stats::cov2cor(solve(sigma_samples[, , i]))
    diag(tmp) <- 0
    pcor_samples[, , i] <- tmp
  }
  
  return(list(
    beta = beta_samples,
    sigma = sigma_samples,
    rho = pcor_samples
  ))
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
