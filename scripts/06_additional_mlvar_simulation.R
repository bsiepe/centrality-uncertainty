#' ---
#' title: "Uncertainty is central for reliable inferences: 
#'         Using dynamic network features as predictors"
#' subtitle: "Investigations into mlVAR performance"
#' author: 
#'  - name: Björn S. Siepe
#'    orcid: 0000-0002-9558-4648
#'    affiliations: University of Marburg
#'  - name: Matthias Kloft
#'    orcid: 0000-0003-1845-6957
#'    affiliations: University of Marburg  
#'  - name: Fridtjof Petersen
#'    orcid: 0000-0002-4913-8532
#'    affiliations: University of Groningen
#'  - name: Yong Zhang
#'    orcid: 0000-0002-6313-2575
#'    affiliations: University of Groningen
#'  - name: Laura F. Bringmann
#'    orcid: 0000-0002-8091-9935
#'    affiliations: University of Groningen
#'  - name: Daniel W. Heck
#'    orcid: 0000-0002-6302-9252
#'    affiliations: University of Marburg
#' date: "`r Sys.Date()`"
#' format:
#'   html:
#'     toc: true
#'     number-sections: true
#'     theme: cosmo
#'     code-fold: true
#'     code-tools: true
#'     code-summary: "Show the code"
#'     fig-width: 7
#'     fig-height: 4.5
#'     fig-align: "center"
#'     embed-resources: true
#' execute:
#'   message: false
#'   warning: false
#'   eval: true # only show code
#' ---
#' 
#' # Background
#' This script contains the `SimDesign` code for an additional simulation investigation into `mlVAR`, where we try various different things to investigate why mlVARs performance for random effects was poor. 
#' 
#' We first load all relevant packages: 
## ----packages----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(SimDesign)
library(mlVAR)
library(graphicalVAR)
library(gimme)
library(here)
library(future)
library(corpcor)
library(gridExtra)
library(nlme)
source(here::here("scripts", "00_functions.R"))


#' 
#' 
#' 
#' 
#' ## Idea 1.2.: Sparser DGP
#' 
#' ### Data-Generating Processes
#' 
#' Create a pretty sparse DGP, considerably sparser than that used in the simulation: 
## ----create-var-sparse-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' var_mat <- create_var_matrix(6, 
#'                           sparse = TRUE,
#'                           sparsity_proportion = 0.5,
#'                           boost_factor = 1.25)
#' 
#' # set some additional elements to zero
#' var_mat[3,1] <- 0
#' var_mat[5,1] <- 0
#' var_mat[1,3] <- 0
#' var_mat[1,5] <- 0
#' 
#' # colSums(var_mat)                            
#' # rowSums(var_mat)
#' 
#' arr_to_latex(var_mat)
#' 
#' #' 
#' #' 
#' #' Then we can create the corresponding innovation covariance matrix:
#' ## ----create-cov-sparse-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' cov_mat <- create_cov_matrix(6, 
#'                              off_diag_val = .15,
#'                               sparse = TRUE, 
#'                               sparsity_proportion = 0.7,
#'                               boost_factor = 1.25)
#' 
#' cov_mat[3,1] <- 0
#' cov_mat[5,1] <- 0
#' cov_mat[1,3] <- 0
#' cov_mat[1,5] <- 0
#' 
#' # colSums(cov_mat)
#' # rowSums(cov_mat)
#' arr_to_latex(cov_mat)

#' 
#' Create the corresponding pcor matrix:
## ----eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pcor_mat <- -stats::cov2cor(solve(cov_mat))

#' 
#' And the corresponding precision matrix: 
## ----eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# kappa_mat <- solve(cov_mat)

#' 
#' 
#' Save the DGP:
## ----eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# graph_sparse_synth_additional_sim <- list(beta = var_mat,
#                     sigma = cov_mat,
#                     kappa = kappa_mat,
#                     pcor = pcor_mat)
# 
# saveRDS(graph_sparse_synth_additional_sim , here::here("data", "graph_sparse_synth_additional_sim.RDS"))

#' 
#' Reload it: 
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

graph_nonsparse <- readRDS(here::here("data", "graph_nonsparse.RDS"))
graph_sparse_synth_additional_sim <- readRDS(here::here("data", "graph_sparse_synth_additional_sim.RDS"))

#' 
#' 
#' ### Setting parameters
#' 
#' We define the conditions and the fixed parameters for the simulation.
#' 
## ----params-sparse, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
dgp <- c("sparse")

# Number of timepoints
n_tp <- c(60, 900)

heterogeneity <- "high"

# Simulation parameters
# Number of individuals
n_id <- c(75, 200)

# Design Conditions
df_design <- SimDesign::createDesign(
  dgp = dgp,
  n_id = n_id,
  heterogeneity = heterogeneity,
  n_tp = n_tp
  )

# Number of variables
n_var <- 6

# for regression
reg_error_sd = 1

# Should GIMME use only nondirected contemporaneous effects?
gimme_var_only = TRUE

# Should GIMME return partial correlations?
gimme_pcor = TRUE

# Should data for non-Bayesian methods be person-mean centered?
mean_center = TRUE

sim_pars <- list(
  n_id = n_id,
  n_var = n_var,
  reg_error_sd = reg_error_sd,
  graph_nonsparse = graph_nonsparse,
  graph_sparse = graph_sparse_synth_additional_sim,
  gimme_var_only = gimme_var_only,
  mean_center = mean_center
)

#' 
#' 
#' We recompute the true SD: 
#' 
## ----eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # number of simulations to obtain sd
# n_sim_sd <- 30000
# sd_results_strength <- sd_results_outstrength <- sd_results_instrength <- vector("list", length = nrow(df_design))
# 
# for(i in 1:nrow(df_design)){
#   condition <- df_design[i,]
# 
#   dgp_graph <- ifelse(condition$dgp == "sparse",
#                       "graph_sparse",
#                       "graph_nonsparse")
#   beta_sd <- ifelse(condition$heterogeneity == "low",
#                     0.05,
#                     0.075)
#   sigma_sd <- ifelse(condition$heterogeneity == "low",
#                     0.05,
#                     0.075)
#   ml_sim <- sim_gvar_loop(
#                      graph = sim_pars[[dgp_graph]],
#                      beta_sd = beta_sd,
#                      kappa_sd = kappa_sd,
#                      sigma_sd = sigma_sd,
#                      n_person = n_sim_sd,
#                      n_time = condition$n_tp,
#                      n_node = sim_pars$n_var,
#                      max_try = 1000,
#                      verbose = FALSE,
#                      listify = TRUE,
#                      sim_pkg = "mlVAR",
#                      sparse_sim = TRUE,
#                      most_cent_diff_temp = TRUE,
#                      most_cent_diff_temp_min = 0.05,
#                      most_cent_diff_cont = TRUE,
#                      most_cent_diff_cont_min = 0.05,
#                      innov_var_fixed_sigma = TRUE)
# 
#   # Obtain true centralities
#   true_cent <- centrality_mlvar_sim(ml_sim,
#                                   sim_fn = "sim_gvar_loop")
# 
#   strength <- sapply(true_cent$strength, `[`, 1)
#   outstrength <- sapply(true_cent$outstrength, `[`, 1)
#   instrength <- sapply(true_cent$instrength, `[`, 1)
#   strength_sd <- sd(strength)
#   outstrength_sd <- sd(outstrength)
#   instrength_sd <- sd(instrength)
#   sd_results_strength[[i]] <- strength_sd
#   sd_results_outstrength[[i]] <- outstrength_sd
#   sd_results_instrength[[i]] <- instrength_sd
# }
# 
# 
# sd_results <- list(sd_results_strength, sd_results_outstrength, sd_results_instrength)
# names(sd_results) <- c("sd_results_strength", "sd_results_outstrength", "sd_results_instrength")
# 
# 
# saveRDS(sd_results, here::here("data", "true_sd_sparse_additional_sim.RDS"))
# 

#' 
#' As we only need to compute the true SD once, we can simply load it from disk to save time: 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
true_sd <- readRDS(here::here("data", "true_sd_sparse_additional_sim.RDS"))
names(true_sd) <- c("sd_results_strength", "sd_results_outstrength", "sd_results_instrength")
df_design$strength_sd <- unlist(true_sd$sd_results_strength)
df_design$outstrength_sd <- unlist(true_sd$sd_results_outstrength)
df_design$instrength_sd <- unlist(true_sd$sd_results_instrength)

#' 
#' 
#' 
#' 
#' ### Simulating Data
#' 
#' 
## ----generate-sparse, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------
sim_generate <- function(condition, fixed_objects = NULL){
  source(here::here("scripts", "00_functions.R"))

  # obtain fixed params
  SimDesign::Attach(fixed_objects)

  dgp_graph <- ifelse(condition$dgp == "sparse",
                      "graph_sparse",
                      "graph_nonsparse")
  beta_sd <- ifelse(condition$heterogeneity == "low",
                    0.05,
                    0.075)
  sigma_sd <- ifelse(condition$heterogeneity == "low",
                    0.005,
                    0.005)

  # scale kappa random effects w.r.t diagonal elements
  mean_diag_kappa <- mean(diag(fixed_objects[[dgp_graph]]$kappa))
  kappa_sd_low <- 0.005 * mean_diag_kappa
  kappa_sd_high <- 0.005 * mean_diag_kappa

  kappa_sd <- ifelse(condition$heterogeneity == "low",
                     kappa_sd_low,
                     kappa_sd_high)

  ml_sim <- sim_gvar_loop(
                     graph = fixed_objects[[dgp_graph]],
                     beta_sd = beta_sd,
                     kappa_sd = kappa_sd,
                     sigma_sd = sigma_sd,
                     n_person = condition$n_id,
                     n_time = condition$n_tp,
                     n_node = n_var,
                     max_try = 10000,
                     listify = TRUE,
                     sim_pkg = "mlVAR",
                     sparse_sim = TRUE,
                     most_cent_diff_temp = TRUE,
                     most_cent_diff_temp_min = 0.05,
                     most_cent_diff_cont = TRUE,
                     most_cent_diff_cont_min = 0.05,
                     innov_var_fixed_sigma = TRUE)
  if(any(is.na(ml_sim$beta))){
    stop("Generation of Betas failed")
  }

  if(any(is.na(ml_sim$pcor))){
    stop("Generation of PCORs failed")
  }


  # Obtain true centralities
  true_cent <- centrality_mlvar_sim(ml_sim,
                                    sim_fn = "sim_gvar_loop")

  # Obtain true centralities
  instrength <- sapply(true_cent$instrength, `[`, 1)
  outstrength <- sapply(true_cent$outstrength, `[`, 1)
  strength <- sapply(true_cent$strength, `[`, 1)


  # Simulate correlated data with true variance of centralities
  # first column is the density to keep track of the density
  # and for compatibility with an old sim function
  sim_correlated <- function(density,
                             r_true,
                             sd_true){
  y_matrix <- matrix(NA, nrow = length(density), ncol = length(r_true) + 1)
  y_matrix[, 1] <- density

  for (i in seq_along(r_true)) {
    sd_error <- sd_true * sqrt(1 - r_true[i]^2)
    y_matrix[, i + 1] <- r_true[i] * density + rnorm(length(density), mean = 0, sd = sd_error)
  }

  return(y_matrix)
  }

  covariate_cont_strength <- sim_correlated(density = strength,
                                            r_true = c(0, 0.2, 0.4),
                                            sd_true = condition$strength_sd)
  covariate_in_strength <- sim_correlated(density = instrength,
                                          r_true = c(0, 0.2, 0.4),
                                          sd_true = condition$instrength_sd)
  covariate_out_strength <- sim_correlated(density = outstrength,
                                           r_true = c(0, 0.2, 0.4),
                                           sd_true = condition$outstrength_sd)


  if(sd(instrength) == 0){
    stop("No variation in instrength")
  }
  if(sd(outstrength) == 0){
    stop("No variation in outstrength")
  }
  if(sd(strength) == 0){
    stop("No variation in strength")
  }

  # Sanity checks
  cor_strength <- cor(covariate_cont_strength[,1], covariate_cont_strength)
  cor_instrength <- cor(covariate_in_strength[,1], covariate_in_strength)
  cor_outstrength <- cor(covariate_out_strength[,1], covariate_out_strength)


  # Return data and true centralities
  ret_data <- list(
    data = ml_sim$data,
    beta = ml_sim$beta_l,
    kappa = ml_sim$kappa_l,
    pcor = ml_sim$pcor_l,
    covariate_cont_strength = covariate_cont_strength,
    covariate_out_strength = covariate_out_strength,
    covariate_in_strength = covariate_in_strength,
    true_cent = true_cent,
    cor_strength = cor_strength,
    cor_instrength = cor_instrength,
    cor_outstrength = cor_outstrength
  )

  return(ret_data)


}

#' 
#' 
#' 
#' 
#' 
#' ### Analysis
#' 
## ----data-analysis-sparse, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
sim_analyse <- function(condition, dat, fixed_objects = NULL){

  #--- Preparation
  SimDesign::Attach(fixed_objects)

  # Check if mean centering option is true, if so, mean center
  if(isTRUE(fixed_objects$mean_center)){
    dat$data_centered <- lapply(dat$data, function(x){
      x |>
        dplyr::as_tibble() |>
        dplyr::mutate(across(everything(), ~ . - mean(.)))

    })
  } else {
    dat$data_centered <- lapply(dat$data, function(x){
      x |>
        dplyr::as_tibble()

    })
  }

  # Concatenate list of data into dataframe with id column
  df_data_centered <- dplyr::bind_rows(purrr::map(dat$data_centered, dplyr::as_tibble)
                              , .id = "ID") |>
    dplyr::mutate(ID = as.factor(ID))

  df_data <- dplyr::bind_rows(purrr::map(dat$data, dplyr::as_tibble)
                              , .id = "ID") |>
    dplyr::mutate(ID = as.factor(ID))

  n_id <- condition$n_id


  # #--- frequentist mlVAR

  # Fit model
  fit_mlvar <- tryCatch({mlVAR::mlVAR(df_data_centered,
                            vars = paste0("V", seq(1:n_var)),
                            idvar = "ID",
                            estimator = "lmer",
                            contemporaneous = "orthogonal",
                            temporal = "orthogonal",
                            nCores = 1,
                            scale = FALSE)}, error = function(e) NA)
  if(any(is.na(fit_mlvar))){
    stop("mlVAR did not converge.")
  }



  # Obtain centralities
  cent_mlvar <- centrality_mlvar(fit_mlvar)
  dens_temp_mlvar <- unlist(cent_mlvar$dens_temp)
  dens_cont_mlvar <- unlist(cent_mlvar$dens_cont)
  outstrength_mlvar <- lapply(cent_mlvar$outstrength, function(x) unname(x))
  strength_mlvar <- lapply(cent_mlvar$strength, function(x) unname(x))
  instrength_mlvar <- lapply(cent_mlvar$instrength, function(x) unname(x))
  outstrength_mlvar_first <- sapply(cent_mlvar$outstrength, function(x) unname(x[1]))
  instrength_mlvar_first <- sapply(cent_mlvar$instrength, function(x) unname(x[1]))
  strength_mlvar_first <- sapply(cent_mlvar$strength, function(x) unname(x[1]))

  # Create list of subject-specific estimates
  beta_mlvar <- lapply(1:n_id, function(i){
    t(mlVAR::getNet(fit_mlvar,
                  subject = i,
                  type = "temporal",
                  nonsig = "show"))
  })
  pcor_mlvar <- lapply(1:n_id, function(i){
    mlVAR::getNet(fit_mlvar,
                  subject = i,
                  type = "contemporaneous",
                  nonsig = "show")
  })

  # # Calculate fixed effects adjacency matrix
  adj_mat_beta <- t(ifelse(
    mlVAR::getNet(fit_mlvar,
       type = "temporal",
       nonsig = "hide") != 0, 1, 0))
  adj_mat_pcor <- ifelse(
    mlVAR::getNet(fit_mlvar,
       type = "contemporaneous",
       nonsig = "hide") != 0, 1, 0)

  # set effects to zero
  beta_mlvar <- lapply(beta_mlvar, function(i){
    i * adj_mat_beta
  })
  pcor_mlvar <- lapply(pcor_mlvar, function(i){
    i * adj_mat_pcor
  })

  # Fit regression
  reg_mlvar_in_strength <- lapply(2:4, function(i) lm(dat$covariate_in_strength[, i] ~ instrength_mlvar_first))
  reg_mlvar_cont_strength <- lapply(2:4, function(i) lm(dat$covariate_cont_strength[, i] ~ strength_mlvar_first))
  reg_mlvar_out_strength <- lapply(2:4, function(i) lm(dat$covariate_out_strength[, i] ~ outstrength_mlvar_first))

  rm(fit_mlvar)

  #--- Return Results
  # Also return true centralities for comparison later
  ret_results <- list(
    mlvar = list(
      fit_mlvar = list(
        beta = beta_mlvar,
        pcor = pcor_mlvar),
      dens_temp = dens_temp_mlvar,
      dens_cont = dens_cont_mlvar,
      outstrength = outstrength_mlvar,
      strength = strength_mlvar,
      instrength = instrength_mlvar,
      reg_instrength = reg_mlvar_in_strength,
      reg_strength = reg_mlvar_cont_strength,
      reg_outstrength = reg_mlvar_out_strength
    ),
    true_cent = dat$true_cent,
    covariate_cont_strength = dat$covariate_cont_strength,
    covariate_in_strength = dat$covariate_in_strength,
    covariate_out_strength = dat$covariate_out_strength,
    data = dat$data,
    beta = dat$beta,
    kappa = dat$kappa,
    pcor = dat$pcor,
    true_cor_strength = dat$cor_strength,
    true_cor_instrength = dat$cor_instrength,
    true_cor_outstrength = dat$cor_outstrength
  )


  return(ret_results)

}


#' 
#' 
#' 
#' ### Summary
#' 
#' 
## ----summarize-sparse, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
sim_summarise <- function(condition, results, fixed_objects = NULL){

  #--- Preparation
  SimDesign::Attach(fixed_objects)
  ret <- list()
  n_id <- condition$n_id

  # Global settings
  use_lm_beta <- FALSE


  #--- Parameter recovery
  # IMPORTANT: Keep in mind structure of sim object (rows vs. cols)
  #-- RMSE
  summary_calc <- function(results, method, fit, measure, func) {
  sapply(seq_along(results), function(i){
    if(is.null(results[[i]][[method]][[fit]])){
      return(NA)
    }
    else{
      unlist(func(results[[i]][[method]][[fit]][[measure]], results[[i]][[measure]]))
    }
  })
}


  methods <- c("mlvar")
  measures <- c("beta", "pcor")
  funcs <- list("rmse" = rmse_mean_list,
                "mse" = mse_mean_list,
                "bias" = bias_mean_list)

  # Loop to get summaries and save them into the ret list
  rmse_list <- list()
  mse_list <- list()
  for (method in methods) {
      for (measure in measures) {
        for (func_name in names(funcs)) {
          result_name <- paste(func_name, measure, method, sep = "_")
          calc_list <- summary_calc(
            results,
            method,
            paste0("fit_", method),
            measure,
            funcs[[func_name]])
          mean_tmp <- colMeans(apply(as.matrix(calc_list), c(1,2), as.numeric))
          ret[[paste0(result_name, "_mean")]] <- mean(mean_tmp)
          if(func_name == "rmse"){
            # calculate MSE, needed for MCSE of RMSE
            mse_tmp <- summary_calc(
            results,
            method,
            paste0("fit_", method),
            measure,
            funcs[["mse"]])
            mse_mean_tmp <- mean(colMeans(apply(as.matrix(mse_tmp), c(1,2), as.numeric)))
            ret[[paste0(result_name, "_mcse")]] <- sqrt(stats::var(mean_tmp) /
                                                          (4 * n_rep * mse_mean_tmp))
          }
          else if(func_name == "mse"){
            ret[[paste0(result_name, "_mcse")]] <- sqrt(stats::var(mean_tmp) / n_rep)
          }
          else if(func_name == "bias"){
            ret[[paste0(result_name, "_mcse")]] <- sqrt(stats::var(mean_tmp) / n_rep)
          }
        }
      }
  }



  #--- Centrality
  #-- Rank-Order between individuals
  calc_correlation <- function(results, method, measure) {
  sapply(seq_along(results), function(i){
    if(is.null(results[[i]][[method]][[paste0("fit_", method)]])){
      return(NA)
    }
    else{
      stats::cor(results[[i]][[method]][[paste0("dens_", measure)]],
                 unlist(results[[i]]$true_cent[[paste0("dens_", measure)]]),
                 method = "spearman")
    }
  })
  }

  # Bootstrap the MCSE for the rank correlation
  bootstrap_rankcor <- function(data, n_boot){
    bootstrap_res <- vector("numeric", n_boot)
    for(i in 1:n_boot){
      ind <- sample(1:n_boot, replace = TRUE)
      bootstrap_res[i] <- mean(data[ind], na.rm = TRUE)
    }
      sd(bootstrap_res, na.rm = TRUE)
  }

  measures <- c("temp", "cont")
  for (method in methods) {
    for (measure in measures) {
      result_name <- paste("rankcor", measure, method, sep = "_")
      calc_list <- calc_correlation(
        results = results,
        method = method,
        measure = measure)
      mean_tmp <- mean(calc_list)
      sd_tmp <- bootstrap_rankcor(calc_list, 1000)
      ret[[paste0(result_name, "_mean")]] <- mean_tmp
      ret[[paste0(result_name, "_mcse")]] <- sd_tmp
    }
  }


  #-- Most central within individuals
  calc_most_cent_ident <- function(results, method, measure) {
    sapply(seq_along(results), function(i){
      if(is.null(results[[i]][[method]][paste0("fit_", method)])){
        return(NA)
      }
      else{
        mci <- most_cent_ident(results[[i]][[method]][[measure]],
                              results[[i]]$true_cent[[measure]])
        # check if calculation worked to prevent error
        if(!is.vector(mci) |is.null(mci) | is.null(unlist(mci))){
          return(NA)
        }
        colSums(matrix(unlist(mci),
                      nrow = n_id,
                      byrow = FALSE))/n_id
      }
    })
  }

  measures <- c("instrength", "outstrength", "strength")

  for(method in methods){
    for(measure in measures){
      calc_list <- calc_most_cent_ident(
        results = results,
        method = method,
        measure = measure)
        mean_tmp <- mean(calc_list, na.rm = TRUE)
        mcse_tmp <- sqrt(mean_tmp * (1 - mean_tmp) / n_rep)
      if(measure == "outstrength"){
        ret[[paste0("mostcent_beta_", method, "_mean")]] <- mean_tmp
        ret[[paste0("mostcent_beta_", method, "_mcse")]] <- mcse_tmp

      } else if(measure == "strength"){
        ret[[paste0("mostcent_pcor_", method, "_mean")]] <- mean_tmp
        ret[[paste0("mostcent_pcor_", method, "_mcse")]] <- mcse_tmp
      }
    }
  }





  #--- Regression
  #-- Power
  regression_power <- function(results,
                               method,
                               measure,
                               true_coef = c(0,.2,.4),
                               lm_beta = use_lm_beta) {
  sapply(seq_along(results), function(i){
    if(is.null(results[[i]][[method]][[paste0("reg_", measure)]])){
      return(NA)
    }
    else{
      reg_pvals <- sapply(results[[i]][[method]][[paste0("reg_", measure)]],
                          function(x){
                            if(isTRUE(lm_beta)){
                              res <- summary(lm.beta::lm.beta(x))$coefficients
                              res[2,5]
                            } else{
                              res <- summary(x)$coefficients
                              res[2,4]
                            }})

      # Calculate Empirical detection rate/power
      ifelse(reg_pvals < .05, 1, 0)
    }
  })
  }


  #-- Summary
  true_coef <- c(0, .2, .4)
  regression_summary <- function(results,
                              method,
                              measure,
                              true_coef = c(0,.2,.4),
                              lm_beta = use_lm_beta,
                              summary = "rmse") {
  sapply(seq_along(results), function(i){
    if(is.null(results[[i]][[method]][[paste0("reg_", measure)]])){
      return(NA)
    }
    else{
      reg_coefs <- sapply(results[[i]][[method]][[paste0("reg_", measure)]],
                          function(x){
                            if(isTRUE(lm_beta)){
                              lm.beta::lm.beta(x)$standardized.coefficients
                            } else{
                              x$coefficients
                            }})[2,]   # 2nd row -> beta coef
      # Calculate RMSE
      if(summary == "rmse"){
        sqrt( (reg_coefs - true_coef )^2 )
      } else if(summary == "mse"){
        (reg_coefs - true_coef)^2
      } else if(summary == "bias"){
        reg_coefs - true_coef
      }


      }
  })
  }




  # exclude bmlvar because of its different data structure
  methods <- c("mlvar")
  measures <- c("instrength", "outstrength", "strength")
  summaries <- c("rmse", "bias", "mse")


  for (method in methods) {
      for (measure in measures) {
        for (summary in summaries) {
          # Calculate summary
          result_name_summary <- paste(summary, "reg", measure, method, sep = "_")
          calc_list <- regression_summary(
            results = results,
            method = method,
            measure = measure,
            summary = summary)
          mean_tmp <- rowMeans(calc_list, na.rm = TRUE)
          # bit of a hack to add the names for the different correlation
          # levels
          names(mean_tmp) <- paste0(summary,
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mean")
          ret[[paste0(result_name_summary, "_mean")]] <- mean_tmp
          if(summary == "rmse"){
            # calculate MSE
            mse_tmp <- regression_summary(
            results,
            method,
            measure,
            summary = "mse")
            mse_mean_tmp <- rowMeans(mse_tmp, na.rm = TRUE)
            mcse_tmp <- sqrt (apply(calc_list, 1, var, na.rm = TRUE) / (4*n_rep*mse_mean_tmp))
            names(mcse_tmp) <- paste0(summary,
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mcse")
            ret[[paste0(result_name_summary, "_mcse")]] <- mcse_tmp
          } else if(summary == "mse"){
            mcse_tmp <- sqrt(
              apply(calc_list, 1, var, na.rm = TRUE) / n_rep)
            names(mcse_tmp) <- paste0(summary,
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mcse")
            ret[[paste0(result_name_summary, "_mcse")]] <- mcse_tmp
          } else if(summary == "bias"){
            mcse_bias <- sqrt(
              apply(calc_list, 1, var, na.rm = TRUE) / n_rep
            )
            names(mcse_bias) <- paste0(summary,
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mcse")
            ret[[paste0(result_name_summary, "_mcse")]] <- mcse_bias
          }
          # Calculate Power
          result_name_power <- paste("power", "reg", measure, method, sep = "_")
          calc_list_pwr <- regression_power(
            results = results,
            method = method,
            measure = measure)
          pwr_mean <- rowMeans(calc_list_pwr, na.rm = TRUE)
          names(pwr_mean) <- paste0("power",
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mean")
          ret[[paste0(result_name_power, "_mean")]] <- pwr_mean
          pwr_mcse_tmp <- sqrt(
            pwr_mean * (1 - pwr_mean) / n_rep)
          names(pwr_mcse_tmp) <- paste0("power",
                                   "_reg",
                                   c(0, 2, 4),
                                   "_",
                                   measure,
                                   "_",
                                   method,
                                   "_mcse")
          ret[[paste0(result_name_power, "_mcse")]] <- pwr_mcse_tmp
        }
      }
  }

  ret_vec <- unlist(ret, use.names = TRUE)
  return(ret_vec)

}



#' 
#' 
#' ### Executing Simulation
#' 
## ----run-sim-sparse, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
n_rep <- 20
future::plan(multisession, workers = n_rep)

sim_results <- SimDesign::runSimulation(
                                    design = df_design,
                                    replications = n_rep,
                                    generate = sim_generate,
                                    analyse = sim_analyse,
                                    summarise = sim_summarise,
                                    fixed_objects = sim_pars,
                                    parallel = "future",
                                    # parallel = TRUE,
                                    max_errors = 2,
                                    packages = c("tidyverse",
                                                 "gimme",
                                                 "mlVAR",
                                                 "graphicalVAR",
                                                 "lm.beta",
                                                 "bayestestR",
                                                 "posterior",
                                                 "rstan",
                                                 "corpcor",
                                                 "Rcpp"),
                                    save_results = TRUE,
                                    # debug = "generate",
                                    filename = "additional_sparser_dgp"
                                    # save_seeds = TRUE
                                    )

plan(sequential)


