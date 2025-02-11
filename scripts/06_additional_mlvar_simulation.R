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
## ----packages----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(SimDesign)
library(mlVAR)
library(graphicalVAR)
library(gimme)
library(here)
library(future)
library(corpcor)
library(gridExtra)
source(here::here("scripts", "00_functions.R"))

#' 
#' 
#' 
#' 
#' # Idea 1: Different Simulation
#' 
#' We have tried various different simulation settings (different random effects variance, different number of time points), but none of them really improved the performance of mlVAR for outcome prediction with centrality indices. Here, we showcase results with a low innovation variance. Note that the association with the outcome must be off here (because we did not resimulate the true standard deviation for the centralities), but the rank order of centrality should still be recovered if mlVAR worked well. 
#' 
#' ## Data-Generating Processes
#' 
#' Load DGP based on estimated network structures:  
## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# non-sparse Graph to simulate from
graph_nonsparse <- readRDS(here::here("data/graph_semisparse_synth.RDS"))

# sparse DGP
graph_sparse <- readRDS(here::here("data/graph_semisparse_synth.RDS"))

#'
#'
#'
#' ## Setting parameters
#'
#' We define the conditions and the fixed parameters for the simulation.
#'
## ----params, eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
  graph_sparse = graph_sparse,
  gimme_var_only = gimme_var_only,
  mean_center = mean_center
)

#'
#'
#' As we only need to compute the true SD once, we can simply load it from disk to save time:
## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
true_sd <- readRDS(here::here("data", "true_sd_semisparse.RDS"))
names(true_sd) <- c("sd_results_strength", "sd_results_outstrength", "sd_results_instrength")
df_design$strength_sd <- unlist(true_sd$sd_results_strength)
df_design$outstrength_sd <- unlist(true_sd$sd_results_outstrength)
df_design$instrength_sd <- unlist(true_sd$sd_results_instrength)

#'
#'
#'
#' ## Simulating Data
#'
#'
## ----generate, eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#' ## Analysis
#'
## ----data-analysis, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#' ## Summary
#'
#'
## ----summarize, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
#' ## Executing Simulation
#'
## ----run-sim, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_rep <- 5
future::plan(multisession, workers = n_rep)

sim_results <- SimDesign::runSimulation(
                                    design = df_design[c(3,4),],
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
                                    filename = "low_kappa_sd_mlvar"
                                    # save_seeds = TRUE
                                    )

plan(sequential)




#'
#' 
#' Investigate the results: 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' 
#' 
#' 
#' 
#' # Idea 2: Investigate Centrality Estimates in More Detail
#' 
#' Load the results of the full simulation, specifically the condition with $n = 200$ and $tp = 120$:
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' sim_res_4 <- readRDS("~/centrality-uncertainty/sim_full.rds-results_pc04798/results-row-4.rds")
#' 
#' #' 
#' #' 
#' #' ## Plot Centrality Recovery
#' #' 
#' #' Plot correlation of centrality estimates with true centrality: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' centrality_list_mlvar <- list()
#' 
#' for (i in seq_along(sim_res_4$results)) {
#'   df_temp <- data.frame(
#'     mlvar_outstrength = sapply(sim_res_4$results[[i]]$mlvar$outstrength, function(x) x[[1]]),
#'     mlvar_instrength = sapply(sim_res_4$results[[i]]$mlvar$instrength, function(x) x[[1]]),
#'     mlvar_strength = sapply(sim_res_4$results[[i]]$mlvar$strength, function(x) x[[1]]),
#'     true_outstrength = sapply(sim_res_4$results[[i]]$true_cent$outstrength, function(x) x[[1]]),
#'     true_instrength = sapply(sim_res_4$results[[i]]$true_cent$instrength, function(x) x[[1]]),
#'     true_strength = sapply(sim_res_4$results[[i]]$true_cent$strength, function(x) x[[1]])
#'   )
#'   
#'   centrality_list_mlvar[[i]] <- df_temp
#' }
#' 
#' df_centrality_mlvar <- do.call(rbind, centrality_list_mlvar)
#' 
#' # add an identifier column to track which element each row came from
#' df_centrality_mlvar$sim_rep <- rep(seq_along(sim_res_4$results), sapply(centrality_list_mlvar, nrow))
#' 
#' 
#' 
#' df_centrality_mlvar |> 
#'   ggplot(aes(x = true_outstrength, y = mlvar_outstrength)) + 
#'   geom_point() + 
#'   geom_smooth() + 
#'   theme_centrality()+
#'   labs(x = "True Outstrength",
#'        y = "mlVAR Outstrength")
#' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' ## Potential extraction errors
#' #' One possible reason for a surprisingly bad performance of a method in a simulation study could be errors in the simulation code or, more specifically, some extraction function. In this section, we work with the results of our full simulation study. We checked this explanation multiple times and could not find any errors that would explain the results. 
#' #' 
#' #' ### Transposing
#' #' One potential error could be wrong transposing, which would mean that instrength and outstrength are confused with one another.
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' df_centrality_mlvar |> 
#'   ggplot(aes(x = true_outstrength, y = mlvar_instrength)) + 
#'   geom_point() + 
#'   geom_smooth() + 
#'   theme_centrality()+
#'   labs(x = "True Outstrength",
#'        y = "mlVAR Instrength")
#' 
#' #' 
#' #' Obviously, the results don't change much, so this is an unlikely issue. 
#' #' 
#' #' 
#' #' ### Overregularization
#' #' 
#' #' Compare the true and estimated variability of centrality estimates: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' df_centrality_mlvar |> 
#'   group_by(sim_rep) |> 
#'   summarize(across(everything(),list(mean = mean, sd = sd))) |> 
#'   ungroup() |> 
#'   pivot_longer(cols = !sim_rep) |> 
#'   separate_wider_delim(cols = name, delim = "_", names = c("est", "centrality", "summary")) |> 
#'   pivot_wider(id_cols = c(sim_rep, est, centrality), names_from = summary, values_from = value) |> 
#'   group_by(est, centrality) |> 
#'   summarize(average_mean = mean(mean),
#'             average_sd = mean(sd)) |> 
#'   knitr::kable()
#' 
#' #' mlVAR slightly underestimates the variability of temporal centrality coefficients, but severely underestimates the variabli
#' #' 
#' #' 
#' #' 
#' #' 
#' #' # Idea 3: Investigate Recovery of Raw Random Effects
#' #' 
#' #' Create plots showing recovery for each of the coefficients in the first column, i.e., for the outgoing edges from variable 1: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' all_results <- lapply(seq_along(sim_res_4[["results"]]), function(i) {
#'   
#'   col_1_true <- lapply(sim_res_4[["results"]][[i]][["beta"]], function(x) x[,1])
#'   
#'   df_true <- do.call(rbind, col_1_true) |> 
#'     as.data.frame() |> 
#'     mutate(param = "true_est", id = row_number(), iteration = i)
#'   
#'   col_1_mlvar <- lapply(sim_res_4[["results"]][[i]][["mlvar"]][["fit_mlvar"]][["beta"]], function(x) x[,1])
#'   
#'   df_mlvar <- do.call(rbind, col_1_mlvar) |> 
#'     as.data.frame() |> 
#'     mutate(param = "mlvar_est", id = row_number(), iteration = i)
#'   
#'   col_1_bmlvar <- lapply(sim_res_4[["results"]][[i]][["bmlvar"]][["fit_bmlvar"]][["beta"]], function(x) x[,1])
#'   
#'   df_bmlvar <- do.call(rbind, col_1_bmlvar) |> 
#'     as.data.frame() |> 
#'     mutate(param = "bmlvar_est", id = row_number(), iteration = i)
#'   
#'   df_all <- rbind(df_true, df_mlvar, df_bmlvar)
#'   
#'   df_all_wide <- df_all |> 
#'     pivot_wider(id_cols = c(id, iteration), names_from = param, values_from = c(V1, V2, V3, V4, V5, V6)) |> 
#'     mutate(
#'       across(starts_with("V"), 
#'              .fns = list(mlvar_bias = ~ get(sub("true_est", "mlvar_est", cur_column())) - ., 
#'                          bmlvar_bias = ~ get(sub("true_est", "bmlvar_est", cur_column())) - .), 
#'              .names = "{.col}_{.fn}")
#'     )
#'   
#'   return(df_all_wide)
#' })
#' 
#' df_combined <- all_results |> 
#'     bind_rows() |> 
#'     select(!contains("mlvar_bias")) |> 
#'     select(!contains("bmlvar_bias")) 
#' 
#' df_combined_2 <- rename_with(df_combined, ~gsub("_est", "", .x))
#' df_ests <- df_combined_2 |>   
#'   pivot_longer(
#'     cols = !c(id, iteration), 
#'     names_to = c("variable", "method"), 
#'     names_sep = "_", 
#'     values_to = "est"
#'   ) |> 
#'   pivot_wider(
#'     id_cols = c(id, iteration, variable),
#'     names_from = method, 
#'     values_from = est
#'   )
#' 
#' df_ests |> 
#'   ggplot(aes(y = mlvar, x = true))+
#'   geom_point()+
#'   geom_smooth()+
#'   facet_grid(~variable) + 
#'   theme_centrality()+
#'   labs(x = "True Estimate",
#'        y = "mlVAR Estimate")
#' 
#' 
#' #' 
#' #' That's how the results look for the Bayesian mlVAR approach: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' df_ests |> 
#'   ggplot(aes(y = bmlvar, x = true))+
#'   geom_point()+
#'   geom_smooth()+
#'   facet_grid(~variable) + 
#'   theme_centrality()+
#'   labs(x = "True Estimate",
#'        y = "BmlVAR Estimate")
#' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' We can additionally create heatmaps indicating the bias for each estimate for each individual. We do not show the output here, as it is a long PDF - but we found no systematic issues in these investigations. 
#' #' 
#' #' 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' generate_heatmaps <- function(true_matrices, estimated_matrices, output_file, n_ind) {
#'   
#'   if (length(true_matrices) != length(estimated_matrices)) {
#'     stop("The lists of matrices must have the same length.")
#'   }
#'   
#'   pdf(file = output_file, width = 8, height = 10)
#'   
#'   # loop over each pair
#'   for (i in 1:n_ind) {
#'     
#'     diff_matrix <- true_matrices[[i]] - estimated_matrices[[i]]
#'     
#'     # convert the difference matrix to a data frame for ggplot
#'     diff_df <- reshape2::melt(diff_matrix)
#'     
#'     p <- ggplot(diff_df, aes(x = Var1, y = Var2, fill = value)) +
#'       geom_tile() +
#'       scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
#'       # plot the numerical value
#'       geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
#'       labs(title = paste("Difference Matrix for Individual", i),
#'            x = "",
#'            y = "",
#'            fill = "Difference") +
#'       theme_minimal()
#'     
#'     print(p)
#'   }
#'   
#'   dev.off()
#' }
#' 
#' 
#' generate_heatmaps(sim_res_4$results[[1]]$beta, sim_res_4$results[[1]]$mlvar$fit_mlvar$beta, "heatmaps.pdf", n_ind = 200)
#' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' # Idea 4: Different DGP, long time series
#' #' Here, we don't simulate any external outcome, but rather simulate from a single DGP with a very long time series to check random effects recovery in this somewhat idealized situation. This implicitly also checks whether our simulation setup might be at fault for the poor performance. 
#' #' 
#' #' 
#' #' ## New DGP
#' #' 
#' #' We use a grpah obtained by estimating mlVAR models on data by Fried et al. (2022). More information on the fitting process can be found in the file `04_dgps.qmd`. 
#' #' 
#' #' Load the graph:
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' graph_sparse <- readRDS(here("data", "graph_sparse.RDS"))
#' 
#' # compute Sigma (covariance matrix)
#' graph_sparse$sigma <- solve(graph_sparse$kappa)
#' 
#' 
#' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' ## Simulate Data
#' #' 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' sim_data <- sim_gvar_loop(
#'                      graph = graph_sparse,
#'                      beta_sd = .075,
#'                      kappa_sd = NA,
#'                      sigma_sd = .075, 
#'                      n_person = 200,
#'                      n_time = 1000,
#'                      n_node = 6,
#'                      max_try = 10000,
#'                      listify = TRUE,
#'                      sim_pkg = "mlVAR",
#'                      sparse_sim = TRUE,
#'                      most_cent_diff_temp = FALSE,
#'                      most_cent_diff_cont = FALSE,
#'                      innov_var_fixed_sigma = FALSE)
#' 
#' #' 
#' #' ## Estimate
#' #' 
#' ## ----message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' df_data <- dplyr::bind_rows(purrr::map(sim_data$data, dplyr::as_tibble)
#'                               , .id = "ID") |> 
#'     dplyr::mutate(ID = as.factor(ID))
#' 
#' fit_mlvar <- mlVAR(df_data,
#'                    vars = paste0("V", seq(1:6)),
#'                    idvar = "ID",
#'                    estimator = "lmer",
#'                    contemporaneous = "correlated",
#'                    temporal = "correlated",
#'                    nCores = 4,
#'                    scale = TRUE)
#' 
#' #' 
#' #' 
#' #' ## Check Recovery
#' #' 
#' #' Compare estimates: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' mlvar_beta <- lapply(fit_mlvar$results$Beta$subject, function(x) x[,,1])
#' true_beta <- sim_data$beta_l
#' 
#' #' 
#' #' Compare them individually:
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' bias_list <- mapply(`-`, mlvar_beta, true_beta, SIMPLIFY = FALSE)
#' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' Plot recovery for each variable: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' agg_data <- do.call(rbind, lapply(1:200, function(i) {
#'   data.frame(True = c(true_beta[[i]]), Estimate = c(mlvar_beta[[i]]),
#'              Variable = rep(1:36, each = 1))
#' }))
#' 
#' plot_list <- lapply(1:36, function(var) {
#'   data <- subset(agg_data, Variable == var)
#'   ggplot(data, aes(x = True, y = Estimate)) +
#'     geom_point(alpha = 0.5) +
#'     geom_smooth(formula = y ~ x, method = loess)+
#'     theme_minimal() +
#'     ggtitle(paste("Variable", var))
#' })
#' 
#' 
#' gridExtra::grid.arrange(grobs = plot_list, ncol = 6)
#' 
#' 
#' plot_list[[3]]
#' 
#' #' 
#' #' Sanity check: What happens if I transpose the matrix:
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' mlvar_beta_transposed <- lapply(mlvar_beta, t)
#' 
#' agg_data <- do.call(rbind, lapply(1:200, function(i) {
#'   data.frame(True = c(true_beta[[i]]), Estimate = c(mlvar_beta_transposed[[i]]),
#'              Variable = rep(1:36, each = 1))
#' }))
#' 
#' plot_list_transposed <- lapply(1:36, function(var) {
#'   data <- subset(agg_data, Variable == var)
#'   ggplot(data, aes(x = True, y = Estimate)) +
#'     geom_point(alpha = 0.5) +
#'     geom_smooth()+
#'     theme_minimal() +
#'     ggtitle(paste("Variable", var))
#' })
#' 
#' gridExtra::grid.arrange(grobs = plot_list_transposed, ncol = 6)
#' 
#' #' Nothing else.
#' #' 
#' #' 
#' #' We can additionally calculate the within-matrix correlation of true and estimates effects for each individual: 
#' ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' beta_cors <- sapply(1:200, function(i) {
#'   cor(c(true_beta[[i]]), c(mlvar_beta[[i]]))
#' })
#' hist(beta_cors)
#' 
#' #' 
#' #' This means that for each individual, the estimated betas show good recovery, but not across individuals. That also explains why recovery of the most central edge works alright, but relating it to an outcome does not work. 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' 
#' #' # Write to standard R script
#' #' 
#' #' To run the simulation on the server, it can be easier to just execute an R script.
#' #' 
#' ## ----eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' # knitr::purl(here::here("scripts", "06_additional_mlvar_simulation.qmd"),
#' #             output = here::here("scripts", "06_additional_mlvar_simulation.R"),
#' #             documentation = 2)
#' 
#' #' 
