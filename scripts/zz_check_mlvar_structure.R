source(here::here("scripts", "00_functions.R"))
library(mlVAR)

# Test the mlVAR extraction functions
graph_sparse <- readRDS(here::here("data/graph_sparse_synth_new.RDS"))

graph_cut <- lapply(graph_sparse, function(x){
  x[1:5, 1:5]
})

ml_sim <- sim_gvar_loop(
  graph = graph_cut,
  beta_sd = .1,
  kappa_sd = .1,
  sigma_sd = .1, 
  n_person = 300,
  n_time = 300,
  n_node = 5,
  max_try = 10000,
  listify = TRUE,
  sim_pkg = "mlVAR",
  sparse_sim = TRUE,
  most_cent_diff_temp = TRUE,
  most_cent_diff_temp_min = 0.1,
  innov_var_fixed_sigma = TRUE)


# Concatenate list of data into dataframe with id column
df_data <- dplyr::bind_rows(purrr::map(ml_sim$data, dplyr::as_tibble)
                            , .id = "ID") |>  
  dplyr::mutate(ID = as.factor(ID))

# estimate the model
fit_mlvar <- mlVAR::mlVAR(
  df_data,
  vars = paste0("V", seq(1:5)),
  idvar = "ID",
  estimator = "lmer",
  contemporaneous = "correlated",
  temporal = "correlated",
  nCores = 1,
  scale = FALSE
)
