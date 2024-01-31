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
  
  #--- Prepare 
  n_var <- fit$n_vars_total
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
  
  
  #--- Density
  dens_temp <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      sum(abs(x[, temp_ind]))
    }
  })
  
  dens_cont <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      x <- x[, cont_ind]
      diag(x) <- 0
      x[lower.tri(x)] <- 0L
      sum(abs(x))
    }
  })

  #--- Centrality
  outstrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      colSums(abs(x[, temp_ind]))
    }
  })
  instrength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      rowSums(abs(x[, temp_ind]))
    }
  })
  strength <- lapply(fit$path_est_mats, function(x){
    if(!is.double(x)){
      NA
    }
    else{
      x <- x[, cont_ind]
      diag(x) <- 0
      x[lower.tri(x)] <- 0L
      colSums(abs(x))
    }
  })
  
  # return list
  return(list(outstrength = outstrength,
              instrength = instrength,
              strength = strength,
              dens_temp = dens_temp,
              dens_cont = dens_cont))
  
  
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
