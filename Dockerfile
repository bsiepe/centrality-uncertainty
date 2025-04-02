## set R version (https://hub.docker.com/r/rocker/verse/tags)
FROM rocker/verse:4.4.2

## set up directories
RUN mkdir -p /home/rstudio/scripts/simulation_viz_figures /home/rstudio/data /home/rstudio/output /home/rstudio/figures
COPY centrality-uncertainty.Rproj /home/rstudio/

## copy specific files
COPY scripts/00_functions.R /home/rstudio/scripts/
COPY scripts/05_simulation_viz.qmd /home/rstudio/scripts/
COPY output/sim_results.RDS /home/rstudio/output/
COPY output/sim_full.rds /home/rstudio/output/


## install R packages from CRAN the last day of the specified R version
## ncpus set to -1 (all available cores)
RUN install2.r --error --skipinstalled --ncpus -1 \
    tidyverse SimDesign here cowplot ggh4x pander MetBrewer sysfonts showtext ggokabeito janitor svglite rmarkdown knitr

