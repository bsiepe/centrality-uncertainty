## Lightweight image for visualization (no Stan/simulation packages).
## Output files are mounted as volumes at runtime — not baked into the image.
## R version: https://hub.docker.com/r/rocker/verse/tags
FROM rocker/verse:4.5.1

## set up directories
RUN mkdir -p /home/rstudio/scripts/models \
             /home/rstudio/scripts/simulation_viz_figures \
             /home/rstudio/data \
             /home/rstudio/output \
             /home/rstudio/figures

COPY centrality-uncertainty.Rproj /home/rstudio/

## copy scripts (data and output are mounted as volumes at runtime)
COPY scripts/00_functions.R /home/rstudio/scripts/
COPY scripts/04_dgps.qmd /home/rstudio/scripts/
COPY scripts/05_simulation_viz.qmd /home/rstudio/scripts/
COPY scripts/05_simulation_viz_heterogeneous.qmd /home/rstudio/scripts/
COPY scripts/09_additional_simulations_viz.qmd /home/rstudio/scripts/
COPY scripts/references.bib /home/rstudio/scripts/

## visualization packages only
RUN install2.r --error --skipinstalled --ncpus -1 \
    tidyverse here cowplot ggh4x pander MetBrewer \
    sysfonts showtext ggokabeito janitor svglite lm.beta quarto \
    SimDesign corpcor 
