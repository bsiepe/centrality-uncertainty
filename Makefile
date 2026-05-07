IMAGE_VIZ=centrality-uncertainty-viz
IMAGE_SIM=centrality-uncertainty-sim

## Default: show available targets
.PHONY: all rstudio dbuild dbuild-sim render render-viz render-dgps render-sim render-additional smoke-test help

all: help

## Print available targets
help:
	@echo "Available targets:"
	@echo "  make rstudio           - Launch interactive RStudio Server at localhost:8787"
	@echo "  make render-viz        - Render visualization scripts (fast, uses pre-computed output)"
	@echo "  make render-dgps       - Render 04_dgps.qmd (generates data/ files)"
	@echo "  make render-sim        - Render 01_centrality_simulation.qmd (WARNING: long-running)"
	@echo "  make render-additional - Render additional simulations (WARNING: long-running)"
	@echo "  make smoke-test        - Check that all R packages load correctly (fast)"
	@echo "  make dbuild            - Build lightweight visualization image (Dockerfile)"
	@echo "  make dbuild-sim        - Build full simulation image (Dockerfile.sim, includes Stan)"

## Build lightweight visualization image (skips if already up to date)
dbuild: Dockerfile
	docker image inspect $(IMAGE_VIZ) > /dev/null 2>&1 || docker build -t $(IMAGE_VIZ) -f Dockerfile .

## Build full simulation image with Stan (skips if already up to date)
dbuild-sim: Dockerfile.sim
	docker image inspect $(IMAGE_SIM) > /dev/null 2>&1 || docker build -t $(IMAGE_SIM) -f Dockerfile.sim .

## Launch interactive RStudio Server using the visualization image
## Open at http://localhost:8787 — run scripts manually from there
rstudio: dbuild
	@echo "Open RStudio Server at http://localhost:8787"
	docker run \
	--rm \
	-ti \
	-e DISABLE_AUTH=true \
	-e ROOT=true \
	-e USERID=$(id -u) \
	-e GROUPID=$(id -g) \
	-p 8787:8787 \
	-v $(CURDIR)/scripts:/home/rstudio/scripts \
	-v $(CURDIR)/data:/home/rstudio/data \
	-v $(CURDIR)/output:/home/rstudio/output \
	-v $(CURDIR)/figures:/home/rstudio/figures \
	$(IMAGE_VIZ)

## Verify visualization packages load correctly (fast smoke test)
smoke-test: dbuild
	docker run \
	--rm \
	$(IMAGE_VIZ) \
	Rscript -e " \
	  pkgs <- c('tidyverse','here','cowplot','ggh4x','pander','MetBrewer', \
	            'sysfonts','showtext','ggokabeito','janitor','svglite','lm.beta'); \
	  invisible(lapply(pkgs, function(p) { \
	    if (!requireNamespace(p, quietly = TRUE)) stop('Missing package: ', p); \
	    message('OK: ', p) \
	  })); \
	  message('All packages available.') \
	"

## Render visualization scripts using pre-computed output files (fast)
render-viz: dbuild
	docker run \
	--rm \
	-v $(CURDIR)/scripts:/home/rstudio/scripts \
	-v $(CURDIR)/data:/home/rstudio/data \
	-v $(CURDIR)/output:/home/rstudio/output \
	-v $(CURDIR)/figures:/home/rstudio/figures \
	$(IMAGE_VIZ) \
	bash -c " \
	  cd /home/rstudio && \
	  quarto render scripts/05_simulation_viz.qmd && \
	  quarto render scripts/05_simulation_viz_heterogeneous.qmd && \
	  quarto render scripts/09_additional_simulations_viz.qmd \
	"

## Render only the DGP script (generates data/ files used by simulations)
render-dgps: dbuild-sim
	docker run \
	--rm \
	-v $(CURDIR)/scripts:/home/rstudio/scripts \
	-v $(CURDIR)/data:/home/rstudio/data \
	$(IMAGE_SIM) \
	bash -c "cd /home/rstudio && quarto render scripts/04_dgps.qmd"

## Render main simulation (WARNING: computationally intensive, HPC recommended)
render-sim: dbuild-sim
	docker run \
	--rm \
	-v $(CURDIR)/scripts:/home/rstudio/scripts \
	-v $(CURDIR)/data:/home/rstudio/data \
	-v $(CURDIR)/output:/home/rstudio/output \
	$(IMAGE_SIM) \
	bash -c "cd /home/rstudio && quarto render scripts/01_centrality_simulation.qmd --execute"

## Render additional simulations (WARNING: computationally intensive, HPC recommended)
render-additional: dbuild-sim
	docker run \
	--rm \
	-v $(CURDIR)/scripts:/home/rstudio/scripts \
	-v $(CURDIR)/data:/home/rstudio/data \
	-v $(CURDIR)/output:/home/rstudio/output \
	$(IMAGE_SIM) \
	bash -c " \
	  cd /home/rstudio && \
	  quarto render scripts/09_additional_simulations_no_bayes.qmd --execute && \
	  quarto render scripts/09_additional_simulations_only_bayes.qmd --execute \
	"
