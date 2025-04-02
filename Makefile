all: docker

FILE=centrality-uncertainty

## build docker image (requires root rights for docker)
dbuild: Dockerfile
	docker build \
	-t $(FILE) .

## run docker container
docker: dbuild
	echo "open RStudio Server at http://localhost:8787"
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
	-v $(CURDIR)/scripts/simulation_viz_figures:/home/rstudio/scripts/simulation_viz_figures \
	$(FILE)
