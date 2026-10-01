#!/bin/bash

cat <<EOF > Dockerfile
## Docker container
FROM rocker/rstudio:4.6.1
MAINTAINER Max Salm <maxsalm3@gmail.com>
LABEL Description="Targets environment" Version="1.0"
### The EXPOSE instruction informs Docker that the container listens on the specified network ports at runtime.
EXPOSE 8787
### Extra system level libs if necessary for R packages
RUN apt-get update && apt-get install -y --no-install-recommends pandoc 
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends libcurl4-openssl-dev 
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends libgit2-dev 
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends libmariadbclient-dev     
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends libxml2-dev 
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends librsvg2-dev 
RUN apt-get update && apt-get --yes --force-yes install --no-install-recommends libv8-dev 
RUN apt-get install -y --no-install-recommends libpng-dev && apt-get --yes --force-yes install --no-install-recommends libmagick++-dev
### Core R section
RUN echo "pkgs_needed <- c('devtools', 'checkmate')" > ~/.Rprofile
RUN echo "repo <- 'https://cran.rstudio.com/'" > ~/.Rprofile
RUN R -e "install.packages(pkgs = 'devtools')"
RUN R -e "install.packages(pkgs = 'remotes')"
RUN R -e "remotes::update_packages(packages = TRUE, upgrade = 'always')"
## Installing Quarto
RUN apt-get update && apt-get install -y --no-install-recommends pandoc pandoc-citeproc curl gdebi-core && rm -rf /var/lib/apt/lists/*
RUN R -e "remotes::install_version(package = 'shiny', version = '>= 1.7.4', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'jsonlite', version = '>= 1.8.4', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'ggplot2', version = '>= 3.4.2', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'renv', version = '>= 0.17.3', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'knitr', version = '>= 1.42', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'rmarkdown', version = '>= 2.21', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN R -e "remotes::install_version(package = 'quarto', version = '>= 1.2', dependencies = TRUE, upgrade = 'always', repos = repo)"
RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb
RUN R -e "remotes::install_version(package = 'targets', version = '>= 1.4.1', dependencies = TRUE, upgrade = 'always', repos = repo)"
## Install project-specific dependencies
COPY renv.lock .
COPY renv/activate.R renv/
COPY .Rprofile .
RUN R -e "renv::restore()"
COPY . .
# Run bash when the container launches
CMD ["bash"]
EOF
cat <<EOF > .dockerignore
# Don't load bash/R code into docker image at build time
./*
EOF
# Build docker image
if [[ "$(docker images -q docker_img 2> /dev/null)" == "" ]]; then
echo "Building Docker image"
docker build -t docker_img:0.1.0 -f ./code/docker/Dockerfile .
docker run -i docker_img:0.1.0 quarto check
else
echo "Docker image exists"
# docker build -t docker_img:0.1.0 -f ./code/docker/Dockerfile .
docker run -i docker_img:0.1.0 quarto check
fi
