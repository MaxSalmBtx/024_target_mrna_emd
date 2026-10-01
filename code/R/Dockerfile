# 1. Start from a versioned base image from the Rocker project
FROM rocker/r-ver:4.6.1

# 2. Install system dependencies
# Add any system libraries your R packages might need (e.g., libcurl, libxml2)
# This layer is cached as long as the system dependencies don't change
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Copy only the renv files to the container
# This is the key step for caching. The subsequent steps will only re-run
# if you change your renv.lock file (i.e., your R package dependencies).
COPY renv.lock .
COPY renv/activate.R renv/
COPY .Rprofile .

# 4. Install renv and restore the R packages from the lockfile
# This is the most time-consuming step, but it will be cached.
RUN R -e "install.packages('renv')"
RUN R -e "renv::restore()"

# 5. Copy the rest of your project code
# This layer is separate from package installation. Changes to your scripts
# will re-run from here, but won't trigger a full package reinstall.
COPY . .

# 6. Define the command to run when the container starts
# For example, to run an R script named 'main.R'
# CMD ["Rscript", "00_run_targets.R"]


# docker build -t 024_target_emd:1.0 .
