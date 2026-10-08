FROM rocker/verse:4.6.1 AS builder
RUN id -u rstudio >/dev/null 2>&1 || useradd -m -d /home/rstudio -s /bin/bash rstudio
ARG RENV_PATHS_CACHE=/home/rstudio/.cache/R/renv
ENV "RENV_PATHS_CACHE"="${RENV_PATHS_CACHE}"
RUN apt-get update -y && apt-get install -y  cmake make libuv1-dev pandoc libx11-dev libcurl4-openssl-dev libssl-dev zlib1g-dev libicu-dev git libsodium-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = FALSE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("renv")'
COPY renv.lock.prod renv.lock
RUN mkdir -p "${RENV_PATHS_CACHE}" && chown -R rstudio:rstudio "${RENV_PATHS_CACHE}"
USER rstudio
RUN --mount=type=cache,id=renv-cache,target=${RENV_PATHS_CACHE} R -e 'renv::restore()'
FROM builder AS final
COPY MS2Search_*.tar.gz /app.tar.gz
RUN R -e 'remotes::install_local("/app.tar.gz",upgrade="never")'
RUN rm /app.tar.gz
EXPOSE 80
USER rstudio
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0',golem.app.prod=TRUE);library(MS2Search);MS2Search::run_app()"
