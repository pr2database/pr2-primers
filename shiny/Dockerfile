# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny-verse

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

WORKDIR /srv/shiny-server

# copy necessary files

COPY *.md ./
COPY *.R ./
COPY /www  ./www
COPY /data  ./data
COPY /R  ./R

# install renv & restore packages

# COPY renv.lock ./
# COPY /renv  ./renv
# RUN Rscript R/renv.R

# install directly the packages

RUN Rscript install_packages.R

# For testing
# CMD Rscript R/test.R


# expose port

EXPOSE 8080

USER shiny

# avoid s6 initialization
# see https://github.com/rocker-org/shiny/issues/79

# CMD ["R", "-e", "shiny::runApp(port = 8080)"]
CMD ["/usr/bin/shiny-server"]
