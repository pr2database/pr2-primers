# Dockerise a shiny app on Google cloud

## Major steps

1. Make sure that the shiny application works on R studio
2. Build docker file 
  * from shiny-verse
  * Install packages as needed (load install_packages.R)
  * Use the shiny-server command
3. Test docker locally
4. Upload to GitHub, this will automatically force compilation on Google cloud and Docker site.

If the Docker application has not been yet uploaded to Google
1. Upload to Google cloud
  * Set up enough memory (2 Go for pr2 primers)
2. Upload to Docker web site

## Building docker image

* Do not use package renv which is really heavy....
* Start Powershell under windows

```
cd C:/daniel.vaulot@gmail.com/Databases/_PR2/pr2-primers


# Buid with cache

docker build . -t pr2-primers

docker run --rm -p 8080:8080 pr2-primers

```

Test locally

* http://localhost:8080/


## Push to Cloud run


Utilize Google Builds to build image on the cloud

```
gcloud auth login

gcloud auth configure-docker

gcloud config set project tactile-bolt-247111

gcloud builds submit --tag asia.gcr.io/tactile-bolt-247111/pr2-primers
```

Deploy to Google Cloud Run (Need only first time)

```
gcloud run deploy --image asia.gcr.io/tactile-bolt-247111/pr2-primers --platform managed --max-instances 1
```

Effectuer ensuite un mappage de domaine sur:

http://app.pr2-primers.org


DO NOT USE - Push image to Google Registry 
gcloud auth login

gcloud auth configure-docker

docker tag pr2-primers asia.gcr.io/tactile-bolt-247111/pr2-primers:v1.0.2

docker push asia.gcr.io/tactile-bolt-247111/pr2-primers
```


## Push container to Docker repository

* Can also be done with Docker Desktop

```
docker images

docker tag pr2-primers vaulot/pr2-primers:v1.0.2

docker push vaulot/pr2-primers:v1.0.2
```

## Docker misc

* List running containers

```
docker container ls
```

* Stop a container
```
docker stop 12a32e8928ef
```

* Remove dangling caches
```
docker builder prune
docker image prune
```

* Buid without cache
```
docker build --no-cache . -t asia.gcr.io/aerial-citron-246112/pr2-primers
```

* Management

```
# see All containers
docker ps -a

# See images and remove

docker images
docker rmi xxxxx

# Remove containers, image, cache and volumes

docker system prune --volumes
```

## Example Dockerfile

```

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


# install directly the packages

RUN Rscript install_packages.R

# For testing
# CMD Rscript R/test.R


# expose port

EXPOSE 8080

USER shiny

# avoid s6 initialization
# see https://github.com/rocker-org/shiny/issues/79

# The next line prevents the application to start on Google
# CMD ["R", "-e", "shiny::runApp(port = 8080)"]

# Better to use
CMD ["/usr/bin/shiny-server"]

```

## Example .dockerignore


```
# .Rprofile is needed to start R studio but not for the shiny app
.Rprofile
rsconnect
archives
sandbox
renv
```