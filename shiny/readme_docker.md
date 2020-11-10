# Building docker image

* Do not use package renv which is really heavy....

```
cd C:/daniel.vaulot@gmail.com/Databases/_PR2/pr2-primers/shiny


# Buid with cache

docker build . -t asia.gcr.io/aerial-citron-246112/pr2-primers
docker build . -t pr2-primers

```

Test locally

* http://localhost:8080/

```
docker run --rm -p 8080:8080 asia.gcr.io/aerial-citron-246112/pr2-primers
```

Push image to Google Registry
```
gcloud auth login

gcloud auth configure-docker

docker push asia.gcr.io/aerial-citron-246112/pr2-primers
```

Alternatively, ultilize Google Builds to build image

```
gcloud auth login

gcloud auth configure-docker

gcloud config set project aerial-citron-246112

gcloud builds submit --tag asia.gcr.io/aerial-citron-246112/pr2-primers
```

Deploy to Google Cloud Run
```
gcloud run deploy --image asia.gcr.io/aerial-citron-246112/pr2-primers --platform managed --max-instances 1
```


# Docker misc

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
```

* Buid without cache
```
docker build --no-cache . -t asia.gcr.io/aerial-citron-246112/pr2-primers
```

* Remove images

docker images
docker rmi #container
