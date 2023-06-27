#!/bin/bash

############################################################
## build and manage docker images
#build docker image, repository name must be lowercase
docker build -t  siftgpu-image-match .

#lists the images you have locally.
docker images

#remove local image, You can refer to an image by its ID or its name  (-f: force)
docker rmi  -f siftgpu-image-match 



############################################################
## run and manage docker containers

# create and run container, --rm: remove the container on exit
docker run --rm -it siftgpu-image-match 

# run, mount /data (-v PATH-on-host-machine:PATH-inside-container)
# on UVic desktop
docker run -v /data/LingcaoHuang:/data -v $HOME:/home/hlc  -it siftgpu-image-match 
docker run --rm -v /data/LingcaoHuang:/data -v $HOME:/home/hlc  -it siftgpu-image-match 
# on my Laptop:
docker run --rm -v ${HOME}/Data:/data -v ${HOME}:/home/hlc  -it siftgpu-image-match 

#  output contain ids
docker ps -a | grep -v hello-world | grep -v IMAGE |  awk '{print $1}' > container_ids.txt


# bulid with tag and ID of hub.docker.com, then push to hub.docker.com
docker build -t yghlc/siftgpu-image-match:v1 .
docker push yghlc/siftgpu-image-match:v1

# using GPU
# Since Docker 19.03, you need to install nvidia-container-toolkit package and then use the --gpus all flag.

### launch a new terminal to the container, e9ef58868d14 is the container by "nvidia-docker ps" or "nvidia-docker ps -a"
# docker exec -it 0659b37a0deb bash

### start the container at the background
#4cc63f4a50d1 is got by "nvidia-docker ps -q -l"
# docker start e9ef58868d14

### attach to the container
# docker attach e9ef58868d14


# use CMD to set default program to run when nothing input. 
