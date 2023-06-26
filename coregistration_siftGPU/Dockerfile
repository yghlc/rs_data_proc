FROM ubuntu:16.04

# Install base utilities
RUN apt-get update && apt-get install -y  \
        build-essential \
        pkg-config \
        zip \
        unzip \
        git \
	    wget \
	    && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*

# Install libary needed for building SiftGPU
RUN apt-get update && apt-get install -y  \
	libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev libdevil-dev libglew-dev && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y  \
	cmake && \
	apt-get clean && rm -rf /var/lib/apt/lists/*


# default install GDAL 1.11.3 for Ubuntu 16.04 (good, we need gdal version 1)
RUN  apt-get update && apt-get install -y gdal-bin libgdal-dev && \
	apt-get clean && rm -rf /var/lib/apt/lists/*

## Install miniconda
#ENV CONDA_DIR /opt/miniconda3
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/sh ~/miniconda.sh -b -p /opt/miniconda3  && rm ~/miniconda.sh
#
## Make RUN commands use the new environment (conda base): 
#SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]
#
## Put conda in path so we can use conda activate
#ENV PATH=$CONDA_DIR/bin:$PATH
#
## t: tarball, i index-cache p: package  y: Do not ask for confirmation
#RUN conda clean -tipy
#
#COPY run.sh .
## WORKDIR /data
#
#CMD ["conda", "run", "--no-capture-output", "-n", "base","/bin/bash","-c","./run.sh"]

