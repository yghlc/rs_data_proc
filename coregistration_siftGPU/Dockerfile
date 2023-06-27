FROM ubuntu:14.04

## for old Ubuntu, need to change the link for apt-get: see https://stackoverflow.com/questions/30316812/ubuntu-apt-get-unable-to-fetch-packages
#RUN cp /etc/apt/sources.list /etc/apt/sources.list.bak  && \
#	sed -i -e s%archive.ubuntu.com%old-releases.ubuntu.com%g /etc/apt/sources.list

#	sed -i -e s%security%old-releases%g /etc/apt/sources.list

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
#
# Install library needed for building SiftGPU and others
RUN apt-get update && apt-get install -y  \
	libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev libdevil-dev libglew-dev  \
	cmake vim && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

#
## default install GDAL 1.11.3 for Ubuntu 16.04 (good, we need gdal version 2.0.3)
#RUN  apt-get update && apt-get install -y gdal-bin libgdal-dev && \
#	apt-get clean && rm -rf /var/lib/apt/lists/*

#install GDAL 2.0.3 (needed by ImageMatchsiftGPU), install to /usr/local
COPY gdal-2.0.3.tar.gz .
RUN tar -zxf gdal-2.0.3.tar.gz&& cd gdal-2.0.3 && ./configure && make -j 4 && make install && \
	cd .. && rm -rf gdal-2.0.3*
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH} 

#install SiftGPU (CPU only, please refer to https://github.com/wangq95/SiftGPU_Linux for how to build with GPU)
RUN cd /usr/local && git clone https://github.com/pitzer/SiftGPU.git && cd SiftGPU && make 

# install ImageMatchsiftGPU, install to /usr/local (ImageMatchsiftGPU.zip already contains Eigen in the include folder)
COPY ImageMatchsiftGPU.zip .
RUN unzip ImageMatchsiftGPU.zip && cd ImageMatchsiftGPU/ImageMatchsiftGPU && cmake CMakeLists.txt && \
	make && \
	cd ../.. && rm -rf ImageMatchsiftGPU*

#not need to install pythton, we run python script outside the container

## Install miniconda
#ENV CONDA_DIR /usr/local/miniconda3
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/sh ~/miniconda.sh -b -p /usr/local/miniconda3  && rm ~/miniconda.sh
#
## Put conda in path so we can use conda activate
#ENV PATH=$CONDA_DIR/bin:$PATH
## need for pip install gdal
#ENV CPLUS_INCLUDE_PATH=/usr/local/include:${CPLUS_INCLUDE_PATH}
#ENV C_INCLUDE_PATH=/usr/local/include:${C_INCLUDE_PATH}
#
### Make RUN commands use the new environment (conda base): 
##SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]
#
## donwgrade setuptools<58, 
#RUN pip install setuptools==57.5.0 && pip install GDAL==2.0.0 
#
#
### t: tarball, i index-cache p: package  y: Do not ask for confirmation
##RUN conda clean -tipy
#
#RUN git clone https://github.com/yghlc/DeeplabforRS.git ${HOME}/codes/PycharmProjects/DeeplabforRS && \
#	git clone https://github.com/yghlc/rs_data_proc.git ${HOME}/codes/PycharmProjects/rs_data_proc
#
#ENV PATH=${HOME}/codes/PycharmProjects/rs_data_proc/coregistration_siftGPU:$PATH
#
#COPY run.sh .
## WORKDIR /data
#
#CMD ["conda", "run", "--no-capture-output", "-n", "base","/bin/bash","-c","./run.sh"]

