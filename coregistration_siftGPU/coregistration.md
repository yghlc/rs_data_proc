## Using SIFTGPU for Co-registration of remote sensing images

I developed a program to co-registrator remote sensing images based on SIFTGPU in 2013. 
The software environment has changed a lot. 
Here, I built a docker container, based on Ubuntu 14.04, then install SIFTGPU and the program 
I developed. 

Update: The legacy opengl need display (X server) to run or CUDA to run, so the container does not work, 
we build to Virtual machine (VirtualBox), with Ubuntu 14.04, to run this. 



### Install SIFTGPU and ImageMatchsiftGPU
[SIFTGPU](https://github.com/pitzer/SiftGPU) is GPU implement of the [SIFT](https://link.springer.com/article/10.1023/B:VISI.0000029664.99615.94) 
algorithm.
We will build a docker container to install it. There is a GitHub repo about installing SiftGPU on Ubuntu 16.04: [link](https://github.com/wangq95/SiftGPU_Linux). 

For the entire installation, please see [ImageMatchsiftGPU](https://github.com/yghlc/ImageMatchsiftGPU) or
the "Dockerfile" in this folder.


### Running ImageMatchsiftGPU (Docker container)

```commandline
docker run --rm -w work_dir -u $UID:$UID  -v /data/LingcaoHuang:/data -v host_work_dir:work_dir \
-it siftgpu-image-match ImageMatchsiftGPU imagelist.txt result 2
```
Notes: 
* please change work_dir and "host_work_dir:work_dir" accordingly. 
"-w" is to set the working directory. If it's not set, the work_dir will be "/root" inside the container. "-v" will mount the working directory. 
* Make sure all the data are inside /data folder and is mounted. 
* "-u" is to set the user id and group id (user_id:group_id). On Ubuntu, $UID output current user id, and group id usually is the same as the user id. 
* We may set username "$(whoami)" to "-u", but may get error if this user is not create inside the docker container: "docker: Error response from daemon: unable to find user lingcao: no matching entries in passwd file." 

If "-u" is not set, the output file will be owned by "root". Please use chown to change the file ownership. 
```commandline
sudo chown NewUser:NewGroup filename
sudo chown -R NewOwner:Newgroup NameOfDirectory
```

add followings to allow docker container to used X server (displace) on the host machine. 
```commandline
--net=host -e DISPLAY=$DISPLAY  # not test yet
```


### TODO:
add GPU support in the docker Container.


    
    
