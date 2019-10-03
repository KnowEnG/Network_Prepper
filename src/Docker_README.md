# Building the Network Prepper Docker Image
The Dockefile in this directory contains all the commands, in order, needed to build the **Network Prepper** docker image.


* Run the "make" command to build the **Network Prepper** docker image (output: docker image called "network_prepper" and a tag with today's date and time):
```
    make build_docker_image
```

* Login to docker hub. When prompted, enter your password and press enter:
```
    make login_to_dockerhub username=(enter your docker login here) email=(enter your email here)
```

* Upload your image to docker hub:
```
    make push_to_dockerhub
```

* * * 
## How to run this docker image
* * * 
### 1. Check on docker.hub to get the latest image tag: 10_03_2019 used here.

### 2. Change directory to the directory where you want to run and start the container.
```
docker run -v \`pwd\`:/home/test/run_dir/ -it knowengdev/network_prepper:10_03_2019 
```
