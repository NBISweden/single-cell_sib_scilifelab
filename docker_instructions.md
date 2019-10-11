# <img border="0" src="https://cdn1.iconfinder.com/data/icons/social-media-2106/24/social_media_social_media_logo_docker-512.png" width="40" height="40"> Docker

***

<br/>

In the last case scenario, if you are having problems installing R packages, please follow these instructions to setup a docker container with R and Rstudio. Then you will need to install packages as usual. See the [pre-course material](/single-cell_sib_scilifelab/precourse.md) for package installations.


<br/>

## <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="40" height="40"> Installation Instructions

***

1. Install docker. You need permission to both install and use `sudo` commands.

* Mac: https://docs.docker.com/docker-for-mac/install/
* Windows: https://docs.docker.com/docker-for-windows/install/
* Ubuntu: https://docs.docker.com/install/linux/docker-ce/ubuntu/

2. Launch docker and got to `Preferences > Advanced`  and set `Memory=8gb, cpu=4, swp=1`. These will set the maximum values possible for your machine.

3. Now go to `Preferences > Security` and tick on "Allow connections from network clients", so that you can communicate with DockerServer.

4. Open a terminal shell.

5. Create a Virtual Machine (VM) named "default" and set the amount of CPU and RAM available for you:

```bash
docker-machine create default
VBoxManage modifyvm default --cpus 4
VBoxManage modifyvm default --memory 8192
```

6. You can now start the machine, the environment (which will give you an IP address to connect to it!) and configure your shell to work with docker.

```bash
docker-machine start
docker-machine env
eval "$(docker-machine env default)"
```

If at any point from now on you get an error asking you to stop the container, you can do so with:

```bash
docker-machine stop
```

7. Now we can download a pre-made docker container containing both the latest version R and RStudio (rocker). Here we need to set a new password for the container, which we can set as "test". The `-p 8787:8787` indicates which ports are visible between your computer and the container for visualizing Rstudio. `-rm` will remove the container after use. `-e` sets the password and

```bash
docker run -d -rm --memory=8g -p 8787:8787 -e PASSWORD=test --name rstudio rocker/verse
```

If no errors are thrown, you will be able to connect to your Rstudio machine using a webbrowser, just type the IP.

* http://192.168.99.100:8787 (this might change, check the output from the command `docker-machine env` on the step above)
* LOGIN=`rstudio`
* PASS=`test`

You can alternativelly check your IP using:

```bash
docker-machine ip
```

8. You can now proceed with using Rstudio and installing packages as usual.

See the [pre-course material](/single-cell_sib_scilifelab/precourse.md)


<br/>


## <img border="0" src="https://cdn4.iconfinder.com/data/icons/proglyphs-computers-and-development/512/Terminal-512.png" width="40" height="40"> Useful commands

***

Some of these commands can be used in

```bash
#get container ip address
docker-machine ip

#list or remove docker containers
docker container ls
docker rm <container_id>

#list or remove docker images
docker image ls
docker rmi <image_id>

#runs a interactive (-ti) shell inside the container.
#"-v" mounts the current directory "$(pwd)" to the folder "/rstudio".
#You can use them separately
docker run -it -v $(pwd):/rstudio <docker_name>
```

<br/>

## [Back to main](README.md)
