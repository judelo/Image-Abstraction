# use one of the images from this repository: https://github.com/centreborelli/ipol-docker-images/
FROM debian:buster

# install additional debian packages
COPY .ipol/packages.txt packages.txt
RUN apt-get update && apt-get install -y $(cat packages.txt) && rm -rf /var/lib/apt/lists/* && rm packages.txt


# copy the code to $bin
ENV bin /workdir/bin/
RUN mkdir -p $bin
WORKDIR $bin
COPY . .

# compile code
RUN qmake -makefile . && make
ENV ex /workdir/exec
COPY image_abstraction ex

# the execution will happen in the folder /workdir/exec
# it will be created by IPOL

# $HOME is writable by the user `ipol`, but 
ENV HOME /home/ipol
# chmod 777 so that any user can use the HOME, in case the docker is run with -u 1001:1001
RUN groupadd -g 1000 ipol && useradd -m -u 1000 -g 1000 ipol -d $HOME && chmod -R 777 $HOME
USER ipol
