#DOCKERFILE to run mash_screen for pATLAS
FROM ubuntu:16.04
MAINTAINER Tiago F. Jesus, tiagojesus@medicina.ulisboa.pt

#INSTALL DEPENDENCIES
RUN apt-get update && apt-get -y install \
    python3\
    wget

WORKDIR /home/bin

RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar
RUN tar -xvf mash-Linux64-v2.0.tar

#adds mash to path so that it can be executed in any pwd
ENV PATH="/home/bin/mash-Linux64-v2.0:$PATH"

WORKDIR /home/data