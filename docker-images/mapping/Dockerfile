# DOCKERFILE to run mash_screen for pATLAS
FROM ubuntu:16.04
MAINTAINER Tiago F. Jesus, tiagojesus@medicina.ulisboa.pt

# INSTALL DEPENDENCIES
RUN apt-get update && apt-get -y install \
    python3\
    wget\
    zip

WORKDIR /home/bin/

# Download and install bowtie2

RUN wget -O bowtie2-2.3.4-linux-x86_64.zip https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.3.4%2Fbowtie2-2.3.4-linux-x86_64.zip%2Fdownload&ts=1517846768

RUN unzip bowtie2-2.3.4-linux-x86_64.zip

ENV PATH="/home/bin/bowtie2-2.3.4-linux-x86_64:$PATH"

# Download and install samtools

RUN wget -O samtools.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2/download

RUN tar jxf samtools.tar.bz2

WORKDIR /home/bin/samtools/

RUN ./configure && make

RUN PATH="/home/bin/samtools:$PATH"

# Download databases from PlasmidCoverage

WORKDIR /home/data/

RUN wget https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/1.0.4/indexes.tar.gz

RUN tar xvzf indexes.tar.gz

# I think that this is it but needs still needs testing