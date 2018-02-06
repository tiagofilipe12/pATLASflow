# DOCKERFILE to run mash_screen for pATLAS
FROM ubuntu:16.04
MAINTAINER Tiago F. Jesus, tiagojesus@medicina.ulisboa.pt

# INSTALL DEPENDENCIES
RUN apt-get update && apt-get -y install \
    python3\
    wget\
    zip\
    bzip2\
    build-essential

WORKDIR /home/bin/

# Download and install bowtie2

RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip

RUN unzip bowtie2-2.3.4-linux-x86_64.zip

ENV PATH="/home/bin/bowtie2-2.3.4-linux-x86_64:$PATH"

# Download and install samtools

RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2

RUN tar jxf samtools-1.7.tar.bz2

WORKDIR /home/bin/samtools/

RUN ./configure && make

RUN PATH="/home/bin/samtools:$PATH"

# Download databases from PlasmidCoverage

WORKDIR /home/data/

RUN wget https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/1.0.4/indexes.tar.gz

RUN tar xvzf indexes.tar.gz

# I think that this is it but needs still needs testing