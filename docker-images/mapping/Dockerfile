# DOCKERFILE to run mash_screen for pATLAS
#FROM ummidock/assemblerflow_base:1.0.0
FROM ubuntu:16.04
MAINTAINER Tiago F. Jesus, tiagojesus@medicina.ulisboa.pt

# INSTALL DEPENDENCIES
RUN apt-get update && apt-get -y install \
    python3 \
    wget \
    zip \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

WORKDIR /home/bin/

# Download and install bowtie2

RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip

RUN unzip bowtie2-2.3.4-linux-x86_64.zip

ENV PATH="/home/bin/bowtie2-2.3.4-linux-x86_64:$PATH"

# Download and install samtools

RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2

RUN tar jxf samtools-1.7.tar.bz2

WORKDIR /home/bin/samtools-1.7/

RUN ./configure --without-curses && make

ENV PATH="/home/bin/samtools-1.7:$PATH"

# Download databases from PlasmidCoverage

WORKDIR /home/data/

RUN wget https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/1.0.4/indexes.tar.gz

RUN tar xvzf indexes.tar.gz

# Fetches json with lenghts used to generate -k parameter for bowtie as well as
# plasmid lenghts
RUN wget https://raw.githubusercontent.com/tiagofilipe12/PlasmidCoverage/master/plasmiduncover/json/reads_sample_result_length.json