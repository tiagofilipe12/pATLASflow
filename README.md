# pATLASflow

[![](https://img.shields.io/badge/nextflow->=0.27.3-blue.svg)](#)
[![](https://img.shields.io/badge/docker_mash-ready-green.svg)](https://hub.docker.com/r/tiagofilipe12/patlasflow_mash_screen/)
[![](https://img.shields.io/badge/docker_mapping-ready-green.svg)](https://hub.docker.com/r/tiagofilipe12/patlasflow_mapping/)
[![](https://img.shields.io/badge/pATLAS-1.x.x-lightgrey.svg)](https://github.com/tiagofilipe12/pATLAS)

A pipeline to run mapping, mash screen and assembly methods for pATLAS.

# TOC

* [Brief description](#brief-description)
* [Requirements](#requirements)
    * [Conda recipe for nextflow](#conda-recipe-for-nextflow)
* [Usage](#usage)
* [Example run](#example-run)
* [TL;DR](#tldr)


## TL;DR

1. Read files must be placed in `<current working dir>/reads/` folder

2. Fasta files must be placed in `<current working dir>/fasta/` folder

3. Run the pipeline `nextflow run tiagofilipe12/pATLASflow` with the options you require:
    * Assembly: `nextflow run tiagofilipe12/pATLASflow --assembly`
    * Mapping: `nextflow run tiagofilipe12/pATLASflow --mapping`
    * Mash screen: `nextflow run tiagofilipe12/pATLASflow --mash_screen`

    Note: you can even run all approaches by doing:
    `nextflow run tiagofilipe12/pATLASflow --assembly --mapping --mash_screen`

## Brief description

This Nextflow script is an implementation of [mash-wrapper](https://github.com/tiagofilipe12/mash_wrapper#mash-screen-for-read-samples)
for mash screen module.
It will output a `JSON` file that can be imported into [pATLAS](http://www.patlas.site).

## Requirements

* [Java 8](http://www.oracle.com/technetwork/java/javase/downloads/index.html) or higher.
* [Docker](https://docs.docker.com/install/) or [Singularity](http://singularity.lbl.gov/install-linux).
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) - follow the two installation steps (there is no need to read anything else).

### Conda recipe for nextflow

If you prefer you can use this conda recipe for nextflow: [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/nextflow/README.html)


## Usage

```
Usage: nextflow run tiagofilipe12/pATLASflow [options] or nextflow run main.nf [options] or ./main.nf [options]

  Nextflow magic options:
       -profile    Forces nextflow to run with docker or singularity.   Default: docker     Choices: standard, singularity
   Main options:
       --help  Opens this help. It will open only when --help is provided. So, yes, this line is pretty useless since you already know that if you reached here.
       --version   Prints the version of the pipeline script.
       --threads   Number of threads that the pipeline will have to run.    Default: 1
       --mash_screen   Enables mash screen run.
       --assembly  Enables mash dist run to use fasta file against plasmid db
       --mapping   Enables mapping pipeline.
   Mash options:
       --kMer  the length of the kmer to be used by mash.   Default: 21
       --pValue    The p-value cutoff. Default: 0.05
   Mash screen exclusive options:
       --identity  The minimum identity value between two sequences. Default: 0.9
       --noWinner  This option allows to disable the -w option of mash screen  Default: false
   Mash dist exclusive options:
       --mash_distance     Provide the maximum distance between two plasmids to be reported.   Default: 0.1
   Reads options:
       --reads The path to the read files. Here users may provide many samples in the same directory. However be assured that glob pattern is unique (e.g. 'path/to/*_{1,2}.fastq').
       --singleEnd Provide this option if you have single-end reads. By default the pipeline will assume that you provide paired-end reads.    Default: false
   Fasta options:
       --fasta     Provide fasta file pattern to be searched by nextflow.  Default: 'fasta/*.fas'
   Bowtie2 options:
       --max_k     Provide the maximum number of alignments allowed per read.  Default: 10949 (the number of plasmids present in pATLAS)
       --trim5     Provide parameter -5 to bowtie2 allowing to trim 5' end.    Default: 0
       --cov_cutoff    Provide a cutoff value to filter results for coverage results.  Default: 0.60

```


## Example run

`nextflow run tiagofilipe12/pATLASflow --assembly`

