# pATLAS_mash_screen.nf


[![](https://img.shields.io/badge/nextflow->=0.27.3-blue.svg)](#)
[![](https://img.shields.io/badge/docker-ready-green.svg)](https://hub.docker.com/r/tiagofilipe12/patlas_mash_screen/)

A pipeline to run mash screen for pATLAS.

## Brief description

This Nextflow script is an implementation of [mash-wrapper](https://github.com/tiagofilipe12/mash_wrapper#mash-screen-for-read-samples)
for mash screen module.
It will output a `JSON` file that can be imported into [pATLAS](http://www.patlas.site).

## Usage

```
Usage:
   nextflow run tiagofilipe12/pATLAS_mash_screen.nf

   Nextflow magic options:
       -profile    Forces nextflow to run with docker or singularity.   Default: docker     Choices: standard, singularity
   Main options:
       --help  Opens this help. It will open only when --help is provided. So, yes, this line is pretty useless since you already know that if you reached here.
       --version   Prints the version of the pipeline script.
       --threads   Number of threads that mash screen will have to run.    Default: 1
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
```


## Example run

`nextflow run tiagofilipe12/pATLAS_auxiliary_scripts --assembly`
