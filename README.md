# pATLAS_mash_screen.nf

A pipeline to run mash screen for pATLAS.

## Brief description

This Nextflow script is an implementation of [mash-wrapper](https://github.com/tiagofilipe12/mash_wrapper#mash-screen-for-read-samples)
for mash screen module.
It will output a `JSON` file that can be imported into [pATLAS](http://www.patlas.site).

## Usage

```
--threads       Number of threads that mash screen will have to run.   
    Default: 1
--kMer      the length of the kmer to be used by mash.   Default: 21
--pValue        The p-value cutoff. Default: 0.05
--identity      The minimum identity value between two sequences. Default: 0.9
--noWinner      This option allows to disable the -w option of mash screen  
    Default: false
--refSketch     The file that has the reference mash screen usedby pATLAS  
    Default: reference/patlas.msh
--reads     The path to the read files. Here users may provide many samples in 
    the same directory. However be assured that glob pattern is unique 
    (e.g. 'path/to/*_{1,2}.fastq'). Default: reads/*_{1,2}.fastq.gz
--singleEnd     Provide this option if you have single-end reads. By default 
    the pipeline will assume that you provide paired-end reads.    Default: false
--help      Opens this help. It will open only when --help is provided. So, 
    yes, this line is pretty useless since you already know that if you 
    reached here.
```

## TODO

* Default inputs for `reads` and `refSketch` are still sketchy and `patlas.msh`
should be downloaded [here](https://github.com/tiagofilipe12/mash_wrapper/releases/download/1.0.5/patlas.msh) for now.