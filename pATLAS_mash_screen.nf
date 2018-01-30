#!/usr/bin/env nextflow

// A pipeline to run mash screen for pATLAS.

// if --help is used, print software usage
if(params.help) {
    log.info "\n----------------------------------------------------------------"
    log.info "mash screen wrapper to run mash screen and convert the output to" +
     " JSON format that can be interpreted by pATLAS (www.patlas.site)"
    log.info "----------------------------------------------------------------"
    log.info "Usage:"
    log.info "  --threads   Number of threads that mash screen will have to" +
         " run.    Default: 1"
    log.info "  --kMer  the length of the kmer to be used by mash.   Default: 21"
    log.info "  --pValue    The p-value cutoff. Default: 0.05"
    log.info "  --identity  The minimum identity value between two sequences." +
        " Default: 0.9"
    log.info "  --noWinner  This option allows to disable the -w option of " +
        "mash screen  Default: false"
    log.info "  --refSketch The file that has the reference mash screen used" +
        "by pATLAS  Default: reference/patlas.msh"
    log.info "  --reads The path to the read files. Here users may provide " +
        "many samples in the same directory. However be assured that glob " +
        "pattern is unique (e.g. 'path/to/*_{1,2}.fastq'). " +
        "Default: reads/*_{1,2}.fastq.gz"
    log.info "  --singleEnd Provide this option if you have single-end reads. " +
        "By default the pipeline will assume that you provide paired-end " +
        "reads.    Default: false"
    log.info "  --help  Opens this help. It will open only when --help is " +
        "provided. So, yes, this line is pretty useless since you already " +
        "know that if you reached here."
    log.info "----------------------------------------------------------------"
    exit 1
}

/* a set of parameters to provide to mash that can be passe through command line
 when using nextflow: nextflow run --threads 2 --kMer 32 --pValue 0.1 --mashDist
  0.05. here we define default values in the case parameters are not passed
  through console.
 */
params.threads = 1
params.kMer = 21
params.pValue = 0.05
params.identity = 0.9
// noWinner allows to disable the -w option for mash, this just have to be
//passed as an argument with nothing else
params.noWinner = false


// fetches reference sketch for pATLAS database
// This parameter can also be passed through shell command
params.refSketch = 'reference/patlas.msh'

/*
this wil search for every file pairs within the reads folder by default.
However users may use another folder but they should use the same pattern: *_{1,
2}.fastq.gz
*/
params.reads = 'reads/*_{1,2}.fastq.gz' // allows to import two files in reads directory
// sets singleEnd reads to false by default. users may want to pass this
//argument to turn it true
params.singleEnd = false

// channel for mash screen inputs. Prints a line for each sample that match a
//glob pattern
mashInputs = Channel
    .fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2, type: 'file')
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nPath" +
     " needs to be enclosed in quotes!\nIf this is single-end data, please " +
     "specify --singleEnd on the command line." }

// channel for reference sketch file
mashRef = Channel
    .fromPath(params.refSketch, type: 'file')


// process to run mashScreen and sort the output into
// sortedMashScreenResults_{sampleId}.txt
process mashScreen {
    input:
    set sample, file(reads) from mashInputs
    file db from mashRef

    output:
    file "sortedMashScreenResults_${sample}.txt" into mashScreenResults

    script:
        if( params.noWinner == false ) {
            """
            echo Processing ${sample}
            mash screen -i ${params.identity} -v ${params.pValue} -p \
            ${params.threads} -w ${db} ${reads} > mashScreenResults.txt
            sort -gr mashScreenResults.txt > sortedMashScreenResults_${sample}.txt
            """
        }
        else {
            """
            echo Processing ${sample}
            mash screen -i ${params.identity} -v ${params.pValue} -p \
            ${params.threads} ${db} ${reads} > mashScreenResults_${sample}.txt
            sort -gr mashScreenResults.txt > sortedMashScreenResults_${sample}.txt
            """
        }
}

// process to parse the output to json format
process mashOutputSort {
    input:
    file mashtxt from mashScreenResults

    script:
    template "mashscreen2json.py"
}