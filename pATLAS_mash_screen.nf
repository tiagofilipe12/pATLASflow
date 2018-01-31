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

// check if noWinner is provided or not
winnerVar = (params.noWinner == false) ? "-w" : ""

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

    """
    echo Processing ${sample}
    mash screen -i ${params.identity} -v ${params.pValue} -p \
    ${params.threads} ${winnerVar} ${db} ${reads} > mashScreenResults.txt
    sort -gr mashScreenResults.txt > sortedMashScreenResults_${sample}.txt
    """
}

// process to parse the output to json format
process mashOutputSort {
    input:
    file mashtxt from mashScreenResults

    script:
    template "mashscreen2json.py"
}