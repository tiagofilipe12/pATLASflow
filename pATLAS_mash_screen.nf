#!/usr/bin/env nextflow

import Helper

// A pipeline to run mash screen for pATLAS.

if (params.help) {
    Help.print_help(params)
    exit 0
}

// check if noWinner is provided or not
winnerVar = (params.noWinner == false) ? "-w" : ""

readInputs = Channel.empty()
mashRef = Channel.empty()
fastaInputs = Channel.empty()
//TODO add channel for mapping index files

if (params.mash_screen) {
    // channel for mash screen inputs. Prints a line for each sample that match a
    //glob pattern
    readInputs = Channel
        .fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2, type: 'file')
        .ifEmpty { exit 0, "reads (fastq) files were not provided" }
    // channel for reference sketch file
    mashRef = Channel
       .fromPath(params.refSketch, type: 'file')
       .ifEmpty { exit 0, "no mash sketch"}
}

if (params.assembly) {
// channel for fasta inputs
    fastaInputs = Channel
        .fromPath(params.fasta, type: 'file')
        .ifEmpty {exit 0, "no fasta file was provided"}
    // channel for reference sketch file
    mashRef2 = Channel
       .fromPath(params.refSketch, type: 'file')
       .ifEmpty { exit 0, "no mash sketch"}
}

//todo implement checks on file type

if (params.mapping == true) {
    //todo something for mapping
    // bowtie2 index and samtools index channels
}

/******************************/
/*** MASH SCREEN processes ****/
/******************************/

// process to run mashScreen and sort the output into
// sortedMashScreenResults_{sampleId}.txt
process mashScreen {
    input:
    set sample, file(reads) from readInputs
    file(db) from mashRef

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
process mashOutputJson {
    input:
    file mashtxt from mashScreenResults

    script:
    template "mashscreen2json.py"
}

/******************************/
/**** MASH DIST processes *****/
/******************************/

process runMashDist {
    input:
    file fasta from fastaInputs
    file db from mashRef2

    output:
    file "${fasta}_mashdist.txt" into mashDistResults

    """
    mash dist -p ${params.threads} -v ${params.pValue} \
    -d ${params.mash_distance} ${db} ${fasta} >"${fasta}_mashdist.txt"
    """
}

process mashDistOutputJson {
    input:
    file mashtxt from mashDistResults

    script:
    template "mashdist2json.py"
}