#!/usr/bin/env nextflow

import Helper

// A pipeline to run mash screen for pATLAS.

if (params.help) {
    Help.print_help(params)
    exit 0
}

// check if noWinner is provided or not
winnerVar = (params.noWinner == false) ? "-w" : ""

// start all optional channels
readInputs = readInputs2 = mashRef = mashRef2 = fastaInputs = Channel.empty()

/**
* The combination of these four channels (two forked) allows to double check for
* the run mode. If mash_screen is provided, then it will search for the two
* channel types. However if mapping is provided but mash_screen don't, it will
* search for the read inputs ignoring mash sketch.
*/
if (params.mash_screen || params.mapping) {
    Channel.fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2, type: 'file')
        .ifEmpty { exit 0, "reads (fastq) files were not provided" }
        .into { readInputs; readInputs2 }
}
if (params.mash_screen || params.assembly) {
    Channel.fromPath(params.refSketch)
        .ifEmpty { exit 0, "no mash sketch"}
        .into { mashRef; mashRef2 }
}

/**
* This channel enables assembly approach, allowing users to import a fasta to
* provide to mash dist
*/
if (params.assembly) {
// channel for fasta inputs
    fastaInputs = Channel
        .fromPath(params.fasta, type: 'file')
        .ifEmpty {exit 0, "no fasta file was provided"}
}

//TODO add channel for mapping index files

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

    tag { "running mash screen for sample: " + sample }

    input:
    set sample, file(reads) from readInputs
    file db from mashRef

    output:
    file "sortedMashScreenResults_${sample}.txt" into mashScreenResults

    """
    mash screen -i ${params.identity} -v ${params.pValue} -p \
    ${params.threads} ${winnerVar} ${db} ${reads} > mashScreenResults.txt
    sort -gr mashScreenResults.txt > sortedMashScreenResults_${sample}.txt
    """
}

// process to parse the output to json format
process mashOutputJson {

    tag { "dumping json file from: " + mashtxt }

    input:
    file mashtxt from mashScreenResults

    script:
    template "mashscreen2json.py"
}

/******************************/
/**** MASH DIST processes *****/
/******************************/

// runs mash dist
process runMashDist {

    tag { "running mash dist for fasta file: " + fasta }

    input:
    file fasta from fastaInputs
    file db from mashRef2

    output:
    file "${fasta}_mashdist.txt" into mashDistResults

    """
    mash dist -p ${params.threads} -v ${params.pValue} \
    -d ${params.mash_distance} ${db} ${fasta} > ${fasta}_mashdist.txt
    """
}

// parses mash dist output to a json file that can be imported into pATLAS
process mashDistOutputJson {

    tag { "dumping json file from: " + mashtxt }

    input:
    file mashtxt from mashDistResults

    script:
    template "mashdist2json.py"
}