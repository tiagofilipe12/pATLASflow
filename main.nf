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
readInputs = Channel.empty()
readInputs2 = Channel.empty()
fastaInputs = Channel.empty()

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
    // this empties one of the channels in the case one is not required to run
    if (params.mash_screen == false) {
        readInputs = Channel.empty()
    }
    else if (params.mapping == false) {
        readInputs2 = Channel.empty()
    }
}
if (params.mash_screen || params.assembly) {
    refSketch = "/home/data/patlas.msh"
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
    //fetch indexes for mapping approach, available in Docker
    bowtie2Index = "/home/data/indexes/bowtie2idx/bowtie2.idx"  // idx_file
    samtoolsIndex = "/home/data/indexes/fasta/samtools.fasta.fai"   // maindb_path
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

    output:
    file "sortedMashScreenResults_${sample}.txt" into mashScreenResults

    """
    mash screen -i ${params.identity} -v ${params.pValue} -p \
    ${params.threads} ${winnerVar} ${refSketch} ${reads} > mashScreenResults_${sample}.txt
    sort -gr mashScreenResults_${sample}.txt > sortedMashScreenResults_${sample}.txt
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

    output:
    file "${fasta}_mashdist.txt" into mashDistResults

    """
    mash dist -p ${params.threads} -v ${params.pValue} \
    -d ${params.mash_distance} ${refSketch} ${fasta} > ${fasta}_mashdist.txt
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

/******************************/
/********** MAPPING ***********/
/******************************/

// process that runs bowtie2
process mappingBowtie {

    tag { "mapping sample: " + sample}

    input:
    set sample, file(reads) from readInputs2

    output:
    set sample, file("mappingBowtie_${sample}.sam") into bowtieResults

    script:

    if (params.singleEnd == true) {
        readsString = "-U ${reads}"
    }
    else {
        readsString = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    """
    bowtie2 -x ${bowtie2Index} ${readsString} -p ${params.threads} -k \
    ${params.max_k} -5 ${params.trim5} -S mappingBowtie_${sample}.sam
    """
}

/**
* samtools faidx is escaped because index file is already provided in docker
* image.
*/
process samtoolsView {
    tag { "samtools view: " +  sample }

    input:
    set sample, file(samtoolsFile) from bowtieResults

    output:
    set sample, file("samtoolsDepthOutput_${sample}.txt") into samtoolsResults

    """
    samtools view -b -t ${samtoolsIndex} -@ ${params.threads} ${samtoolsFile} | \
    samtools sort -@ ${params.threads} -o samtoolsSorted_${sample}.bam
    samtools index samtoolsSorted_${sample}.bam
    samtools depth samtoolsSorted_${sample}.bam > samtoolsDepthOutput_${sample}.txt
    """
}

// process to dump txt depth file to json