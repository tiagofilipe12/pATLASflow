#!/usr/bin/env nextflow

import Helper

// A pipeline to run mash screen for pATLAS.

if (params.help){
    Help.print_help(params)
    exit 0
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