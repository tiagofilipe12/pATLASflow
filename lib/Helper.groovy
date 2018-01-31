// Help class cortesy of Diogo N. Silva (https://github.com/ODiogoSilva/innuca-nf)

class Help{

    static def print_help(params) {
        println("\n===========================================================")
        println("           p A T L A S auxiliary pipelines")
        println("===========================================================\n")
        println("Usage:")
        println("   nextflow run tiagofilipe12/pATLAS_mash_screen.nf\n")
        println("   Main options:")
        println("       --help  Opens this help. It will open only when --help is " +
            "provided. So, yes, this line is pretty useless since you already " +
            "know that if you reached here.")
        println("       --version   Prints the version of the pipeline script.")
        println("       --threads   Number of threads that mash screen will have to" +
            " run.    Default: 1")
        println("   Mash options:")
        println("       --kMer  the length of the kmer to be used by mash.   Default: 21")
        println("       --pValue    The p-value cutoff. Default: 0.05")
        println("       --refSketch The file that has the reference mash screen used" +
            " by pATLAS  Default: reference/patlas.msh")
        println("   Mash screen options")
        println("       --identity  The minimum identity value between two sequences." +
            " Default: 0.9")
        println("       --noWinner  This option allows to disable the -w option of " +
            "mash screen  Default: false")
        println("   Reads options:")
        println("       --reads The path to the read files. Here users may provide " +
            "many samples in the same directory. However be assured that glob " +
            "pattern is unique (e.g. 'path/to/*_{1,2}.fastq').")
        println("       --singleEnd Provide this option if you have single-end reads. " +
            "By default the pipeline will assume that you provide paired-end " +
            "reads.    Default: false")
    }

}