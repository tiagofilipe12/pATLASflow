// Help class courtesy of Diogo N. Silva (https://github.com/ODiogoSilva/innuca-nf)

class Help{

    static def print_help(params) {
        println("\n===========================================================")
        println("                   p A T L A S f l o w")
        println("===========================================================\n")
        println("Version: 1.1.0\n")
        println("Usage:")
        println("   nextflow run tiagofilipe12/pATLAS_mash_screen.nf\n")
        println("   Nextflow magic options:")
        println("       -profile    Forces nextflow to run with docker or " +
            "singularity.   Default: standard     Choices: standard, singularity," +
            "slurm")
        println("   Main options:")
        println("       --help  Opens this help. It will open only when --help is " +
            "provided. So, yes, this line is pretty useless since you already " +
            "know that if you reached here.")
        println("       --version   Prints the version of the pipeline script.")
        println("       --threads   Number of threads that the pipeline will have to" +
            " run.    Default: 1")
        println("       --mash_screen   Enables mash screen run.")
        println("       --assembly  Enables mash dist run to use fasta file against " +
            "plasmid db")
        println("       --mapping   Enables mapping pipeline.")
        println("   Mash options:")
        println("       --kMer  the length of the kmer to be used by mash.   Default: 21")
        println("       --pValue    The p-value cutoff. Default: 0.05")
        println("   Mash screen exclusive options:")
        println("       --identity  The minimum identity value between two sequences." +
            " Default: 0.9")
        println("       --noWinner  This option allows to disable the -w option of " +
            "mash screen  Default: false")
        println("   Mash dist exclusive options:")
        println("       --mash_distance     Provide the maximum distance between " +
            "two plasmids to be reported.   Default: 0.1")
        println("       --shared_hashes     Provide a percentage for the hashes" +
            "shared between the reference and the query sequence(s).   Default: 0.8")
        println("   Reads options:")
        println("       --reads The path to the read files. Here users may provide " +
            "many samples in the same directory. However be assured that glob " +
            "pattern is unique (e.g. 'path/to/*_{1,2}.fastq').")
        println("       --singleEnd Provide this option if you have single-end reads. " +
            "By default the pipeline will assume that you provide paired-end " +
            "reads.    Default: false")
        println("   Fasta options:")
        println("       --fasta     Provide fasta file pattern to be searched by " +
            "nextflow.  Default: 'fasta/*.fas'")
        println("   Bowtie2 options:")
        println("       --max_k     Provide the maximum number of alignments allowed " +
            "per read.  Default: 10949 (the number of plasmids present in pATLAS)")
        println("       --trim5     Provide parameter -5 to bowtie2 allowing to trim " +
            "5' end.    Default: 0")
        println("       --cov_cutoff    Provide a cutoff value to filter results for " +
            "coverage results.  Default: 0.60")
    }

}