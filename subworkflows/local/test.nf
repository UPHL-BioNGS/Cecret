include { DOWNLOAD } from '../modules/local/local'

workflow TEST {
    take:
        ch_sra_accessions
        ch_genome_accessions

    main:
        DOWNLOAD_FASTQ(ch_sra_accessions)
        DONWLOAD_GENOME(ch_genome_accessions)

    //TODO: add meta to reads and fastas
    // TODO : add datasets to download genomes

    emit:
        reads    = ch_reads
        fastas   = ch_fastas
        versions = Channel.empty()
}
