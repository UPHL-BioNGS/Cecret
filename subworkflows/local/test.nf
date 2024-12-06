include { DOWNLOAD_FASTQ              } from '../../modules/local/local'
include { DATASETS as DOWNLOAD_GENOME } from '../../modules/local/datasets'

workflow TEST {
    take:
    ch_sra_accessions // channel: list
    ch_genome_accessions // channel: list

    main:
    ch_reads   = Channel.empty()
    ch_fastas  = Channel.empty()
    ch_versions = Channel.empty()

    // downloads fastq files from the ENA
    // note: works with GA, but not locally
    if ( ! params.sra_accessions.isEmpty() ) {
        DOWNLOAD_FASTQ(ch_sra_accessions)

        DOWNLOAD_FASTQ.out.paired
            .mix(DOWNLOAD_FASTQ.out.single)
            .map { it ->
                meta = [id:it[0], single_end:it[2]]
                tuple( meta, it[1]) }
            .set { ch_downloaded_reads }

        ch_reads = ch_reads.mix(ch_downloaded_reads)
    }

    // downloads genomes with datasets
    if ( ! params.genome_accessions.isEmpty() ) {
        DOWNLOAD_GENOME(ch_genome_accessions)
        ch_fastas = ch_fastas.mix(DOWNLOAD_GENOME.out.fasta)
        ch_versions = ch_versions.mix(DOWNLOAD_GENOME.out.versions)
    }

    emit:
        reads    = ch_reads // channel: [meta, reads]
        fastas   = ch_fastas // channel: fasta
        versions = ch_versions // channel: values
}
