include { DOWNLOAD_FASTQ              } from '../../modules/local/local'
include { DATASETS as DOWNLOAD_GENOME } from '../../modules/local/datasets'

workflow TEST {
    take:
        ch_sra_accessions
        ch_genome_accessions

    main:
        DOWNLOAD_FASTQ(ch_sra_accessions)

        DOWNLOAD_FASTQ.out.paired
            .mix(DOWNLOAD_FASTQ.out.single)
            .map { it ->
                meta = [id:it[0], single_end:it[2]]
                tuple( meta, it[1]) }
            .set { ch_reads }

        DOWNLOAD_GENOME(ch_genome_accessions)

    emit:
        reads    = ch_reads
        fastas   = DOWNLOAD_GENOME.out.fasta
        versions = DOWNLOAD_GENOME.out.versions
}
