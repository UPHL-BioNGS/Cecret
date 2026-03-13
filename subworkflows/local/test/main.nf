include { ENA                         } from '../../../modules/local/ena'
include { DATASETS as DOWNLOAD_GENOME } from '../../../modules/local/datasets'

workflow TEST {
    take:
    ch_sra_accessions // channel: list
    ch_genome_accessions // channel: list

    main:
    log.info """

Running public data retrieval. This workflow fetches raw sequencing reads (FASTQ) 
and reference genomes (FASTA) from public repositories like ENA and NCBI for 
testing or downstream analysis.

Relevant params and their values:
- 'params.sra_accessions' : ${params.sra_accessions}
    - Array of SRA accessions to download from the ENA
- 'params.genome_accessions' : ${params.genome_accessions}
    - Array of genome accessions to download from NCBI genomes

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ ENA                ┃ Downloads single or paired-end FASTQ files using SRA accessions.  ┃
┃ DATASETS           ┃ Downloads reference genome FASTA files using NCBI Datasets.       ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

    ch_reads   = channel.empty()
    ch_fastas  = channel.empty()
    ch_versions = channel.empty()

    // downloads fastq files from the ENA
    // note: works with GA, but not locally
    if ( ! params.sra_accessions.isEmpty() ) {
        ENA(ch_sra_accessions)

        ENA.out.paired
            .mix(ENA.out.single)
            .map { it ->
                def meta = [id:it[0], single_end:it[2]]
                tuple( meta, it[1]) }
            .set { ch_downloaded_reads }

        ch_reads = ch_reads.mix(ch_downloaded_reads)
        ch_versions = ch_versions.mix(ENA.out.versions)
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
