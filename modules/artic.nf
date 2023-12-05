process artic {
    tag        "${sample}"
    label      "process_high"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/uphl/artic:1.2.4-1.11.2-2023-12-05'
  
    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    input:
    tuple val(sample), file(fastq), file(reference), file(bed)

    output:
    tuple val(sample), file("artic/${sample}.primertrim.sorted.bam"), file("artic/${sample}.primertrim.sorted.bam.bai"), emit: bam, optional: true
    path "consensus/${sample}.consensus.fa", emit: consensus, optional: true
    tuple val("artic"), env(artic_version), emit: artic_version
    path "artic/${sample}*"
    path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  
    shell:
    '''
        mkdir -p artic consensus schema/cecret/V1 logs/!{task.process}
        log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > $log
        artic --version >> $log
        artic_version=$(artic --version)

        cp !{reference} schema/cecret/V1/cecret.reference.fasta
        cp !{bed}       schema/cecret/V1/cecret.scheme.bed
        samtools faidx  schema/cecret/V1/cecret.reference.fasta

        artic minion !{params.artic_options} \
            --threads !{task.cpus} \
            --read-file !{fastq} \
            --scheme-directory schema \
            --scheme-version 1 \
            cecret \
            artic/!{sample} \
            | tee -a $log

        if [ -f "artic/!{sample}.consensus.fasta" ]
        then
            echo ">!{sample}"                            > consensus/!{sample}.consensus.fa
            grep -v ">" artic/!{sample}.consensus.fasta >> consensus/!{sample}.consensus.fa
        fi

        if [ -f "artic/!{sample}.primertrimmed.rg.sorted.bam" ]
        then
            cp artic/!{sample}.primertrimmed.rg.sorted.bam artic/!{sample}.primertrim.sorted.bam
            samtools index artic/!{sample}.primertrim.sorted.bam
        fi
    '''
}

process artic_read_filtering {
    tag        "${sample}"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/uphl/artic:1.2.4-1.11.1'
    label      "process_single"
  
    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
    //#UPHLICA memory 1.GB
    //#UPHLICA cpus 3
    //#UPHLICA time '45m'

    input:
    tuple val(sample), file(fastq)

    output:
    tuple val(sample), file("artic/${sample}_filtered.fastq.gz"), emit: fastq
    path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  
    shell:
    '''
        mkdir -p artic logs/!{task.process}
        log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > $log
        artic --version >> $log

        artic guppyplex !{params.artic_read_filtering_options} \
            --directory . \
            --prefix !{fastq} \
            --output artic/!{sample}_filtered.fastq.gz \
            | tee -a $log
    '''
}
