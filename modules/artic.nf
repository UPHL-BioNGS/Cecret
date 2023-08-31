process artic {
    tag        "${sample}"
    label      "maxcpus"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/artic:1.2.3--pyhdfd78af_0'
  
    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    input:
    tuple val(sample), file(fastq), file(bed), file(reference)

    output:
    path "aci/amplicon_depth.csv",     emit: cov
    path "aci/amplicon_depth.png"
    path "aci/amplicon_depth_mqc.png", emit: for_multiqc
    tuple val("artic"), env(artic_version),                 emit: aligner_version
    path "logs/${task.process}/aci.${workflow.sessionId}.log"
  
    shell:
    '''
        mkdir -p logs/!{task.process}
        log=logs/!{task.process}/aci.!{workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > $log
        artic --version >> $log
        artic_version=$(artic --version)

        cp !{reference} cecret.fasta
        cp !{bed} cecret.bed

        artic minion ${params.artic_options} \
            --threads !{task.cpus} \
            --read-file !{fastq} \
            cecret \
            !{sample}

        exit 1
    '''
}

process artic_read_filtering {
    tag        "${sample}"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/artic:1.2.3--pyhdfd78af_0'
  
    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    input:
    tuple val(sample), file(fastq)

    output:
    path "aci/amplicon_depth.csv",     emit: cov
    path "aci/amplicon_depth.png"
    path "aci/amplicon_depth_mqc.png", emit: for_multiqc
    path "logs/${task.process}/aci.${workflow.sessionId}.log"
  
    shell:
    '''
        mkdir -p logs/!{task.process}
        log=logs/!{task.process}/aci.!{workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > $log
        artic --version >> $log

        artic guppyplex ${params.artic_read_filtering_options} \
            --directory . \
            --prefix !{fastq} \
            --output !{sample}_filtered.fastq.gz

        exit 1
    '''
}
