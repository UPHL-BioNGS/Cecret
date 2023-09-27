process aci {
    tag        "Graphing amplicon depths"
    label      "process_high"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/erinyoung/aci:0.1.20230815'
    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
    //#UPHLICA maxForks      10
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    when:
    params.aci

    input:
    tuple file(bam), file(bed)

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
        aci --version >> $log

        aci \
            --bam !{bam} \
            --bed !{bed} \
            --threads !{task.cpus} \
            --out aci

        cp aci/amplicon_depth.png aci/amplicon_depth_mqc.png
    '''
}
