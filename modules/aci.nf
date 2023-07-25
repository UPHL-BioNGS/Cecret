process aci {
    tag        "Graphing ampicon depths"
    label      "maxcpus"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/erinyoung/aci:0.0.20230715'
  
    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
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
        echo !{bam}
        echo !{bed}

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
