process aci {
    tag        "${sample}"
    label      "process_high"
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/uphl/aci:1.4.20240116-2024-01-17'
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}
  
    //#UPHLICA maxForks      10
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    when:
    params.aci

    input:
    tuple val(sample), file(bam), file(bed)

    output:
    path "aci/${sample}/${sample}_amplicon_depth.csv", emit: cov
    path "aci/${sample}/${sample}_amplicon_depth.png", emit: for_multiqc
    path "aci/${sample}/*"
    path "logs/${task.process}/aci.${workflow.sessionId}.log"
  
    shell:
    '''
        mkdir -p aci/!{sample} logs/!{task.process}
        log=logs/!{task.process}/aci.!{workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > $log
        aci --version >> $log

        aci \
            --bam !{bam} \
            --bed !{bed} \
            --threads !{task.cpus} \
            --out aci/!{sample} \
            | tee -a $log
        
        if [ -f "aci/!{sample}/amplicon_depth.csv" ] ; then cp aci/!{sample}/amplicon_depth.csv aci/!{sample}/!{sample}_amplicon_depth.csv ; fi
        if [ -f "aci/!{sample}/amplicon_depth.png" ] ; then cp aci/!{sample}/amplicon_depth.png aci/!{sample}/!{sample}_amplicon_depth.png ; fi

    '''
}
