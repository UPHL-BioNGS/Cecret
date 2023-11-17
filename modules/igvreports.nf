process igv_reports {
  tag         "${sample}"
  label       "process_high"
  publishDir  path: "${params.outdir}", mode: 'copy'
  container   'quay.io/biocontainers/igv-reports:1.9.1--pyh7cba7a3_0'
  
  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.igv_reports

  input:
  tuple val(sample), file(bam), file(reference_genome), file(vcf)

  output:
    path "igv_reports/${sample}/igvjs_viewer.html"
    path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p igv_reports/!{sample} logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    create_report -v >> $log

    create_report !{params.igv_reports_options} \
        --fasta !{reference_genome} \
        --tracks !{vcf} !{bam} \
        --output igv_reports/!{sample}/igvjs_viewer.html \
        !{vcf} \
        | tee -a $log
  '''
}
