process bbnorm {
    tag           "${sample}"
    label         'process_medium'
    tag           "${sample}"
    publishDir    "${params.outdir}", mode: 'copy'
    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    container     'staphb/bbtools:39.01'

    input:
    tuple val(sample), file(reads), val(paired_single)

    output:
    tuple val(sample), path("bbnorm/*.fastq.gz"), val(paired_single), emit: fastq
    tuple val(sample), path("logs/*/*.log"), emit: log

    when:
    params.bbnorm

    shell:
    if ( paired_single == 'paired' ) {
        '''
        mkdir -p bbnorm logs/!{task.process}
        log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        
        bbnorm.sh !{params.bbnorm_options} \
            in=!{reads[0]} \
            in2=!{reads[1]} \
            out=bbnorm/!{sample}_norm_R1.fastq.gz \
            out2=bbnorm/!{sample}_norm_R2.fastq.gz \
            threads=!{task.cpus} \
            | tee -a $log
        '''
    } else if ( paired_single == 'single' ) {
        '''
        mkdir -p bbnorm logs/!{task.process}
        log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        
        bbnorm.sh ${params.bbnorm_options} \
            in=!{reads} \
            out=bbnorm/!{sample}_norm.fastq.gz \
            threads=!{task.cpus} \
            | tee -a $log
        '''
    }
}

