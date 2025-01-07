process BBNORM {
    tag           "${meta.id}"
    label         'process_medium'
    container     'staphb/bbtools:39.13'

    input:
    tuple val(meta), file(reads)

    output:
    tuple val(meta), path("bbnorm/*.fastq.gz"), emit: fastq, optional: true
    path("logs/*/*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: "${params.bbnorm_options}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( meta.single_reads ) {
        """
        mkdir -p bbnorm logs/${task.process}
        log=logs/${task.process}/${prefix}.${workflow.sessionId}.log
            
        bbnorm.sh ${args} \
            in=${reads} \
            out=bbnorm/${prefix}_norm.fastq.gz \
            threads=${task.cpus} \
            | tee -a \$log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbnorm: \$(bbnorm.sh --version 2>&1 | grep version | tail -n 1 | awk '{print \$NF}')
            container: ${task.container}
        END_VERSIONS
        """
    } else {
        """
        mkdir -p bbnorm logs/${task.process}
        log=logs/${task.process}/${prefix}.${workflow.sessionId}.log
        
        bbnorm.sh ${args} \
            in=${reads[0]} \
            in2=${reads[1]} \
            out=bbnorm/${prefix}_norm_R1.fastq.gz \
            out2=bbnorm/${prefix}_norm_R2.fastq.gz \
            threads=${task.cpus} \
            | tee -a \$log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbnorm: \$(bbnorm.sh --version 2>&1 | grep version | tail -n 1 | awk '{print \$NF}')
            container: ${task.container}
        END_VERSIONS
        """
    }
}

