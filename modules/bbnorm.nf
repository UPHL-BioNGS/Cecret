process bbnorm {
    tag           "${sample}"
    label         'process_medium'
    tag           "${sample}"
    publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    container     'staphb/bbtools:39.01'

    input:
    tuple val(sample), file(reads), val(paired_single)

    output:
    tuple val(sample), path("bbnorm/*.fastq.gz"), val(paired_single), emit: fastq, optional: true
    tuple val(sample), path("logs/*/*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    params.bbnorm && (task.ext.when == null || task.ext.when)

    shell:
    def args   = task.ext.args   ?: "${params.bbnorm_options}"
    def prefix = task.ext.prefix ?: "${sample}"
    if ( paired_single == 'paired' ) {
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
    } else if ( paired_single == 'single' ) {
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
    }
}

