process ARTIC_FILTER {
    tag        "${meta.id}"
    container  'staphb/artic:1.10.0'
    label      "process_low"

    input:
    tuple val(meta), file(fastq)

    output:
    tuple val(meta), file("artic/*_filtered.fastq"), emit: fastq
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: "${params.artic_read_filtering_options}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
        mkdir -p artic logs/${task.process}
        log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > \$log
        artic --version >> \$log

        artic guppyplex ${args} \
            --directory . \
            --prefix ${fastq} \
            --output artic/${prefix}_filtered.fastq \
            | tee -a \$log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(artic --version | awk '{print \$NF}')
            container: ${task.container}
        END_VERSIONS
    """
}
