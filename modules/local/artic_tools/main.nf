process ARTIC_TOOLS {
    tag        "${bed}"
    container  'staphb/artic-tools:0.3.1'
    label      "process_low"

    input:
    file(bed)

    output:
    path "artic-tools/*_insert.bed", emit: insert_bed
    path bed, emit: bed
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: bed.baseName
    """
    mkdir -p artic-tools logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    artic-tools validate_scheme \
        ${args} \
        ${bed} \
        --outputInserts artic-tools/${prefix}_insert.bed |
        tee -a >> \log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic-tools: \$(artic-tools --version | awk '{print \$NF}')
        container: ${task.container}
    END_VERSIONS
    """
}
