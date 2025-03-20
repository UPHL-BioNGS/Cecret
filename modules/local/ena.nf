process ENA_DOWNLOAD {
    tag           "${SRR}"
    label         "process_single"
    container     'staphb/enabrowsertools:1.7.1'
    
    input:
    val(SRR)

    output:
    tuple val(sra), file("sra/paired/${sra}*fastq.gz"), val(false), optional: true, emit: paired
    tuple val(sra), file("sra/single/${sra}*fastq.gz"), val(true),  optional: true, emit: single
    path "logs/*/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${SRR}"
    """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
        
    enaDataGet \
        ${args} \
        -f fastq \
        ${SRR} \
        | tee -a \$log_file

    if [ -f "\${sra}_2.fastq.gz" ]
    then
        mv *.fastq.gz sra/paired/.
    elif [ -f "\${sra}.fastq.gz" ]
    then
        mv *.fastq.gz sra/single/.
    else
        echo "Could not download file for SRA accession ${sra}"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        enaDataGet: \$( enaDataGet -v | awk '{print \$NF}' )
    END_VERSIONS
    """
}
