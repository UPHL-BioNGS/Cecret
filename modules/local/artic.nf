process ARTIC {
    tag        "${meta.id}"
    label      "process_high"
    container  'staphb/artic:1.2.4-1.12.0'

    input:
    tuple val(meta), file(fastq), file(reference), file(bed)

    output:
    tuple val(meta), file("artic/*.primertrim.sorted.bam"), file("artic/*.primertrim.sorted.bam.bai"), emit: bam, optional: true
    path "consensus/*.consensus.fa", emit: consensus, optional: true
    tuple val("artic"), env(artic_version), emit: artic_version
    path "artic/*", emit: files
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions
  
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: "${params.artic_options}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p artic consensus schema/cecret/V1 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    artic --version >> \$log
    artic_version=\$(artic --version | awk '{print \$NF}')

    cp ${reference} schema/cecret/V1/cecret.reference.fasta
    cp ${bed}       schema/cecret/V1/cecret.scheme.bed
    samtools faidx  schema/cecret/V1/cecret.reference.fasta

    artic minion ${args} \
        --threads ${task.cpus} \
        --read-file ${fastq} \
        --scheme-directory schema \
        --scheme-version 1 \
        cecret \
        artic/${prefix} \
        | tee -a \$log

    if [ -f "artic/${prefix}.consensus.fasta" ]
    then
        echo ">${prefix}"                            > consensus/${prefix}.consensus.fa
        grep -v ">" artic/${prefix}.consensus.fasta >> consensus/${prefix}.consensus.fa
    fi

    if [ -f "artic/${prefix}.primertrimmed.rg.sorted.bam" ]
    then
        cp artic/${prefix}.primertrimmed.rg.sorted.bam artic/${prefix}.primertrim.sorted.bam
        samtools index artic/${prefix}.primertrim.sorted.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic --version | awk '{print \$NF}')
        medaka: \$( medaka --version | awk '{print \$NF}')
        container: ${task.container}
    END_VERSIONS
    """
}

process ARTIC_FILTER {
    tag        "${meta.id}"
    container  'staphb/artic:1.2.4-1.12.0'
    label      "process_single"

    input:
    tuple val(meta), file(fastq)

    output:
    tuple val(meta), file("artic/*_filtered.fastq.gz"), emit: fastq
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
            --output artic/${prefix}_filtered.fastq.gz \
            | tee -a \$log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(artic --version | awk '{print \$NF}')
            container: ${task.container}
        END_VERSIONS
    """
}
