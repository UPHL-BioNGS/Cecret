process artic {
    tag        "${sample}"
    label      "process_high"
    publishDir params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    container  'quay.io/uphl/artic:1.2.4-1.11.3-2023-12-19'

    //#UPHLICA maxForks      10
    //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
    //#UPHLICA memory 60.GB
    //#UPHLICA cpus 14
    //#UPHLICA time '45m'

    when:
    task.ext.when == null || task.ext.when

    input:
    tuple val(sample), file(fastq), file(reference), file(bed)

    output:
    tuple val(sample), file("artic/${sample}.primertrim.sorted.bam"), file("artic/${sample}.primertrim.sorted.bam.bai"), emit: bam, optional: true
    path "consensus/${sample}.consensus.fa", emit: consensus, optional: true
    tuple val("artic"), env(artic_version), emit: artic_version
    path "artic/${sample}*"
    path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
    path "versions.yml", emit: versions
  
    shell:
    def args   = task.ext.args   ?: "${params.artic_options}}"
    def prefix = task.ext.prefix ?: "${sample}"
    """
    mkdir -p artic consensus schema/cecret/V1 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    artic --version >> \$log
    artic_version=\$(artic --version)

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

process artic_read_filtering {
    tag        "${sample}"
    publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    container  'quay.io/uphl/artic:1.2.4-1.11.3-2023-12-19'
    label      "process_single"
  
    //#UPHLICA maxForks      10
    //#UPHLICA \\errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
    //#UPHLICA memory 1.GB
    //#UPHLICA cpus 3
    //#UPHLICA time '45m'

    input:
    tuple val(sample), file(fastq)

    output:
    tuple val(sample), file("artic/${sample}_filtered.fastq.gz"), emit: fastq
    path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
    path "versions.yml", emit: versions
  
    shell:
    def args   = task.ext.args   ?: "${params.artic_read_filtering_options}"
    def prefix = task.ext.prefix ?: "${sample}"
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
