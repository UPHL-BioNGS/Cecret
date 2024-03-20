process samtools_stats {
  tag        "${sample}"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_stats && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_stats/${sample}.stats.txt", emit: samtools_stats_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_stats_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_stats logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools stats ${args} ${bam} > samtools_stats/${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_coverage {
  tag        "${sample}"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_coverage && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_coverage/${sample}.cov.{txt,hist}", emit: files
  path "samtools_coverage/${sample}.cov.txt", emit: samtools_coverage
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_coverage_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_coverage logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools coverage ${args} ${bam} -m -o samtools_coverage/${prefix}.cov.hist | tee -a \$log
    samtools coverage ${args} ${bam} | awk -v sample=${prefix} '{print sample "\\t" \$0 }' | sed '0,/${prefix}/s//sample/' > samtools_coverage/${prefix}.cov.txt  | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_flagstat {
  tag        "${sample}"
  label      "process_single"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_flagstat && (task.ext.when == null || task.ext.when)

  output:
  path "samtools_flagstat/${sample}.flagstat.txt", emit: samtools_flagstat_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_flagstat_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_flagstat logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools flagstat ${args} \
      ${bam} | \
      tee samtools_flagstat/${prefix}.flagstat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_depth {
  tag        "${sample}"
  label      "process_single"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_depth && (task.ext.when == null || task.ext.when)

  output:
  path "samtools_depth/${sample}.depth.txt", emit: file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_depth_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_depth logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools depth ${args} \
      ${bam} > samtools_depth/${prefix}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_ampliconstats {
  tag        "${sample}"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_ampliconstats && ( params.trimmer != 'none' ) && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("samtools_ampliconstats/${sample}_ampliconstats.txt"), emit: samtools_ampliconstats_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_ampliconstats_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_ampliconstats logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools ampliconstats ${args} \
      ${primer_bed} \
      ${bam} > samtools_ampliconstats/${prefix}_ampliconstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_plot_ampliconstats {
  tag           "${sample}"
  label         "process_single"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/samtools:1.19'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_plot_ampliconstats && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(ampliconstats)

  output:
  path "samtools_plot_ampliconstats/${sample}*", emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_plot_ampliconstats_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p samtools_plot_ampliconstats/${prefix} logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    plot-ampliconstats ${args} \
      samtools_plot_ampliconstats/${prefix} \
      ${ampliconstats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_sort {
  tag        "${sample}"
  label      "process_high"
  publishDir  path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p aligned logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -@ ${task.cpus} ${sam} | \
      samtools view -F 4 -o aligned/${prefix}.sorted.bam | tee -a \$log

    # indexing the bams
    samtools index aligned/${prefix}.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_filter {
  tag        "${sample}"
  label      "process_single"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.filter && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("filter/${sample}_filtered_{R1,R2}.fastq.gz"), optional: true, emit: filtered_reads
  path "filter/${sample}_filtered_unpaired.fastq.gz", optional: true, emit: unpaired
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
    def args   = task.ext.args   ?: '-F 4'
    def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p filter logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -n ${sam} | \
      samtools fastq ${args} ${params.filter_options} \
      -s filter/${prefix}_filtered_unpaired.fastq.gz \
      -1 filter/${prefix}_filtered_R1.fastq.gz \
      -2 filter/${prefix}_filtered_R2.fastq.gz \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_ampliconclip {
  tag        "${sample}"
  label      "process_single"
  publishDir  path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("ampliconclip/${sample}.primertrim.sorted.bam"), file("ampliconclip/${sample}.primertrim.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"                                                         
  tuple val("samtools ampliconclip"), env(trimmer_version), emit: trimmer_version
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.samtools_ampliconclip_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p ampliconclip logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log
    trimmer_version="samtools ampliconclip : \$(samtools --version | head -n 1)"

    # trimming the reads
    samtools ampliconclip ${args} -b ${primer_bed} ${bam} | \
      samtools sort |  \
      samtools view -F 4 -o ampliconclip/${prefix}.primertrim.sorted.bam | tee -a \$log

    samtools index ampliconclip/${prefix}.primertrim.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process samtools_markdup {
  tag        "${sample}"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/samtools:1.19'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(sample), val(type), file(sam) 

  output:
  tuple val(sample), file("markdup/${sample}.markdup.sorted.bam"), file("markdup/${sample}.markdup.sorted.bam.bai"), emit: bam_bai
  path "markdup/${sample}_markdupstats.txt", emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions
  
  shell:
  def args   = task.ext.args   ?: "${params.samtools_markdup_options}"
  def fmargs = task.ext.fmargs ?: "${params.samtools_fixmate_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  if ( type == 'single' ) {
  """
    mkdir -p markdup logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort ${sam} | \
      samtools markdup ${args} -@ ${task.cpus} -s -f markdup/${prefix}_markdupstats.txt - markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    samtools index markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
  } else if (type == 'paired') {
  """
    mkdir -p markdup logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -n ${sam} | \
      samtools fixmate ${fmargs} -m -@ ${task.cpus} - - | \
      samtools sort | \
      samtools markdup ${args} -@ ${task.cpus} -s -f markdup/${prefix}_markdupstats.txt - markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    samtools index markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
  }
}

