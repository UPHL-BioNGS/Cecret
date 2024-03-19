process multiqc_combine {
  tag        "multiqc"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/multiqc:1.19'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.multiqc && (task.ext.when == null || task.ext.when)

  input:
  file(input)
  file(script)

  output:
  path "multiqc/multiqc_report.html", optional: true, emit: html
  path "multiqc/multiqc_data/*", optional: true, emit: files
  path "multiqc/multiqc_data", optional: true, emit: multiqc_data
  path "software_versions.yml", optional: true, emit: versions
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  def args = task.ext.args ?: "${params.multiqc_options}"
  """
    mkdir -p multiqc logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    multiqc --version >> \$log

    touch whatever.png
    pngs=\$(ls *png | grep -v mqc.png\$ )

    for png in \${pngs[@]}
    do
      file=\$(echo \$png | sed 's/.png\$//g' )
      mv \${file}.png \${file}_mqc.png
    done

    rm whatever_mqc.png

    cp collated_versions.yml versions.yml

    cat <<-END_VERSIONS >> versions.yml
    "CECRET":
      workflow: ${workflow.manifest.version}
    "${task.process}":
      multiqc: \$( multiqc --version | awk '{print \$NF}' )
      container: ${task.container}
    END_VERSIONS

    python3 ${script}

    multiqc ${args} \
      --outdir multiqc \
      . \
      | tee -a \$log
  """
}
