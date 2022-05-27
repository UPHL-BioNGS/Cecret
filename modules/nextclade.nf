process nextclade {
  tag "Clade Determination"
  label "medcpus"

  when:
  params.nextclade || params.msa == 'nextalign'

  input:
  file(fasta)

  output:
  path "nextclade/nextclade.csv",                                              emit: nextclade_file
  path "nextclade/*",                                                          emit: results
  path "nextclade/nextclade.aligned.fasta",                                    emit: nextclade_aligned_fasta
  path "dataset",                                                              emit: prepped_nextalign
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p nextclade dataset logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file
    nextclade_version=$(nextclade --version)

    nextclade dataset get --name !{params.nextclade_dataset} --output-dir dataset

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    nextclade !{params.nextclade_options} \
      --input-fasta=ultimate_fasta.fasta \
      --input-dataset dataset \
      --output-json=nextclade/nextclade.json \
      --output-csv=nextclade/nextclade.csv \
      --output-tsv=nextclade/nextclade.tsv \
      --output-tree=nextclade/nextclade.auspice.json \
      --output-dir=nextclade \
      --output-basename=nextclade \
      2>> $err_file >> $log_file
    cp ultimate_fasta.fasta nextclade/combined.fasta
  '''
}
