#/bin/bash
# ~/sandbox/UPHL/COVID/files_for_submission.sh $(pwd)

out=$1

covid_samples=$out/covid_samples.txt
covid_summary=$out/covid/summary.txt
run_results=$out/run_results.txt

if [ -n "$out" ] && [ -f "$covid_samples" ] && [ -f "$covid_summary" ] && [ -f "$run_results" ]
then
  mkdir -p $out/covid/submission_files/
else
  echo "Usage: files_for_submission.sh /home/IDGenomics_NAS/WGS_Serotyping/UT-M70330-200313"
  echo "Needs $covid_samples with tab delimiated columns of lab_accession  submission_id collection_date"
  echo "Needs $covid_summary from Illumina_covid_V3.nf pipeline"
  echo "Needs $run_results from UPHL_reference_free pipeline + URF_scripts/run_results_urf.sh"
  exit
fi

rm -f $out/covid/submission_files/*_submission.fasta

while read line
do
  sample_id=$(echo $line       | awk '{print $1}' )
  submission_id=$(echo $line   | awk '{print $2}' )
  collection_date=$(echo $line | awk '{print $3}' )
  sample_run_id=$(grep $sample_id $out/run_results.txt | awk '{print $2}' )
  run_id=${sample_run_id: -6}
  num_n=$(grep $sample_id $covid_summary | cut -f 7 -d ',')
  if [ -z "$num_n" ] ; then num_n=0 ; fi

  # getting the consensus fasta file
  # changing the fasta header
  cat $out/covid/consensus/$sample_run_id.consensus.fa | sed "s/>.*/>$submission_id/g" > $out/covid/submission_files/$submission_id.consensus.fa

### old and commented out code that may need to be brought back
##  # replacing all degenerate bases with N, removing leading Ns, and folding sequence to 75 bp wide
#  grep -v ">" $out/covid/submission_files/$submission_id.consensus.fa | sed 's/[BDEFHIJKLMOPQRSUVWXYZ]/N/g' | sed 's/^N*N//g' | fold -w 75 >> $out/covid/submission_files/$submission_id.nondegenerate.fa

  # removing leading Ns, and folding sequencing to 75 bp wide
  echo ">$submission_id" > $out/covid/submission_files/$submission_id.nondegenerate.fa
  grep -v ">" $out/covid/submission_files/$submission_id.consensus.fa | sed 's/^N*N//g' | fold -w 75 >> $out/covid/submission_files/$submission_id.nondegenerate.fa
  cat $out/covid/submission_files/$submission_id.nondegenerate.fa | awk -v id=$submission_id -v cldt=$collection_date '{ if ($0 ~ ">") print $0 " [organism=Severe acute respiratory syndrome coronavirus 2][isolate=SARS-CoV-2/Human/USA/" id "/2020][host=Human][country=USA][collection_date=" cldt "]" ; else print $0 }' > $out/covid/submission_files/$submission_id.genbank.fa

  # copying fastq files and changing the file name
  cp $out/Sequencing_reads/Raw/$sample_run_id*R1_001.fastq.gz $out/covid/submission_files/$submission_id.R1.fastq.gz
  cp $out/Sequencing_reads/Raw/$sample_run_id*R2_001.fastq.gz $out/covid/submission_files/$submission_id.R2.fastq.gz

  # preparing fasta for gisaid and genbank submission
  if [ "$num_n" -lt 14952 ]
  then
    cat $out/covid/submission_files/$submission_id.genbank.fa >> $out/covid/submission_files/$run_id.genbank_submission.fasta
    if [ "$num_n" -lt 4903 ]
    then
      cat $out/covid/submission_files/$submission_id.consensus.fa | sed "s/^>/>hCoV-19\/USA\//g" | awk '{ if ($0 ~ ">") print $0 "/2020" ; else print $0}' >> $out/covid/submission_files/$run_id.gisaid_submission.fasta
    fi
  fi
done < <( grep -v "Submission_ID" $covid_samples )
