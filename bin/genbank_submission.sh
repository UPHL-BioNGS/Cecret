#!/bin/bash

version="0.20210430"
USAGE="
\ngenbank_submission.sh is meant to combine consensus fastas
\nwith enough metadata to make submitting to GenBank easy. It
\ntakes the metadata found in covid_samples.txt and adds them
\nto the header line. The following headers are required:\n
\n
\tSample_id\t\t(required, must match sample_id*.fa)\n
\tSubmission_id\t\t(if file needs renaming)\n
\tCollection_Date\n
\nMany other headers are accepted.
\n
./genbank_submission.sh\n
\t-f The file with metadata. Default is covid_samples.csv\n
\t-c Directory with consensus fastas.\n
\t-d Directory with fastq files (will be gzip to fastq.gz).\n
\t-s The threshold of non-ambiguous basees for GISAID submission\n
\t-g The threshold of non-ambiguous basees for GENBANK submission\n
\t-o The directory where the final files will be.\n
\t-h print this help message\n
\t-v print version and exist\n

Usage:
genbank_submission.sh -f covid_samples.csv -c . -d . -s 25000 -g 15000 -o submission_files
"
############################################################

sample_file=covid_samples.csv
consensus_directory=cecret/consensus
fastq_directory=Sequencing_reads/Raw
gisaid_threshold=25000
genbank_threshold=15000
out=submission_files

while getopts 'f:c:d:s:g:o:hv' OPTION
do
  case "$OPTION" in
    f)
    echo "The file with sample metadata is $OPTARG"
    sample_file=$OPTARG
    ;;
    c)
    echo "The directory with the consensus fasta files is $OPTARG"
    consensus_directory=$OPTARG
    ;;
    d)
    echo "The directory with the fastq files is $OPTARG"
    fastq_directory=$OPTARG
    ;;
    s)
    echo "The number of non-ambiguous bases for GISAID submission is $OPTARG"
    gisaid_threshold=$OPTARG
    ;;
    g)
    echo "The number of non-ambiguous bases for GENBANK submission is $OPTARG"
    genbank_threshold=$OPTARG
    ;;
    o)
    echo "The output directory for renamed files $OPTARG"
    out=$OPTARG
    ;;
    h)
    echo -e $USAGE
    exit 0
    ;;
    v)
    echo "Version $version"
    exit 0
    ;;
    :)
    echo "Invalid option: $OPTARG requires an argument"
    echo -e $USAGE
    exit 1
    ;;
    \?)
    echo -e $USAGE
    exit 1
    ;;
  esac
done
shift "$(($OPTIND -1))"

if [ ! -f "$sample_file" ]
then
  echo "File with metadata for genbank header not found at $metadata_file"
  echo -e $USAGE
  exit 1
fi

if [ ! -d "$consensus_directory" ]
then
  echo "The directory with fasta files not found at $consensus_directory"
  echo -e $USAGE
  exit 1
fi

if [ ! -d "$fastq_directory" ]
then
  echo "The directory with fastq files not found at $consensus_directory"
  echo -e $USAGE
  exit 1
fi

mkdir -p $out
if [ ! -d "$out" ]
then
  echo "Could not create directory for final files at $out"
  echo -e $USAGE
  exit 1
fi

sample_id_column=$(head -n 1 $sample_file | tr ',' '\n' | grep -inw "Sample_ID" | cut -f 1 -d ':' )
submission_id_column=$(head -n 1 $sample_file | tr ',' '\n' | grep -inw "Submission_id" | cut -f 1 -d ':' )
collection_date_column=$(head -n 1 $sample_file | tr ',' '\n' | grep -inw "Collection_Date" | cut -f 1 -d ':'  )

if [ -z $sample_id_column ] || [ -z "$submission_id_column" ] || [ -z "$collection_date_column" ]
then
  echo "$sample_file is not the correct format"
  echo "Sorry to be overly picky, but this file needs to be a plain text file with values separated by commas (and no commas in the values)"
  echo "Required headers are 'Sample_ID','Submission_ID','Collection_Date'"
  echo "Please read documentation at https://github.com/UPHL-BioNGS/Cecret"
  exit 1
fi
sample_file_header_reduced=$(head -n 1 $sample_file | tr "," '\n' | grep -iv "Sample_ID" | grep -iv "Collection_Date" | grep -vi "Submission_ID" | tr '\n' ' ' )

gzip $fastq_directory/*.fastq*

while read line
do
  sample_id=$(echo $line | cut -f $sample_id_column -d "," )
  submission_id=$(echo $line | cut -f $submission_id_column -d ",")
  collection_date=$(echo $line | cut -f $collection_date_column -d ",")
  if [ -z "$sample_id" ]
  then
    echo "FATAL : Could not find sample ids for all samples."
    echo "Please read documentation at https://github.com/UPHL-BioNGS/Cecret"
    exit 1
  fi
  if [ -z "$submission_id" ] ; then submission_id=$sample_id ; fi
  if [ -z "$collection_date" ] ; then collection_date="missing" ; fi
  genbank_fasta_header=">$submission_id "
  for column in ${sample_file_header_reduced[@]}
  do
    column_number=$(head -n 1 $sample_file | tr "," "\n" | grep -n "$column" | cut -f 1 -d ':')
    column_value=$(echo $line | cut -f $column_number -d ',')
    if [ -z "$column_value" ] ; then column_value="missing" ; fi
    genbank_fasta_header=$genbank_fasta_header"["$column"="$column_value"]"
  done

  if [ "$collection_date" == "missing" ]
  then
    year=$(date "+%Y")
    echo "The collection date is $collection_date for $sample_id" | tee -a $log_file
  else
    collection_date=$(date -d "$collection_date" "+%Y-%m-%d") || echo "Invalid date format. Try something like yyyy-mm-dd and '-resume' the workflow."
    year=$(date -d "$collection_date" "+%Y")
    echo "The collection date is $collection_date for $sample_id" | tee -a $log_file
    genbank_fasta_header=$genbank_fasta_header"[Collection_Date="$collection_date"]"
  fi

  country_check=$(echo $sample_file_header_reduced | grep -wi "country" | head -n 1 )
  if [ -z "$country_check" ]
  then
    genbank_fasta_header=$genbank_fasta_header"[Country=USA]"
    country="USA"
  else
    column_number=$(head -n 1 $sample_file | tr "," "\n" | grep -in "country" | cut -f 1 -d ':')
    country=$(echo $line | cut -f $column_number -d ',')
    if [ -z "$country" ] ; then country="missing" ; fi
  fi

  host_check=$(echo $sample_file_header_reduced | grep -wi "host"  | head -n 1 )
  if [ -z "$host_check" ]
  then
    genbank_fasta_header=$genbank_fasta_header"[Host=Human]"
    host="Human"
  else
    column_number=$(head -n 1 $sample_file | tr "," "\n" | grep -in "host" | cut -f 1 -d ':')
    host=$(echo $line | cut -f $column_number -d ',')
    if [ -z "$host" ] ; then host="missing" ; fi
  fi

  isolate_check=$(echo $sample_file_header_reduced | grep -wi "isolate"  | head -n 1 )
  if [ -z "$isolate_check" ]
  then
    organism_check=$(head -n 1 $sample_file | tr ',' '\n' | grep -i "organism" | head -n 1 )
    if [ -z "$organism_check" ]
    then
      genbank_organism='SARS-CoV-2'
      gisaid_organism='hCoV-19'
    else
      column_number=$(head -n 1 $sample_file | tr "," "\n" | grep -in "organism" | cut -f 1 -d ':')
      genbank_organism=$(echo $line | cut -f $column_number -d ',')
      gisaid_organism=$(echo $line  | cut -f $column_number -d ',')
      if [ -z "$genbank_organism" ] ; then genbank_organism="missing" ; fi
      if [ -z "$gisaid_organism" ] ; then gisaid_organism="missing" ; fi
    fi
    genbank_fasta_header=$genbank_fasta_header"[Isolate="$genbank_organism"/"$host"/"$country"/"$submission_id"/"$year"]"
    gisaid_fasta_header=">$gisaid_organism/$country/$submission_id/$year"
  fi

  consensus=($(ls $consensus_directory/*.fa | grep "/$sample_id"))
  echo "Found ${#consensus[@]} consensus fastas."
  if [ "${#consensus[@]}" == 1 ]
  then
    num_ACTG=$(grep -v ">" ${consensus[0]} | grep -o -E "C|A|T|G" | wc -l )
    if [ "$num_ACTG" -gt "$gisaid_threshold" ]
    then
      echo $gisaid_fasta_header > $out/$submission_id.gisaid.fa
      grep -v ">" ${consensus[0]} | fold -w 75 >> $out/$submission_id.gisaid.fa
    fi

    if [ "$num_ACTG" -gt "$genbank_threshold" ]
    then
      echo $genbank_fasta_header > $out/$submission_id.genbank.fa
      grep -v ">" ${consensus[0]} | sed 's/^N*N//g' | fold -w 75 >> $out/$submission_id.genbank.fa
    fi
  else
    echo "FATAL : sample_id wasn't specific enough and too many files match."
    echo "Consensus fasta files : ${consensus[@]}"
  fi

  fastq_files=($(ls $fastq_directory/*.fastq.gz | grep -v "filter" | grep -v "unpaired" | grep "/$sample_id" ))
  filtered_fastq_files=($(ls $fastq_directory/*.fastq.gz | grep "filter" | grep -v "unpaired" | grep "/$sample_id" ))
  echo "Found ${#fastq_files[@]} fastq files and ${#filtered_fastq_files[@]} filtered fastq files"
  if [ "${#filtered_fastq_files[@]}" == 2 ]
  then
    echo "Found filtered fastq"
    cp ${filtered_fastq_files[0]} $out/${submission_id}_filtered_R1.fastq.gz
    cp ${filtered_fastq_files[1]} $out/${submission_id}_filtered_R2.fastq.gz
  elif [ "${#filtered_fastq_files[@]}" == 1 ]
  then
    cp ${filtered_fastq_files[0]} $out/${submission_id}_filtered.fastq.gz
  elif [ "${#fastq_files[@]}" == 2 ]
  then
    cp ${fastq_files[0]} $out/${submission_id}_R1.fastq.gz
    cp ${fastq_files[1]} $out/${submission_id}_R2.fastq.gz
  elif [ "${#fastq_files[@]}" == 1 ]
  then
    cp ${fastq_files[0]} $out/${submission_id}.fastq.gz
  elif [ "${#fastq_files[@]}" == 0 ] && [ "${#fastq_unpaired_files[@]}" == 0 ]
  then
    echo "WARN : No fastq files found for $sample_id"
  else
    echo "WARN : Could not determine which fastqs were for $sample_id"
    echo "Fastq files : ${fastq_files[@]} "
  fi
done < <(grep -v "Collection_Date" $sample_file )

#cat $out/*genbank.fa > $out/genbank_submission.fasta
#cat $out/*gisaid.fa > $out/gisaid_submission.fasta
