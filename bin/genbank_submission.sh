USAGE="
genbank_submission.sh is meant to combine consensus fastas
with enough metadata to make submitting to GenBank easy. It
takes the metadata found in covid_samples.txt and adds them
to the header line. Accession and submission ids as well as
collection dates are required.

./genbank_submission \\
  -d <directory where covid/consensus.fasta and covid_samples.txt are> \\
  -f <file with metadata. Default is covid_samples.txt> \\
  -y <4 digit year. Default is current year ($(date +"%Y"))> \\
  -s <(optional) specific sample. Default is to create file for all listed in covid_samples.txt>
"
############################################################

if [ -z "$workdir" ] ; then workdir=$(pwd) ; fi

while getopts 'd:y:s:f:hv' OPTION
do
  case "$OPTION" in
    d)
    echo "The directory is $OPTARG"
    workdir=$OPTARG
    ;;
    f)
    echo "The file with metadata for samples is $OPTARG"
    metadata_file=$OPTARG
    ;;
    y)
    echo "The year is $OPTARG"
    year=$OPTARG
    ;;
    s)
    echo "The sample is $OPTARG"
    samples=($OPTARG)
    ;;
    h)
    echo "$USAGE"
    exit 0
    ;;
    v)
    echo "Version 0.20200902"
    exit 0
    ;;
    :)
    echo "Invalid option: $OPTARG requires an argument"
    echo "$USAGE"
    exit 1
    ;;
    \?)
    echo "$USAGE"
    exit 1
    ;;
  esac
done
shift "$(($OPTIND -1))"

if [ ! -f "$metadata_file" ]
then
  echo "File with metadata for genbank header not found at $metadata_file"
  echo $USAGE
  exit 1
fi

if [ -z "$year" ] ; then year=$(date +"%Y") ; fi
if [ -z "${samples[0]}" ] ; then samples=($(tail -n+2 $metadata_file | awk '{print $1}' | sort | uniq )) ; fi
metadata_samples=($(tail -n+2 $metadata_file | awk '{print $1}' | sort | uniq ))

submission_file_directory=$workdir/covid/submission_files
mkdir -p "$submission_file_directory"
if [ -z "$submission_file_directory" ]
then
  echo "Could not create final directory for files at $submission_file_directory"
  echo $USAGE
  exit 1
fi

header_columns=($(head -n 1 $metadata_file ))
submission_id_columns=($(history -p ${header_columns[@]} | grep -in "submission" | head -n 1 | cut -f 1 -d ":" ))
if [ -z "$submission_id_columns" ] ; then submission_id_columns=1 ; fi
array_location=$(( $submission_id_columns - 1 ))

for sample in ${samples[@]}
do
  (
  for metadata_sample in ${metadata_samples[@]}
  do
    if [[ "$sample" == "$metadata_sample"* ]] ; then sample=$metadata_sample ; fi
  done
  line=($(grep $sample $metadata_file | head -n 1 ))
  submission_id=${line[$array_location]}
  if [ -z "$submission_id" ] ; then submission_id=$sample ; fi

  fasta_header=$(echo ">$submission_id ")
  for (( i=0 ; i<${#header_columns[@]}; i++ ))
  do
    header_test=$(echo ${header_columns[$i]} | grep -i -e "accession" -e "EpiTrax" -e "confidential" -e "zip" -e "id" )
    if [ -z "$header_test" ]
    then
      fasta_header="$fasta_header[${header_columns[$i]}=${line[$i]}]"
    fi
  done
  organism_test=$(echo $fasta_header | grep -i "organism")
  if [ -z "$organism_test" ] ; then fasta_header="$fasta_header[organism=Severe acute respiratory syndrome coronavirus 2]" ; fi
  host_test=$(echo $fasta_header | grep -i "host")
  if [ -z "$host_test" ] ; then fasta_header="$fasta_header[host=human]" ; fi
  country_text=$(echo $fasta_header | grep -i "country")
  if [ -z "$country_text" ] ; then fasta_header="$fasta_header[country=USA]" ; fi
  isolate_test=$(echo $fasta_header | grep -i "isolate")
  if [ -z "$isolate_test" ]
  then
    host=$(echo $fasta_header | sed 's/.*host=//g' | sed 's/\].*//g' )
    country=$(echo $fasta_header | sed 's/.*country=//g' | sed 's/\].*//g' )
    if [ -z "$year" ] ; then year=$(echo $fasta_header | sed 's/.*year=//g' | sed 's/\].*//g' ) ; fi
    fasta_header="$fasta_header[isolate=SARS-CoV-2/$host/$country/$submission_id/$year]"
  fi
  echo $fasta_header > $submission_file_directory/$submission_id.genbank.fa
  if [ -d "$workdir/covid/consensus" ]
  then
    consensus_file=$(ls $workdir/covid/consensus/$sample*.consensus.fa)
  else
    consensus_file=$(ls $sample*.consensus.fa)
  fi
  grep -v ">" $consensus_file | sed 's/^N*N//g' | fold -w 75 >> $submission_file_directory/$submission_id.genbank.fa
  ) &
  if [ $( jobs | wc -l ) -ge $( nproc ) ] ; then wait ; fi
done

exit 0
