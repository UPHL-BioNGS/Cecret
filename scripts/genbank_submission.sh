workdir=$1
samples=$2

if [ ! -f "$workdir/covid_samples.txt" ] ; then echo "File with metadata for genbank header not found" ; exit ; fi
if [ -z "$sample" ]
then
  samples=($(tail -n+2 $workdir/covid_samples.txt | awk '{print $1}' ))
fi
echo ${samples[@]}

header_columns=($(head -n 1 $workdir/covid_samples.txt ))
echo ${#header_columns[@]}

for (( i=0 ; i<${#header_columns[@]}; i++ ))
do
  echo $i
  echo ${header_columns[$i]}
done

for sample in ${samples[@]}
do
  fasta_header=$(echo ">$sample ")
  for (( i=0 ; i<${#header_columns[@]}; i++ ))
  do
    echo $i
    echo "addition check $((i + 1 ))"
    echo ${header_columns[$i]}
    grep $sample $workdir/covid_samples.txt | cut -f $((i + 1 ))
    value=$(grep $sample $workdir/covid_samples.txt | cut -f $((i + 1 )) | head -n 1)
    fasta_header="$fasta_header[${header_columns[$i]}=$value]"
    echo $fasta_header
    i=$((i + 1 ))
  done
  echo $fasta_header
done

exit
if [ -z $sample ]
then

  $workdir/covid_samples.txt


fi



echo ">!{submission_id} [organism=Severe acute respiratory syndrome coronavirus 2][isolate=SARS-CoV-2/Human/USA/!{submission_id}/!{params.year}][host=Human][country=USA][collection_date=!{collection_date}]" > covid/submission_files/!{submission_id}.genbank.fa  2>> $err_file
grep -v ">" !{consensus} | sed 's/^N*N//g' | fold -w 75 >> covid/submission_files/!{submission_id}.genbank.fa  2>> $err_file
