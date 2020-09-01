workdir=$1
samples=($2)

if [ -z "$workdir" ] ; then workdir=$(pwd) ; fi
if [ ! -f "$workdir/covid_samples.txt" ] ; then echo "File with metadata for genbank header not found at $workdir/covid_samples.txt" ; exit ; fi
mkdir -p $workdir/covid/submission_files

if [ -z "$sample" ]; then samples=($(tail -n+2 $workdir/covid_samples.txt | awk '{print $1}' | sort | uniq )) ; fi
header_columns=($(head -n 1 $workdir/covid_samples.txt ))
submission_id_columns=($(history -p ${header_columns[@]} | grep -in "submission" | head -n 1 | cut -f 1 -d ":" ))
if [ -z "$submission_id_columns" ] ; then submission_id_columns=1 ; fi

for sample in ${samples[@]}
do
  line=($(grep $sample $workdir/covid_samples.txt | head -n 1 ))
  submission_id=${line[$(($submission_id_columns - 1))]}

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
  year_test=$(echo $fasta_header | grep -i "year" )
  if [ -z "$year_test" ] ; then year=$(date +"%Y") ; fi
  isolate_test=$(echo $fasta_header | grep -i "isolate")
  if [ -z "$organism_test" ]
  then
    host=$(echo $fasta_header | sed 's/.*host=//g' | sed 's/\].*//g' )
    country=$(echo $fasta_header | sed 's/.*country=//g' | sed 's/\].*//g' )
    if [ -z "$year" ] ; then year=$(echo $fasta_header | sed 's/.*year=//g' | sed 's/\].*//g' ) ; fi
    fasta_header="$fasta_header[isolate=SARS-CoV-2/$host/$country/$submission_id/$year]"
  fi
  echo $fasta_header > $workdir/covid/submission_files/$submission_id.genbank.fa
  grep -v ">" $workdir/covid/consensus/$sample*fa* | sed 's/^N*N//g' | fold -w 75 >> $workdir/covid/submission_files/$submission_id.genbank.fa
done
