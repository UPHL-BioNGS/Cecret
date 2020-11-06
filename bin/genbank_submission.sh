USAGE="
\ngenbank_submission.sh is meant to combine consensus fastas
\nwith enough metadata to make submitting to GenBank easy. It
\ntakes the metadata found in covid_samples.txt and adds them
\nto the header line. The following headers are accepted:\n
\n
\tSample_id\t\t(required, must match sample_id*.fa*)\n
\tSubmission_id\t\t(if file needs renaming)\n
\tCountry\t\t\t(default is 'USA')\n
\tHost\t\t\t(default is 'Human')\n
\tIsolate\t\t\t(default is submission_id)\n
\tCollection_Date\n
\tIsolation_Source\t(default is 'SARS-CoV-2/host/location/isolate/date')\n
\tClone\n
\tCollected_By\n
\tFwd_Primer_Name\n
\tFwd_Primer_Seq\n
\tLatitude_Longitude\n
\tRev_Primer_Name\n
\tRev_Primer_Seq\n
\tNote\n
\tBioproject\n
\tBiosample\n
\tSra\n
\n
./genbank_submission.sh\n
\t-f <fasta with sample_id in first column of metadata file and in filename like sample_id*.fa*>\n
\t-m <file with metadata. Default is covid_samples.txt>\n
\t-o <output fasta for genbank submission>\n
\t-y <(optional) 4 digit year for isolate label. Default is current year ($(date +"%Y"))>\n
"
############################################################

if [ -z "$workdir" ] ; then workdir=$(pwd) ; fi

while getopts 'f:m:o:y:hv' OPTION
do
  case "$OPTION" in
    f)
    echo "The fasta is $OPTARG"
    fasta_file=$OPTARG
    ;;
    m)
    echo "The file with metadata for samples is $OPTARG"
    metadata_file=$OPTARG
    ;;
    o)
    echo "The final fasta is $OPTARG"
    genbank_fasta=$OPTARG
    ;;
    y)
    echo "The year is $OPTARG"
    year=$OPTARG
    ;;
    h)
    echo -e $USAGE
    exit 0
    ;;
    v)
    echo "Version 0.20201005"
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

if [ ! -f "$metadata_file" ]
then
  echo "File with metadata for genbank header not found at $metadata_file"
  echo -e $USAGE
  exit 1
fi

if [ ! -f "$fasta_file" ]
then
  echo "Fasta file not found at $fasta_file"
  echo -e $USAGE
  exit 1
fi
fasta=$(basename $fasta_file)

if [ -z "$genbank_fasta" ] && [ ! -d "$(pwd)/covid/submission_files" ]
then
  echo "-o specifies the name of the final file for genbank submission"
  echo -e $USAGE
  exit 1
fi

if [ -z "$year" ] ; then year=$(date +"%Y") ; fi

date

header_columns=($(head -n 1 $metadata_file ))
sample_id_column=$(echo ${header_columns[@]}          | tr " " "\n" | tr "\t" "\n" | grep -in "sample_id"                    | cut -f 1 -d ":" | head -n 1 )
if [ -z "$sample_id_column" ] ; then echo -e "No sample_id column. Could not determine ids for samples.\nFix $metadata_file so that one of the column headers is sample_id.\n$USAGE" ; exit 1 ; fi
submission_id_column=$(echo ${header_columns[@]}      | tr " " "\n" | tr "\t" "\n" | grep -in "submission"                   | cut -f 1 -d ":" | head -n 1 )
country_column=$(echo ${header_columns[@]}            | tr " " "\n" | tr "\t" "\n" | grep -in "country"                      | cut -f 1 -d ":" | head -n 1 )
host_column=$(echo ${header_columns[@]}               | tr " " "\n" | tr "\t" "\n" | grep -in "host"                         | cut -f 1 -d ":" | head -n 1 )
isolate_column=$(echo ${header_columns[@]}            | tr " " "\n" | tr "\t" "\n" | grep -in "isolate"                      | cut -f 1 -d ":" | head -n 1 )
collection_date_column=$(echo ${header_columns[@]}    | tr " " "\n" | tr "\t" "\n" | grep -in "collection" | grep -i "date"  | cut -f 1 -d ":" | head -n 1 )
isolation_source_column=$(echo ${header_columns[@]}   | tr " " "\n" | tr "\t" "\n" | grep -in "isolation" | grep -i "source" | cut -f 1 -d ":" | head -n 1 )
clone_column=$(echo ${header_columns[@]}              | tr " " "\n" | tr "\t" "\n" | grep -in "clone"                        | cut -f 1 -d ":" | head -n 1 )
collected_by_column=$(echo ${header_columns[@]}       | tr " " "\n" | tr "\t" "\n" | grep -in "collected" | grep -i "by"     | cut -f 1 -d ":" | head -n 1 )
fwd_primer_name_column=$(echo ${header_columns[@]}    | tr " " "\n" | tr "\t" "\n" | grep -in "Fwd-Primer-Name"              | cut -f 1 -d ":" | head -n 1 )
fwd_primer_seq_column=$(echo ${header_columns[@]}     | tr " " "\n" | tr "\t" "\n" | grep -in "Fwd-Primer-Seq"               | cut -f 1 -d ":" | head -n 1 )
latitude_longitude_column=$(echo ${header_columns[@]} | tr " " "\n" | tr "\t" "\n" | grep -in "Latitude-Longitude"           | cut -f 1 -d ":" | head -n 1 )
rev_primer_name_column=$(echo ${header_columns[@]}    | tr " " "\n" | tr "\t" "\n" | grep -in "Rev-Primer-Name"              | cut -f 1 -d ":" | head -n 1 )
rev_primer_seq_column=$(echo ${header_columns[@]}     | tr " " "\n" | tr "\t" "\n" | grep -in "Rev-Primer-Seq"               | cut -f 1 -d ":" | head -n 1 )
note_column=$(echo ${header_columns[@]}               | tr " " "\n" | tr "\t" "\n" | grep -in "Note"                         | cut -f 1 -d ":" | head -n 1 )
bioproject_column=$(echo ${header_columns[@]}         | tr " " "\n" | tr "\t" "\n" | grep -in "Bioproject"                   | cut -f 1 -d ":" | head -n 1 )
biosample_column=$(echo ${header_columns[@]}          | tr " " "\n" | tr "\t" "\n" | grep -in "Biosample"                    | cut -f 1 -d ":" | head -n 1 )
sra_column=$(echo ${header_columns[@]}                | tr " " "\n" | tr "\t" "\n" | grep -in "Sra"                          | cut -f 1 -d ":" | head -n 1 )

while read line
do
  line=$(echo $line | sed 's/ /\t/g')
  sample_id=$(echo $line | cut -f $sample_id_column -d " " | head -n 1)
  if [[ "$fasta" == "$sample_id"*.fa* ]]
  then
    if [ -n "$submission_id_column" ] ; then submission_id=$(echo $line | cut -f $submission_id_column -d " "  ) ; fi
    if [ -z "$submission_id" ] ; then submission_id=$sample_id ; fi
    fasta_header=">$submission_id "

    if [ -n "$country_column" ]
    then
      country=$(echo $line | cut -f $country_column -d " "  )
    else
      country="USA"
    fi
    fasta_header="$fasta_header[Country=$country]"

    if [ -n "$host_column" ]
    then
      host=$(echo $line | cut -f $host_column -d " "  )
    else
      host="Human"
    fi
    fasta_header="$fasta_header[Host=$host]"

    if [ -n "$isolate_column" ]
    then
      isolate=$(echo $line | cut -f $isolate_column -d " "  )
      fasta_header="$fasta_header[Isolate=$isolate]"
    else
      fasta_header="$fasta_header[Isolate=SARS-CoV-2/$host/$country/$submission_id/$year]"
    fi

    if [ -n "$collection_date_column" ]
    then
      collection_date=$(echo $line | cut -f $collection_date_column -d " "  )
      fasta_header="$fasta_header[Collection-Date=$collection_date]"
    fi

    if [ -n "$isolation_source_column" ]
    then
      isolation_source=$(echo $line | cut -f $isolation_source_column -d " "  )
      fasta_header="$fasta_header[Isolation-Source=$isolation_source]"
    fi

    if [ -n "$clone_column" ]
    then
      clone=$(echo $line | cut -f $clone_column -d " "  )
      fasta_header="$fasta_header[Clone=$clone]"
    fi

    if [ -n "$collected_by_column" ]
    then
      collected_by=$(echo $line | cut -f $collected_by_column -d " "  )
      fasta_header="$fasta_header[Collected-By=$collected_by]"
    fi

    if [ -n "$fwd_primer_name_column" ]
    then
      fwd_name=$(echo $line | cut -f $fwd_primer_name_column -d " "  )
      fasta_header="$fasta_header[Fwd-Primer-Name=$fwd_name]"
    fi

    if [ -n "$fwd_primer_seq_column" ]
    then
      fwd_seq=$(echo $line | cut -f $fwd_primer_seq_column -d " "  )
      fasta_header="$fasta_header[Fwd-Primer-Seq=$fwd_seq]"
    fi

    if [ -n "$latitude_longitude_column" ]
    then
      lat_lon=$(echo $line | cut -f $latitude_longitude_column -d " "  )
      fasta_header="$fasta_header[Latitude-Longitude=$lat_lon]"
    fi

    if [ -n "$rev_primer_name_column" ]
    then
      rev_name=$(echo $line | cut -f $rev_primer_name_column -d " "  )
      fasta_header="$fasta_header[Rev-Primer-Name=$rev_name]"
    fi

    if [ -n "$rev_primer_seq_column" ]
    then
      rev_seq=$(echo $line | cut -f $rev_primer_seq_column -d " "  )
      fasta_header="$fasta_header[Rev-Primer-Seq=$rev_seq]"
    fi

    if [ -n "$note_column" ]
    then
      note=$(echo $line | cut -f $note_column -d " "  )
      fasta_header="$fasta_header[Note=$note]"
    fi

    if [ -n "$bioproject_column" ]
    then
      bioproject=$(echo $line | cut -f $bioproject_column -d " "  )
      fasta_header="$fasta_header[Bioproject=$bioproject]"
    fi

    if [ -n "$biosample_column" ]
    then
      biosample=$(echo $line | cut -f $biosample_column -d " "  )
      fasta_header="$fasta_header[Biosample=$biosample]"
    fi

    if [ -n "$sra_column" ]
    then
      sra=$(echo $line | cut -f $sra_column -d " "  )
      fasta_header="$fasta_header[Sra=$sra]"
    fi

    if [ -z "$genbank_fasta" ]
    then
      mkdir -p $(pwd)/covid/submission_files
      final_fasta="$(pwd)/covid/submission_files/$submission_id.genbank.fa"
    else
      final_fasta=$genbank_fasta
    fi

    echo $fasta_header > $final_fasta
    grep -v ">" $fasta_file | sed 's/^N*N//g' | fold -w 75 >> $final_fasta
  fi
done < $metadata_file

if [ -f "$final_fasta" ]
then
  echo "$fasta_file has been converted to genbank submission file"
  echo "The file for easy genbank submission is located at"
  echo $final_fasta
else
  echo "File could not be made for $fasta_file"
  echo "Ensure sample is found in $metadata_file"
fi

exit 0
