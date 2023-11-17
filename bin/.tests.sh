#/bin/bash
# nextflow run ~/sandbox/Cecret -profile singularity --reads /home/eriny/sandbox/test_files/cecret/reads --outdir tests -with-tower -resume
# nextflow run ~/sandbox/Cecret -profile singularity,mpx --reads /home/eriny/sandbox/test_files/cecret/mpx --outdir tests --cleaner 'fastp' -with-tower -resume
# nextflow run ~/sandbox/Cecret -profile singularity --nanopore /home/eriny/sandbox/test_files/cecret/nanopore --outdir tests -with-tower -resume
# nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret -profile singularity --reads /home/eriny/sandbox/test_files/cecret/reads --outdir tests -with-tower -resume

echo "usage : bash .test.sh {small,primers,else}"

test=$1

if [ -z "$test" ]; then test="small" ; fi

if [ "$test" == "small" ]
then

  # sample sheet
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/cecret/sample_sheet.csv \
    --outdir singularity_defaults_sample_sheet \
    -resume \
    -with-tower

  # using included nextclade data
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/cecret/sample_sheet.csv \
    --outdir singularity_defaults_sample_sheet_included_nextclade \
    --download_nextclade_dataset false \
    -resume \
    -with-tower

  options=("reads" "single_reads" "fastas" "nanopore")

  for option in ${options[@]}
  do
    # defaults
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir singularity_defaults_$option \
      -resume \
      -with-tower

    # removed test for bamsnap and rename because of lack of interest
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir all_on_$option \
      -with-tower \
      --filter true \
      -resume

    # removing primer trimming
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir nontrimmed_$option \
      --trimmer 'none' \
      -with-tower \
      -resume

    # changing the cleaner, aligner, and trimmer
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir toggled_$option \
      --cleaner 'fastp' \
      --trimmer 'samtools' \
      --aligner 'minimap2' \
      --markdup true \
      -with-tower \
      -resume

    # with UPHL's config
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile uphl,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir uphl_$option \
      -with-tower \
      -resume
  done

  # multifasta
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/test_files/cecret/reads \
    --single-reads /home/eriny/sandbox/test_files/cecret/single-reads \
    --fastas /home/eriny/sandbox/test_files/cecret/fastas \
    --multifastas /home/eriny/sandbox/test_files/cecret/multifasta \
    --outdir kitchen_sink  \
    -with-tower

  # empty
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,artic_V3 \
    --reads doesntexit \
    --single-reads willnotexist \
    --fastas shouldntexit \
    --outdir empty  \
    -with-tower

elif [ "$test" == "primers" ]
then

  nanopore_primers=("midnight_idt_V1" "midnight_ont_V1" "midnight_ont_V2" "midnight_ont_V3")
  for primer in ${nanopore_primers[@]}
  do
    nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
      -profile singularity \
      --nanopore /home/eriny/sandbox/test_files/cecret/nanopore \
      --outdir primer_${primer}_artic \
      --primer_set $primer \
      --bcftools_variants false \
      --fastqc false \
      --ivar_variants false \
      --samtools_stats false \
      --samtools_coverage false \
      --samtools_depth false \
      --nextclade false \
      --pangolin false \
      --freyja false \
      --vadr false \
      -resume \
      -with-tower
  done

  trimmers=("samtools" "ivar")
  for trimmer in ${trimmers[@]}
  do
    illumina_primers=("ncov_V3" "ncov_V4" "ncov_V4.1" "ncov_V5.3.2")
    for primer in ${illumina_primers[@]}
    do
      nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
        -profile singularity \
        --reads /home/eriny/sandbox/test_files/cecret/reads \
        --outdir primer_${primer}_${trimmer} \
        --trimmer $trimmer \
        --cleaner 'fastp' \
        --primer_set $primer \
        --bcftools_variants false \
        --fastqc false \
        --ivar_variants false \
        --samtools_stats false \
        --samtools_coverage false \
        --samtools_depth false \
        --nextclade false \
        --pangolin false \
        --freyja false \
        --vadr false \
        -resume \
        -with-tower
    done

    mpx_primers=("mpx_primalseq" "mpx_idt")
    for primer in ${mpx_primers[@]}
    do
      nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
        -profile singularity \
        --reads /home/eriny/sandbox/test_files/cecret/mpx_idt \
        --outdir primer_${primer}_${trimmer} \
        --trimmer $trimmer \
        --cleaner 'fastp' \
        --primer_set $primer \
        --bcftools_variants false \
        --fastqc false \
        --ivar_variants false \
        --samtools_stats false \
        --samtools_coverage false \
        --samtools_depth false \
        --nextclade false \
        --pangolin false \
        --freyja false \
        --vadr false \
        -with-tower \
        -resume
    done
  done

else
  # CDC's test data with relatedness
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir default_datasets \
    --relatedness true  \
    -with-tower

  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile uphl,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir uphl_datasets \
    -with-tower \
    -resume \
    --relatedness true

  # CDC's test data with relatedness using nextalign
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir toggled_datasets \
    --cleaner 'fastp' \
    --trimmer 'samtools' \
    --aligner 'minimap2' \
    --markdup true \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume

  # MPX
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,mpx \
    --reads /home/eriny/sandbox/test_files/cecret/mpx \
    --outdir mpx \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume

  # MPX with idt primers
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity,mpx_idt \
    --reads /home/eriny/sandbox/test_files/cecret/mpx_idt \
    --outdir mpx_idt \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume

  # other
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/cecret/mpx \
    --outdir mpx \
    --cleaner 'fastp' \
    --trimmer 'none' \
    --species 'other' \
    --nextclade_dataset  'hMPXV' \
    --vadr_options '--split --glsearch -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin' \
    --vadr_reference 'mpxv' \
    --vadr_trim_options '--minlen 50 --maxlen 210000' \
    --kraken2_organism 'Monkeypox virus' \
    -with-tower \
    -resume

  # MPX with kraken
  nextflow run /Volumes/IDGenomics_NAS/Bioinformatics/eriny/Cecret \
    -profile uphl,mpx \
    --reads /home/eriny/sandbox/test_files/cecret/mpx \
    --outdir uphl_mpx \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume
fi
