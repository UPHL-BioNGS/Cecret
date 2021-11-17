#/bin/bash

test=$1

if [ -z "$test" ]; then test="small" ; fi

if [ "$test" == "small" ]
then
  options=("reads" "single_reads" "fastas")

  for option in ${options[@]}
  do
    # defaults with docker
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev \
      -profile docker \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir docker_defaults_$option \
      -with-tower

    # defaults
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev \
      -profile singularity \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir singularity_defaults_$option \
      -with-tower

    # removed test for bamsnap and rename because of lack of interest
    # attempted bcftools and filter
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev \
      -profile singularity \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir all_on_$option \
      -with-tower \
      --bcftools_variants true \
      --filter true \
      -resume

    # changing the cleaner, aligner, and trimmer
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev \
      -profile singularity \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir toggled_$option \
      --cleaner 'fastp' \
      --trimmer 'samtools' \
      --aligner 'minimap2' \
      -with-tower \
      -resume

    # removing primer trimming
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev \
      -profile singularity \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir nontrimmed_$option \
      --trimmer 'none' \
      -with-tower \
      -resume

    # with UPHL's config
    nextflow run UPHL-BioNGS/Cecret \
      -r erin-dev -profile uphl \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir uphl_$option \
      -with-tower \
      -resume
  done

  #################
  ##### empty #####
  #################
  nextflow run UPHL-BioNGS/Cecret \
    -r erin-dev \
    -profile singularity \
    --reads doesntexit \
    --single-reads willnotexist \
    --fastas shouldntexit \
    --outdir empty  \
    -with-tower
else
  # CDC's test data with relatedness
  nextflow run UPHL-BioNGS/Cecret \
    -r erin-dev \
    -profile singularity \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads/ \
    --outdir default_datasets \
    --relatedness true  \
    -with-tower

  nextflow run UPHL-BioNGS/Cecret \
    -r erin-dev -profile uphl \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads/ \
    --outdir uphl_datasets \
    -with-tower \
    -resume \
    --relatedness true

  # CDC's test data with relatedness using nextalign
  nextflow run UPHL-BioNGS/Cecret \
    -r erin-dev \
    -profile singularity \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads/ \
    --outdir toggled_datasets \
    --cleaner 'fastp' \
    --trimmer 'samtools' \
    --aligner 'minimap2' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume
fi
