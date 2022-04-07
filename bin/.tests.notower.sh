#/bin/bash
#nextflow ~/sandbox/Cecret/Cecret.nf -profile singularity --reads /home/eriny/sandbox/test_files/cecret/reads --outdir tests -with-tower -resume

test=$1

if [ -z "$test" ]; then test="small" ; fi

if [ "$test" == "small" ]
then
  options=("reads" "single_reads" "fastas")

  for option in ${options[@]}
  do
    # defaults
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir singularity_defaults_$option

    # removed test for bamsnap and rename because of lack of interest
    # attempted bcftools and filter
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir all_on_$option \
      --bcftools_variants true \
      --filter true \
      -resume

    # removing primer trimming
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir nontrimmed_$option \
      --trimmer 'none' \
      -resume

    # changing the cleaner, aligner, and trimmer
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir toggled_$option \
      --cleaner 'fastp' \
      --trimmer 'samtools' \
      --aligner 'minimap2' \
      -resume

    # with UPHL's config
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile uphl,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir uphl_$option \
      -resume
  done

  # multifasta
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/test_files/cecret/reads \
    --single-reads /home/eriny/sandbox/test_files/cecret/single-reads \
    --fastas /home/eriny/sandbox/test_files/cecret/fastas \
    --multifastas /home/eriny/sandbox/test_files/cecret/multifasta \
    --outdir kitchen_sink

  # empty
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads doesntexit \
    --single-reads willnotexist \
    --fastas shouldntexit \
    --outdir empty

else
  # CDC's test data with relatedness
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir default_datasets \
    --relatedness true

  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile uphl,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir uphl_datasets \
    -resume \
    --relatedness true

  # CDC's test data with relatedness using nextalign
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir toggled_datasets \
    --cleaner 'fastp' \
    --trimmer 'samtools' \
    --aligner 'minimap2' \
    --relatedness true \
    --msa 'nextalign' \
    -resume
fi
