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
      --outdir singularity_defaults_$option \
      -with-tower

    # removed test for bamsnap and rename because of lack of interest
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir all_on_$option \
      -with-tower \
      --filter true \
      -resume

    # removing primer trimming
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir nontrimmed_$option \
      --trimmer 'none' \
      -with-tower \
      -resume

    # changing the cleaner, aligner, and trimmer
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir toggled_$option \
      --cleaner 'fastp' \
      --trimmer 'samtools' \
      --aligner 'minimap2' \
      -with-tower \
      -resume

    # with UPHL's config
    nextflow ~/sandbox/Cecret/Cecret.nf \
      -profile uphl,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir uphl_$option \
      -with-tower \
      -resume
  done

  # multifasta
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/test_files/cecret/reads \
    --single-reads /home/eriny/sandbox/test_files/cecret/single-reads \
    --fastas /home/eriny/sandbox/test_files/cecret/fastas \
    --multifastas /home/eriny/sandbox/test_files/cecret/multifasta \
    --outdir kitchen_sink  \
    -with-tower

  # empty
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads doesntexit \
    --single-reads willnotexist \
    --fastas shouldntexit \
    --outdir empty  \
    -with-tower

else
  # CDC's test data with relatedness
  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir default_datasets \
    --relatedness true  \
    -with-tower

  nextflow ~/sandbox/Cecret/Cecret.nf \
    -profile uphl,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir uphl_datasets \
    -with-tower \
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
    -with-tower \
    -resume
fi
