#/bin/bash
#nextflow ~/sandbox/Cecret/main.nf -profile singularity --reads /home/eriny/sandbox/test_files/cecret/reads --outdir tests -with-tower -resume
# nextflow ~/sandbox/Cecret/main.nf -profile singularity,mpx --reads /home/eriny/sandbox/test_files/cecret/mpx --outdir tests --cleaner 'fastp' -with-tower -resume


test=$1

if [ -z "$test" ]; then test="small" ; fi

if [ "$test" == "small" ]
then
  options=("reads" "single_reads" "fastas")

  for option in ${options[@]}
  do
    # defaults
    nextflow ~/sandbox/Cecret/main.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir singularity_defaults_$option \
      -with-tower

    # removed test for bamsnap and rename because of lack of interest
    nextflow ~/sandbox/Cecret/main.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir all_on_$option \
      -with-tower \
      --filter true \
      -resume

    # removing primer trimming
    nextflow ~/sandbox/Cecret/main.nf \
      -profile singularity,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir nontrimmed_$option \
      --trimmer 'none' \
      -with-tower \
      -resume

    # changing the cleaner, aligner, and trimmer
    nextflow ~/sandbox/Cecret/main.nf \
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
    nextflow ~/sandbox/Cecret/main.nf \
      -profile uphl,artic_V3 \
      --$option /home/eriny/sandbox/test_files/cecret/$option \
      --outdir uphl_$option \
      -with-tower \
      -resume
  done

  # multifasta
  nextflow ~/sandbox/Cecret/main.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/test_files/cecret/reads \
    --single-reads /home/eriny/sandbox/test_files/cecret/single-reads \
    --fastas /home/eriny/sandbox/test_files/cecret/fastas \
    --multifastas /home/eriny/sandbox/test_files/cecret/multifasta \
    --outdir kitchen_sink  \
    -with-tower

  # empty
  nextflow ~/sandbox/Cecret/main.nf \
    -profile singularity,artic_V3 \
    --reads doesntexit \
    --single-reads willnotexist \
    --fastas shouldntexit \
    --outdir empty  \
    -with-tower

else
  # CDC's test data with relatedness
  nextflow ~/sandbox/Cecret/main.nf \
    -profile singularity,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir default_datasets \
    --relatedness true  \
    -with-tower

  nextflow ~/sandbox/Cecret/main.nf \
    -profile uphl,artic_V3 \
    --reads /home/eriny/sandbox/sars-cov-2-datasets/reads \
    --outdir uphl_datasets \
    -with-tower \
    -resume \
    --relatedness true

  # CDC's test data with relatedness using nextalign
  nextflow ~/sandbox/Cecret/main.nf \
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
  nextflow ~/sandbox/Cecret/main.nf \
    -profile singularity,mpx \
    --reads /home/eriny/sandbox/test_files/cecret/mpx \
    --outdir mpx \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume

  # MPX with idt primers
  nextflow ~/sandbox/Cecret/main.nf \
    -profile singularity,mpx_idt \
    --reads /home/eriny/sandbox/test_files/cecret/mpx_idt \
    --outdir mpx_idt \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume

  # other
  nextflow ~/sandbox/Cecret/main.nf \
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
  nextflow ~/sandbox/Cecret/main.nf \
    -profile uphl,mpx \
    --reads /home/eriny/sandbox/test_files/cecret/mpx \
    --outdir uphl_mpx \
    --cleaner 'fastp' \
    --relatedness true \
    --msa 'nextalign' \
    -with-tower \
    -resume
fi
