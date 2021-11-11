#/bin/bash

options=("reads" "single_reads" "fastas")

for option in ${options[@]}
do
  # defaults with docker
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile docker --$option /home/eriny/sandbox/test_files/cecret/$option --outdir docker_defaults_$option -with-tower
  cp docker_defaults_$option/summary.csv docker_defaults_summary_$option.csv

  # defaults
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity --$option /home/eriny/sandbox/test_files/cecret/$option --outdir defaults_$option -with-tower
  cp defaults_$option/summary.csv defaults_summary_$option.csv

  # attempted bcftools, filter, and bamsnap
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity  --$option /home/eriny/sandbox/test_files/cecret/$option --outdir all_on_$option --bamsnap true --bcftools_variants true --filter true  -with-tower -resume
  cp all_on_$option/summary.csv all_on_summary_$option.csv

  # changing the cleaner, aligner, and trimmer
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity  --$option /home/eriny/sandbox/test_files/cecret/$option --outdir toggled_$option --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2' -with-tower -resume
  cp toggled_$option/summary.csv toggled_summary_$option.csv

  # removing primer trimming
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity  --$option /home/eriny/sandbox/test_files/cecret/$option --outdir nontrimmed_$option --trimmer 'none' -with-tower -resume
  cp nontrimmed_$option/summary.csv nontrimmed_summary_$option.csv

  # with UPHL's config
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile uphl --$option /home/eriny/sandbox/test_files/cecret/$option --outdir uphl_$option -with-tower -resume
  cp uphl_$option/summary.csv uphl_summary_$option.csv

  # CDC's test data with relatedness
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_$option --relatedness true  -with-tower
  cp datasets_$option/summary.csv datasets_summary_$option.csv

  # CDC's test data with relatedness using nextalign
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_nextalign_$option --relatedness true --msa nextalign -with-tower -resume
  cp datasets_nextalign_$option/summary.csv datasets_nextalign_summary_$option.csv

  # CDC's test data with relatedness using nextalign
  nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile uphl --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_uphl_$option -with-tower -resume
  cp datasets_uphl_$option/summary.csv datasets_uphl_summary_$option.csv
done

#################
##### empty #####
#################
nextflow run UPHL-BioNGS/Cecret -r erin-dev -profile singularity --reads doesntexit --single-reads willnotexist --fastas shouldntexit --outdir empty  -with-tower
