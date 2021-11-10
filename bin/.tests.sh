#/bin/bash

sand=/home/eriny/sandbox/Cecret

######################
##### paired-end #####
######################

options=("reads" "single_reads" "fastas")

for option in ${options[@]}
do
  # defaults
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option $sand/data --outdir defaults_$option -with-tower
  cp defaults_$option/summary.csv defaults_summary_$option.csv

  # attempted bcftools, filter, and bamsnap
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option $sand/data --outdir all_on_$option --bamsnap true --bcftools_variants true --filter true  -with-tower -resume
  cp all_on_$option/summary.csv all_on_summary_$option.csv

  # changing the cleaner, aligner, and trimmer
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option $sand/data --outdir toggled_$option --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2' -with-tower -resume
  cp toggled_$option/summary.csv toggled_summary_$option.csv

  # removing primer trimming
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option $sand/data --outdir nontrimmed_$option --trimmer 'none' -with-tower -resume
  cp nontrimmed_$option/summary.csv nontrimmed_summary_$option.csv

  # with UPHL's config
  nextflow $sand/Cecret.nf -c $sand/configs/UPHL.config --$option $sand/data --outdir uphl_$option -with-tower -resume
  cp uphl_$option/summary.csv uphl_summary_$option.csv

  # CDC's test data with relatedness
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_$option --relatedness true  -with-tower -resume
  cp datasets_$option/summary.csv datasets_summary_$option.csv

  # CDC's test data with relatedness using nextalign
  nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_nextalign_$option --relatedness true --msa nextalign -with-tower -resume
  cp datasets_nextalign_$option/summary.csv datasets_nextalign_summary_$option.csv
done

#################
##### empty #####
#################
nextflow $sand/Cecret.nf -c $sand/configs/singularity.config --reads doesntexit --single-reads willnotexist --fastas shouldntexit --outdir empty  -with-tower
