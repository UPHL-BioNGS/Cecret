#/bin/bash

sand_dir=/home/eriny/sandbox/Cecret

######################
##### paired-end #####
######################

option=reads

# defaults
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir defaults -with-tower -resume
cp defaults/summary.csv defaults_summary.csv

# attempted bcftools, filter, and bamsnap
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir all_on --bamsnap true --bcftools_variants true --filter true  -with-tower -resume
cp all_on/summary.csv all_on_summary.csv

# changing the cleaner, aligner, and trimmer
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir toggled --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2' -with-tower -resume
cp toggled/summary.csv toggled_summary.csv

# removing primer trimming
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir nontrimmed --trimmer 'none' -with-tower -resume
cp nontrimmed/summary.csv nontrimmed_summary.csv

# with UPHL's config
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/UPHL.config --$option $sand_dir/data --outdir uphl  -with-tower -resume
cp uphl/summary.csv uphl_summary.csv

# CDC's test data with relatedness
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets --relatedness true  -with-tower -resume
cp datasets/summary.csv datasets_summary.csv

# CDC's test data with relatedness using nextalign
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_nextalign --relatedness true --msa nextalign -with-tower -resume
cp datasets_nextalign/summary.csv datasets_nextalign_summary.csv

######################
##### single-end #####
######################

option=single_reads

# defaults
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir defaults_single  -with-tower -resume
cp defaults_single/summary.csv defaults_single_summary.csv

# attempted bcftools, filter, and bamsnap
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir all_on_single --bamsnap true --bcftools_variants true --filter true  -with-tower -resume
cp all_on_single/summary.csv all_on_single_summary.csv

# changing the cleaner, aligner, and trimmer
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir toggled_single --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2'  -with-tower -resume
cp toggled_single/summary.csv toggled_single_summary.csv

# with UPHL's config
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/UPHL.config --$option $sand_dir/data --outdir uphl_single -with-tower -resume
cp uphl_single/summary.csv uphl_single_summary.csv

# CDC's test data with relatedness
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_single --relatedness true -with-tower -resume
cp datasets_single/summary.csv datasets_single_summary.csv

######################
##### empty #####
######################
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --reads doesntexit --single-reads willnotexist --fastas shouldntexit --outdir empty  -with-tower -resume
