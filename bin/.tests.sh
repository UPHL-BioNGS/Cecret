#/bin/bash

sand_dir=/home/eriny/sandbox/Cecret

######################
##### paired-end #####
######################

option=reads

# defaults
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir defaults  -with-tower
cp defaults/summary.csv defaults_summary.csv

# attempted bcftools, filter, and bamsnap
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir all_on --bamsnap true --bcftools true --filter true  -with-tower
cp all_on/summary.csv all_on_summary.csv

# changing the cleaner, aligner, and trimmer
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir toggled --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2' -with-tower
cp toggled/summary.csv toggled_summary.csv

# with UPHL's config
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/UPHL.config --$option $sand_dir/data --outdir uphl  -with-tower
cp uphl/summary.csv uphl_summary.csv

# CDC's test data with relatedness
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets --relatedness true  -with-tower
cp datasets/summary.csv datasets_summary.csv

######################
##### single-end #####
######################

option=single-reads

# defaults
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir defaults_single  -with-tower
cp defaults_single/summary.csv defaults_single_summary.csv

# attempted bcftools, filter, and bamsnap
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir all_on_single --bamsnap true --bcftools true --filter true  -with-tower
cp all_on_single/summary.csv all_on_single_summary.csv

# changing the cleaner, aligner, and trimmer
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option $sand_dir/data --outdir toggled_single --cleaner 'fastp' --trimmer 'samtools' --aligner 'minimap2'  -with-tower
cp toggled_single/summary.csv toggled_single_summary.csv

# with UPHL's config
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/UPHL.config --$option $sand_dir/data --outdir uphl_single -with-tower
cp uphl_single/summary.csv uphl_single_summary.csv

# CDC's test data with relatedness
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --$option /home/eriny/sandbox/sars-cov-2-datasets/reads/ --outdir datasets_single --relatedness true -with-tower
cp datasets_single/summary.csv datasets_single_summary.csv

######################
##### empty #####
######################
nextflow $sand_dir/Cecret.nf -c $sand_dir/configs/singularity.config --reads doesntexit --single-reads willnotexist --fastas shouldntexit --outdir empty  -with-tower
