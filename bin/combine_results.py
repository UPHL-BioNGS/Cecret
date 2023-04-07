#!/bin/python

import pandas as pd
import glob
import os
from os.path import exists

# Entire sample set already in one file
ampliconstats_file      = 'ampliconstats.summary'
samtools_coverage_file  = 'samtools_coverage_summary.tsv'
pangolin_file           = 'lineage_report.csv'
nextclade_file          = 'nextclade.csv'
vadr_file               = 'vadr.csv'
freyja_file             = 'aggregated-freyja.tsv'
seqyclean_file          = 'Combined_SummaryStatistics.tsv'
seqycln_file            = 'Combined_seqyclean_SummaryStatistics.tsv'

# One sample per file files

# samtools ampliconstats


#num_failed_amplicons=$(grep ^FREADS samtools_ampliconstats/!{sample}_ampliconstats.txt | cut -f 2- | tr '\t' '\n' | awk '{ if ($1 < 20) print $0 }' | wc -l)
#    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi


#if (len(ampliconstats_files) > 0 ):

#depth=$(awk '{ if ($3 > !{params.minimum_depth} ) print $0 }' samtools_depth/!{sample}.depth.txt | wc -l )
#    if [ -z "$depth" ] ; then depth="0" ; fi

# 2249693-IA-M05216-230323.depth.txt
# 2249693-IA-M05216-230323.multicov.txt
# 2249693-IA-M05216-230323.stats.txt
# 2249693-IA-M05216-230323.variants.tsv
# 2249693-IA-M05216-230323.vcf

depth = list()
depth_files = glob.glob("*.depth.txt")
# TO DO : change 20 to minimum depth set by workflow
for file in depth_files:
    print(file)
    if not (os.stat(file).st_size==0) : 
        ind_depth_df = pd.read_table(file, header=None)
        depth = ind_depth_df[2][ind_depth_df[2] > 20].count()
        sample = str(file).replace('.depth.txt', '')
        print(str(depth) + " and " + sample)
        new_values = [sample, depth]
        print(new_values)
        depth.append(new_values)

print(depth)

if exists(samtools_coverage_file) : 
    print("Getting results from samtools coverage file " + samtools_coverage_file)
    scov_df = pd.read_table(samtools_coverage_file)
    #print(scov_df)

if exists(ampliconstats_file) :
    print("Getting results from samtools ampliconstats file " + ampliconstats_file)
    amp_df = pd.read_table(ampliconstats_file, header=None)
    amp_df = amp_df.rename(columns = {1:'sample'})
    amp_df['failed_amplicons'] = amp_df.iloc[:,2:].lt(20).sum(axis=1)
    amp_df = amp_df.add_prefix('ampliconstats_')
    #print(amp_df)


# if exists(vadr_file) :
#     print("Getting results from vadr file " + vadr_file)
#     vadr_df = pd.read_csv(vadr_file, dtype = str, usecols = ['name', 'p/f', 'model', 'alerts'], index_col= False)
#     vadr_df=vadr_df.add_prefix('vadr_')
#     vadr_columns = list(vadr_df.columns)
#     vadr_columns.remove('vadr_name')
#     vadr_columns.remove('vadr_p/f')

#     summary_df = pd.merge(summary_df, vadr_df, left_on = 'fasta_line', right_on = 'vadr_name', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['vadr_name'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['vadr_name'], inplace=True)
#     summary_df.drop('vadr_name', axis=1, inplace=True)
#     columns = ['vadr_p/f'] + columns + vadr_columns

# if exists(nextclade_file) :
#     print("Getting results from nextclade file " + nextclade_file)

#     use_cols = ['seqName', 'clade', 'qc.overallStatus', 'qc.overallScore']

#     first = pd.read_csv(nextclade_file, sep = ';' , dtype = str, nrows=1)
#     if 'clade_who' in first.columns:
#         use_cols.append('clade_who')
#     if 'outbreak' in first.columns:
#         use_cols.append('outbreak')
#     if 'lineage' in first.columns:
#         use_cols.append('lineage')

#     nextclade_df = pd.read_csv(nextclade_file, sep = ';' , dtype = str, usecols = use_cols)
#     nextclade_df=nextclade_df.add_prefix('nextclade_')
#     nextclade_columns = list(nextclade_df.columns)
#     nextclade_columns.remove('nextclade_seqName')
#     nextclade_columns.remove('nextclade_clade')

#     summary_df = pd.merge(summary_df, nextclade_df, left_on = 'fasta_line', right_on = 'nextclade_seqName', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['nextclade_seqName'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['nextclade_seqName'], inplace=True)
#     summary_df.drop('nextclade_seqName', axis=1, inplace=True)
#     columns = ['nextclade_clade'] + columns + nextclade_columns

# if exists(pangolin_file) :
#     print("Getting results from pangolin file " + pangolin_file)

#     pangolin_df = pd.read_csv(pangolin_file, dtype = str)
#     pangolin_df=pangolin_df.add_prefix('pangolin_')
#     pangolin_columns = list(pangolin_df.columns)
#     pangolin_columns.remove('pangolin_taxon')
#     pangolin_columns.remove('pangolin_lineage')

#     summary_df = pd.merge(summary_df, pangolin_df, left_on = 'fasta_line', right_on = 'pangolin_taxon', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['pangolin_taxon'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['pangolin_taxon'], inplace=True)
#     summary_df.drop('pangolin_taxon', axis=1, inplace=True)
#     columns = ['pangolin_lineage'] + columns + pangolin_columns

# if exists(freyja_file) :
#     print("Getting results from freyja file " + freyja_file)
#     freyja_df = pd.read_table(freyja_file, dtype = str, sep="\t")
#     freyja_df = freyja_df.add_prefix('freyja_')
#     freyja_df['freyja_Unnamed: 0'] = freyja_df['freyja_Unnamed: 0'].str.replace("_variants.tsv", "")
#     freyja_columns = list(freyja_df.columns)
#     freyja_columns.remove('freyja_Unnamed: 0')

#     summary_df = pd.merge(summary_df, freyja_df, left_on = 'sample_id', right_on = 'freyja_Unnamed: 0', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['freyja_Unnamed: 0'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['freyja_Unnamed: 0'], inplace=True)
#     summary_df.drop('freyja_Unnamed: 0', axis=1, inplace=True)
#     columns = columns + freyja_columns

# if exists(seqyclean_file) :
#     print("Getting results from seqyclean file " + seqyclean_file)
#     seqyclean_df = pd.read_table(seqyclean_file, dtype = str, sep="\t", usecols = ['OutputPrefix', 'PairsKept', 'Perc_Kept'])
#     seqyclean_df = seqyclean_df.add_prefix('seqyclean_')
#     seqyclean_df['seqyclean_OutputPrefix'] = seqyclean_df['seqyclean_OutputPrefix'].str.replace("seqyclean/", "")
#     seqyclean_df['seqyclean_OutputPrefix'] = seqyclean_df['seqyclean_OutputPrefix'].str.replace("_clean", "")
#     seqyclean_columns = list(seqyclean_df.columns)
#     seqyclean_columns.remove('seqyclean_OutputPrefix')

#     summary_df = pd.merge(summary_df, seqyclean_df, left_on = 'sample_id', right_on = 'seqyclean_OutputPrefix', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['seqyclean_OutputPrefix'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['seqyclean_OutputPrefix'], inplace=True)
#     summary_df.drop('seqyclean_OutputPrefix', axis=1, inplace=True)
#     columns = columns + seqyclean_columns

# if exists(seqycln_file) :
#     print("Getting results from seqyclean file " + seqycln_file)
#     seqycln_df = pd.read_table(seqycln_file, dtype = str, sep="\t", usecols = ['OutputPrefix', 'SEReadsKept', 'Perc_Kept'])
#     seqycln_df = seqycln_df.add_prefix('seqycln_')
#     seqycln_df['seqycln_OutputPrefix'] = seqycln_df['seqycln_OutputPrefix'].str.replace("seqyclean/", "")
#     seqycln_df['seqycln_OutputPrefix'] = seqycln_df['seqycln_OutputPrefix'].str.replace("_clean", "")
#     seqycln_columns = list(seqycln_df.columns)
#     seqycln_columns.remove('seqycln_OutputPrefix')

#     summary_df = pd.merge(summary_df, seqycln_df, left_on = 'sample_id', right_on = 'seqycln_OutputPrefix', how = 'outer')
#     summary_df['sample_id'].fillna(summary_df['seqycln_OutputPrefix'], inplace=True)
#     summary_df['fasta_line'].fillna(summary_df['seqycln_OutputPrefix'], inplace=True)
#     summary_df.drop('seqycln_OutputPrefix', axis=1, inplace=True)
#     columns = columns + seqycln_columns

# #summary_df.dropna(axis=1, how='all', inplace=True)
# summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True, inplace=True)
# summary_df.sort_values(by=['sample_id'], ascending=True)
# summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True, inplace=True)
# summary_df.drop_duplicates(keep='first', inplace=True)
# summary_df.to_csv('cecret_results.csv', columns = ['sample_id','sample'] + columns, index=False)
