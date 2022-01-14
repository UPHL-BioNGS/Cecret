#!/bin/python

import pandas as pd
from os.path import exists

pangolin_file='lineage_report.csv'
nextclade_file='nextclade.csv'
vadr_file='vadr.csv'
summary_file='combined_summary.csv'

summary_df = pd.read_csv(summary_file, dtype = str)
columns = list(summary_df.columns)
columns.remove('sample_id')
columns.remove('sample')
columns.remove('fasta_line')

if exists(vadr_file) :
    vadr_df = pd.read_csv(vadr_file, dtype = str, usecols = ['name', 'p/f', 'model', 'alerts'], index_col= False)
    vadr_df=vadr_df.add_prefix('vadr_')
    vadr_columns = list(vadr_df.columns)
    vadr_columns.remove('vadr_name')
    vadr_columns.remove('vadr_p/f')

    summary_df = pd.merge(summary_df, vadr_df, left_on = 'fasta_line', right_on = 'vadr_name', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['vadr_name'], inplace=True)
    summary_df['fasta_line'].fillna(summary_df['vadr_name'], inplace=True)
    summary_df.drop('vadr_name', axis=1, inplace=True)
    columns = ['vadr_p/f'] + columns + vadr_columns

if exists(nextclade_file) :
    nextclade_df = pd.read_csv(nextclade_file, sep = ';' , dtype = str, usecols = ['seqName', 'clade', 'qc.overallStatus', 'qc.overallScore'])
    nextclade_df=nextclade_df.add_prefix('nextclade_')
    nextclade_columns = list(nextclade_df.columns)
    nextclade_columns.remove('nextclade_seqName')
    nextclade_columns.remove('nextclade_clade')

    summary_df = pd.merge(summary_df, nextclade_df, left_on = 'fasta_line', right_on = 'nextclade_seqName', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['nextclade_seqName'], inplace=True)
    summary_df['fasta_line'].fillna(summary_df['nextclade_seqName'], inplace=True)
    summary_df.drop('nextclade_seqName', axis=1, inplace=True)
    columns = ['nextclade_clade'] + columns + nextclade_columns

if exists(pangolin_file) :
    pangolin_df = pd.read_csv(pangolin_file, dtype = str, usecols = ['taxon', 'lineage', 'status', 'scorpio_call', 'version', 'pangolin_version', 'pango_version', 'pangoLEARN_version'])
    pangolin_df=pangolin_df.add_prefix('pangolin_')
    pangolin_columns = list(pangolin_df.columns)
    pangolin_columns.remove('pangolin_taxon')
    pangolin_columns.remove('pangolin_lineage')

    summary_df = pd.merge(summary_df, pangolin_df, left_on = 'fasta_line', right_on = 'pangolin_taxon', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['pangolin_taxon'], inplace=True)
    summary_df['fasta_line'].fillna(summary_df['pangolin_taxon'], inplace=True)
    summary_df.drop('pangolin_taxon', axis=1, inplace=True)
    columns = ['pangolin_lineage'] + columns + pangolin_columns

summary_df.to_csv('cecret_results.csv', columns = ['sample_id','sample'] + columns, index=False)
