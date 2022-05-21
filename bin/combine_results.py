#!/bin/python

import pandas as pd
from os.path import exists

pangolin_file='lineage_report.csv'
nextclade_file='nextclade.csv'
vadr_file='vadr.csv'
summary_file='combined_summary.csv'
seqyclean_file

summary_df = pd.read_csv(summary_file, dtype = str)

if exists(pangolin_file) :
    pangolin_df = pd.read_csv(pangolin_file, dtype = str, usecols = ['taxon', 'status', 'scorpio_call', 'version', 'pangolin_version', 'pango_version', 'pangoLEARN_version'])
    pangolin_df['simple_name'] = pangolin_df['taxon'].replace('Consensus_','', regex=True).replace('.consensus.*', '', regex=True)
    pangolin_df=pangolin_df.add_prefix('pangolin_')

    summary_df = pd.merge(summary_df, pangolin_df, left_on = 'sample', right_on = 'pangolin_simple_name', how = 'outer')

    summary_df['sample_id'] = summary_df.loc[summary_df['sample_id'].ne('null'), 'pangolin_taxon']
    summary_df['sample'] = summary_df.loc[summary_df['sample'].ne('null'), 'pangolin_simple_name']

    summary_df.drop('pangolin_simple_name', axis=1, inplace=True)
    summary_df.drop('pangolin_taxon', axis=1, inplace=True)

if exists(nextclade_file) :
    nextclade_df = pd.read_csv(nextclade_file, sep = ';' , dtype = str, usecols = ['seqName', 'clade', 'qc.overallStatus', 'qc.overallScore'])
    nextclade_df['simple_name'] = nextclade_df['seqName'].replace('Consensus_','', regex=True).replace('.consensus.*', '', regex=True)
    nextclade_df=nextclade_df.add_prefix('nextclade_')

    summary_df = pd.merge(summary_df, nextclade_df, left_on = 'sample', right_on = 'nextclade_simple_name', how = 'outer')

    summary_df['sample_id'] = summary_df.loc[summary_df['sample_id'].ne('null'), 'nextclade_seqName']
    summary_df['sample'] = summary_df.loc[summary_df['sample'].ne('null'), 'nextclade_simple_name']

    summary_df.drop('nextclade_simple_name', axis=1, inplace=True)
    summary_df.drop('nextclade_seqName', axis=1, inplace=True)

if exists(vadr_file) :
    vadr_df = pd.read_csv(vadr_file, dtype = str, usecols = ['name', 'p/f', 'model', 'alerts'])
    vadr_df['simple_name'] = vadr_df['name'].replace('Consensus_','', regex=True).replace('.consensus.*', '', regex=True)
    vadr_df=vadr_df.add_prefix('vadr_')

    summary_df = pd.merge(summary_df, vadr_df, left_on = 'sample', right_on = 'vadr_simple_name', how = 'outer')

    # This currently erases the columns that were there for some reason
    #summary_df['sample_id'] = summary_df.loc[summary_df['sample_id'].ne('null'), 'vadr_name']
    #summary_df['sample'] = summary_df.loc[summary_df['sample'].ne('null'), 'vadr_simple_name']

    summary_df.drop('vadr_simple_name', axis=1, inplace=True)
    summary_df.drop('vadr_name', axis=1, inplace=True)

summary_df.to_csv('cecret_results.csv')
