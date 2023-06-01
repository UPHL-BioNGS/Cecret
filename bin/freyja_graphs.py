#!/bin/python

##############################################################################
# written by Pooja Gupta at https://github.com/UPHL-BioNGS/Cecret/issues/176 #
# adapted by Erin Young                                                      #
##############################################################################

import pandas as pd

freyja_file = 'freyja/aggregated-freyja.tsv'

# Function to make pie charts for each sample
def makePieCharts_simple(agg_df, lineages, outputFnBase):
    # Loop over all samples in the DataFrame
    for k in range(0, agg_df.shape[0]):
        sample                  = agg_df['sample'].iloc[k]
        sample_df               = agg_df[agg_df['sample' ] == sample ].copy()
        try:
            sample_df['lineages']   = sample_df['lineages'].str.split(" ")
            sample_df['abundances'] = sample_df['abundances'].str.split(" ")
            sample_df               = sample_df.explode(['lineages', 'abundances'])
        except:
            print(sample_df)
        
        sample_df['abundances'] = sample_df['abundances'].astype(float)
        
        # reducing noise to lineages with 3% or more
        sample_df               = sample_df[sample_df['abundances' ] >= 0.03 ]
        other_adundance         = 1 - sample_df['abundances' ].sum()
        if other_adundance > 0.001 :
            other_dict          = {'sample': [sample], 'lineages': ['other'], 'abundances': [other_adundance] }
            other_df            = pd.DataFrame(other_dict)
            sample_df           = pd.concat([sample_df, other_df ], ignore_index = True)
        
        # adding in some spacing
        sample_df['explode']    = 0.025
       
        # creating the pie chart
        pie_df                  = sample_df.drop('sample', axis = 1)
        pie_df                  = pie_df.set_index('lineages')
        plot = pie_df.plot.pie(y='abundances', figsize=(5, 5), autopct='%.1f%%', startangle=90, legend=False, explode=pie_df['explode'].tolist())
        plot.yaxis.set_visible(False)
        plot.get_figure().savefig('freyja/' + sample + '_freyja_lineages.png')
        plot.get_figure().savefig('freyja/' + sample + '_freyja_lineages_mqc.png')

# Read the current file into a DataFrame
df           = pd.read_table(freyja_file)
df['sample'] = df['Unnamed: 0'].str.replace("_variants.tsv","", regex=False)
df           = df.drop('Unnamed: 0', axis=1)
df           = df.drop('summarized', axis=1)
df           = df.drop('resid',      axis=1)
df           = df.drop('coverage',   axis=1)
df           = df[['sample', 'lineages', 'abundances']]

# Create pie charts for each sample in the DataFrame
makePieCharts_simple(df, 'lineages', 'freyja_sublin_pie_')