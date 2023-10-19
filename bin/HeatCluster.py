#!/usr/bin/python3

###########################################
# HeatCluster-0.4.10                      #
# written by Stephen Beckstrom-Sternberg  #
# Creates SNP heat/cluster maps           #
# from SNP matrices                       #
###########################################

import argparse
import logging
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%y-%b-%d %H:%M:%S', level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='input SNP matrix', default='snp-dists.txt')
parser.add_argument('-o', '--out', type=str, help='final file name', default='SNP_matrix')
parser.add_argument('-t', '--type', type=str, help='file extension for final image', default = 'pdf')
parser.add_argument('-v', '--version', help='print version and exit', action='version', version='%(prog)s ' + '0.4.10')
args = parser.parse_args()

def read_snp_matrix(file):
    logging.debug('Determining if file is comma or tab delimited')
    tabs   = pd.read_csv(file, nrows=1, sep='\t').shape[1]
    commas = pd.read_csv(file, nrows=1, sep=',').shape[1]
    if tabs > commas:
        logging.debug('The file is probably tab-delimited')
        df = pd.read_csv(file, sep='\t', index_col= False)
    else:
        logging.debug('The file is probably comma-delimited')
        df = pd.read_csv(file, sep=',', index_col= False)

    return df

def clean_and_read_df(df):
    """
    Clean and read DataFrame from lines.
    
    Args:
        lines (list): List of strings representing lines of data.
        
    Returns:
        df (DataFrame): Cleaned DataFrame.
    """
     # Define consensus patterns
    consensus_patterns = ['snp-dists 0.8.2', '.consensus_threshold_0.6_quality_20', 'Consensus_', 'Unnamed: 0']
    
    # Replace consensus patterns in column names
    df.columns = df.columns.str.replace('|'.join(consensus_patterns), '', regex=True)
    # Replace consensus patterns in entire dataframe to change row names
    df = df.replace(consensus_patterns, '', regex=True) 

    # Keep only numeric columns
    df = df.set_index(df.columns[0])
    df.dropna(axis=0, inplace=True)
    df.dropna(axis=1, inplace=True)  
    return df


def main():
    try:
        path = Path('./snp-dists.txt')
        path.resolve(strict=True)
    except FileNotFoundError:
        path = Path('./snp_matrix.txt')

    print("Using file path:", path)

    lines = read_snp_matrix(path)
    numSamples = len(lines) - 1

    df = clean_and_read_df(lines)


    if (numSamples) >= 140:
        fontSize = 2
    elif (numSamples) >=100:
        fontSize = 4
    elif (numSamples) >=60:
        fontSize = 6
    else:
        fontSize=8    
       
    df = df.loc[df.sum(axis=1).sort_values(ascending=True).index]
    df.replace([np.inf, -np.inf], np.nan)
    df.dropna()

    df = df.reindex(columns=df.index)
    print("df after re-indexing columns:\n\n",df,"\n\n")
    heatmap = sns.clustermap(
        df,
        xticklabels=True,
        yticklabels=True,
        vmin=0,
        vmax=80,
        center=20,
        annot=True,
        annot_kws={'size': fontSize},
        cbar_kws={"orientation": "vertical", "pad": 0.5},
        cmap='Reds_r',
        linecolor="white",
        linewidths=.1,
        fmt='d',
        col_cluster=False, 
        row_cluster=False
    )
    
# Set orientation of axes labels
    plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right',fontsize=fontSize)
    plt.setp(heatmap.ax_heatmap.get_yticklabels(), rotation='horizontal', fontsize=fontSize)
    
    plt.title('SNP matrix visualized via HeatCluster')
        
    heatmap.ax_row_dendrogram.set_visible(False)
    heatmap.ax_col_dendrogram.set_visible(False)

    heatmap.savefig('SNP_matrix.pdf')

    plt.show()
    print("Done")

if __name__ == "__main__":
    main()
