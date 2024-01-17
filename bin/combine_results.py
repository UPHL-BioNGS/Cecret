#!/bin/python

import pandas as pd
import glob
import os
import sys
from os.path import exists

# Entire sample set already in one file
aci_file                = 'aci_coverage_summary.csv'
ampliconstats_file      = 'ampliconstats.summary'
samtools_coverage_file  = 'samtools_coverage_summary.tsv'
pangolin_file           = 'multiqc_data/multiqc_pangolin.txt'
pango_collapse_file     = 'pango_collapse.csv'
nextclade_file          = 'multiqc_data/multiqc_nextclade.txt'
vadr_file               = 'vadr.csv'
fastp_file              = 'multiqc_data/multiqc_general_stats.txt'
fastq_names_file        = 'fastq_names.csv'
fastqc_file             = 'multiqc_data/multiqc_fastqc.txt'
freyja_file             = 'aggregated-freyja.tsv'
seqyclean_file          = 'multiqc_data/multiqc_seqyclean.txt'
versions_file           = 'versions.csv'

columns = []
summary_df = pd.DataFrame(columns=['sample_id'])
min_depth = int(sys.argv[1])

if exists(fastqc_file) and exists(fastq_names_file):
    print("Getting results from fastqc file " + fastqc_file)
    fastqc_df = pd.read_table(fastqc_file)
    names_df  = pd.read_csv(fastq_names_file)
    fastq_R1_df = pd.merge(names_df, fastqc_df, left_on = 'fastq_1', right_on = 'Filename', how = 'inner' )
    fastq_R1_df['fastqc_raw_reads_1'] = fastq_R1_df['Total Sequences']
    names_tmp_df = names_df.dropna(subset = ['fastq_2'])
    fastq_R2_df = pd.merge(names_tmp_df, fastqc_df, left_on = 'fastq_2', right_on = 'Filename', how = 'inner' )
    fastq_R2_df['fastqc_raw_reads_2'] = fastq_R1_df['Total Sequences']
    fastq_both_df = pd.merge(fastq_R1_df, fastq_R2_df, left_on='sample',right_on='sample', how = 'left')
    fastqc_tmp_df = fastq_both_df[['sample', 'fastqc_raw_reads_1', 'fastqc_raw_reads_2']]
    
    summary_df = pd.merge(summary_df, fastqc_tmp_df, left_on = 'sample_id', right_on = 'sample', how = 'outer' )
    summary_df['sample_id'].fillna(summary_df['sample'], inplace=True)
    summary_df.drop('sample', axis=1, inplace=True)
    columns = columns + ['fastqc_raw_reads_1', 'fastqc_raw_reads_2']

fasta_df = pd.DataFrame(columns=['fasta_sample', "fasta_line", "num_N", "num_total"])
fasta_files = glob.glob("*.fa") + glob.glob("*.fasta") + glob.glob("*.fna")
for file in fasta_files :
    print("Getting basic information from fasta " + file)
    if not ( os.stat(file).st_size==0 ) : 
        sample              = str(file).replace(".fasta", '').replace(".fna", '').replace('.fa', '').replace('.consensus','')
        with open(file) as fasta : 
            fasta_line = ''
            sequence = ''
            for line in fasta :
                if ">" in line :
                    fasta_line = line.replace(">","").strip().replace(".consensus_threshold*", "")
                else :
                    sequence = sequence + line.strip()
                    
            num_N           = sequence.count('N') + sequence.count('n')
            num_total       = len(sequence)
            tmp_fasta_df    = pd.DataFrame({'fasta_sample': [sample],  'fasta_line': [fasta_line], 'num_N': [num_N], 'num_total': [num_total]})
            fasta_df        = pd.concat([fasta_df, tmp_fasta_df], axis=0 )

if not fasta_df.empty :
    summary_df              = pd.merge(summary_df, fasta_df, left_on = 'sample_id', right_on = 'fasta_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['fasta_sample'])
    summary_df              = summary_df.drop('fasta_sample', axis=1)
    columns                 = ['fasta_line'] + columns + ['num_N', 'num_total']

summary_df['sample_id'] = summary_df['sample_id'].astype(object)

if exists(fastp_file) :
    with open(fastp_file) as file:
        contents = file.read()
        if "fastp_mqc" in contents :
            print("Getting results for fastp from " + fastp_file )
            
            fastp_df = pd.read_table(fastp_file, dtype = str, sep="\t")
            fastp_df['fastp_sample_match']  = fastp_df['Sample'].str.replace('Consensus','').str.split('_').str[0]
            fastp_df['fastp_passed_reads']  = fastp_df['fastp_mqc-generalstats-fastp-filtering_result_passed_filter_reads']
            fastp_df['fastp_pct_surviving'] = fastp_df['fastp_mqc-generalstats-fastp-pct_surviving']

            fastp_tmp_df = fastp_df[['fastp_sample_match', 'fastp_passed_reads', 'fastp_pct_surviving']]
            fastp_tmp1_df = fastp_tmp_df[~fastp_tmp_df['fastp_passed_reads'].isna()]

            fastp_columns = list(fastp_tmp1_df.columns)
            fastp_columns.remove('fastp_sample_match')

            summary_df = pd.merge(summary_df, fastp_tmp1_df, left_on = 'sample_id', right_on = 'fastp_sample_match', how = 'outer')
            summary_df['sample_id'].fillna(summary_df['fastp_sample_match'], inplace=True)
            summary_df.drop('fastp_sample_match', axis=1, inplace=True)
            columns = columns + fastp_columns

if exists(seqyclean_file) :
    print("Getting results from seqyclean file " + seqyclean_file)

    use_cols = ['Sample', 'Perc_Kept']

    first = pd.read_table(seqyclean_file, sep = '\t' , dtype = str, nrows=1)
    if 'SEReadsKept' in first.columns:
        use_cols.append('SEReadsKept')
    if 'PairsKept' in first.columns:
        use_cols.append('PairsKept')

    seqyclean_df = pd.read_table(seqyclean_file, dtype = str, sep="\t", usecols = use_cols)
    seqyclean_df = seqyclean_df.add_prefix('seqyclean_')
    seqyclean_df['seqyclean_name'] = seqyclean_df['seqyclean_Sample'].str.replace("seqyclean/", "").str.replace('_clean', '').str.replace("_cln", "")
    seqyclean_columns = list(seqyclean_df.columns)
    seqyclean_columns.remove('seqyclean_Sample')
    seqyclean_columns.remove('seqyclean_name')

    summary_df = pd.merge(summary_df, seqyclean_df, left_on = 'sample_id', right_on = 'seqyclean_name', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['seqyclean_name'], inplace=True)
    summary_df.drop('seqyclean_name', axis=1, inplace=True)
    columns = columns + seqyclean_columns

kraken2_df = pd.DataFrame(columns=['kraken2_sample', '%_human_reads', 'top_organism', 'percent_reads_top_organism'])
kraken2_files = glob.glob("*_kraken2_report.txt")
for file in kraken2_files :
    print("Getting species information from " + file)
    if not ( os.stat(file).st_size==0 ) : 
        sample = str(file).replace('_kraken2_report.txt', '')
        percent_human_reads = 0
        top_organism = 'none'
        percent_reads_top_organism = 0
        with open(file) as report : 
            kraken2_sample_df = pd.DataFrame(columns=['percent', 'reads', 'species'])
            for line in report :
                if ("S" == line.split()[3]):
                    per=line.split()[0].strip()
                    reads=int(line.split()[1].strip())
                    species=' '.join(line.split()[5:]).strip()
                    tmp_kraken2_df    = pd.DataFrame({'percent': [per], 'reads': [reads], 'species': [species]})
                    kraken2_sample_df        = pd.concat([kraken2_sample_df, tmp_kraken2_df], axis=0 )
        kraken2_sample_df = kraken2_sample_df.sort_values(by=['reads'], ascending=False)
        kraken2_human_df   = kraken2_sample_df[kraken2_sample_df['species'].isin(['Homo sapiens','Human'])]
        if not kraken2_human_df.empty :
            percent_human_reads = kraken2_human_df['percent'].iloc[0]
        
        kraken2_species_df   = kraken2_sample_df[~kraken2_sample_df['species'].isin(['Homo sapiens','Human'])]
        if not kraken2_species_df.empty :
            top_organism = kraken2_species_df['species'].iloc[0]
            percent_reads_top_organism = kraken2_species_df['percent'].iloc[0]
        
        kraken2_tmp_df = pd.DataFrame({'kraken2_sample' : [sample], '%_human_reads' : [percent_human_reads], 'top_organism' : [top_organism], 'percent_reads_top_organism' : [percent_reads_top_organism]})
        kraken2_df = pd.concat([kraken2_df, kraken2_tmp_df], axis=0 ) 

if not kraken2_df.empty :
    summary_df              = pd.merge(summary_df, kraken2_df, left_on = 'sample_id', right_on = 'kraken2_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['kraken2_sample'])
    summary_df              = summary_df.drop('kraken2_sample', axis=1)
    columns                 = columns + ['top_organism', 'percent_reads_top_organism', '%_human_reads' ]

depth_df = pd.DataFrame(columns=['samtools_sample', "num_pos_" + str(min_depth) + "X"])
depth_files = glob.glob("*.depth.txt")
for file in depth_files :
    print("Finding bases above " + str(min_depth) + " in " + file)
    if not ( os.stat(file).st_size==0 ) : 
        ind_depth_df        = pd.read_table(file, header=None)
        depth               = ind_depth_df[2][ind_depth_df[2] > min_depth].count()
        sample              = str(file).replace('.depth.txt', '')
        tmp_depth_df        = pd.DataFrame({'samtools_sample': [sample],  "num_pos_" + str(min_depth) + "X": [depth]})
        depth_df            = pd.concat([depth_df, tmp_depth_df], axis=0 )

if not depth_df.empty :
    summary_df              = pd.merge(summary_df, depth_df, left_on = 'sample_id', right_on = 'samtools_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['samtools_sample'])
    summary_df              = summary_df.drop('samtools_sample', axis=1)
    depth_columns           = list(depth_df.columns)
    depth_columns.remove('samtools_sample')
    columns                 = columns + depth_columns

aci_df = pd.DataFrame(columns=['aci_name', 'aci_num_failed_amplicons'])
if exists(aci_file) :
    print("Getting results from ACI " + aci_file)
    aci_df = pd.read_csv(aci_file, dtype = str)
    aci_df = aci_df.add_prefix('aci_')    
    tmp_df = aci_df.iloc[:,aci_df.columns != 'aci_bam'].astype('float').copy()
    aci_df['aci_num_failed_amplicons'] = tmp_df.apply(lambda x: x[x < min_depth ].count(), axis=1)
    aci_df['aci_name'] = aci_df['aci_bam'].str.replace(".primertrim.sorted.bam", "", regex = False).str.replace('.sorted.bam', '', regex = False)

    summary_df = pd.merge(summary_df, aci_df[['aci_name', 'aci_num_failed_amplicons']], left_on = 'sample_id', right_on = 'aci_name', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['aci_name'], inplace=True)
    summary_df.drop('aci_name', axis=1, inplace=True)
    columns = columns + ['aci_num_failed_amplicons']

if not ( os.stat(ampliconstats_file).st_size==0 ) : 
    print("Getting results from samtools ampliconstats file " + ampliconstats_file)
    amp_df = pd.read_table(ampliconstats_file, header=None)
    amp_df = amp_df.rename(columns = {1:'sample'})
    amp_df['samtools_num_failed_amplicons'] = amp_df.iloc[ : , 2: ].lt(20).sum( axis=1 )
    amp_df['sample_name'] = amp_df['sample'].str.replace('.primertrim.sorted', '', regex = False)
    amp_tmp_df = amp_df[['sample_name', 'samtools_num_failed_amplicons']]
    
    summary_df = pd.merge(summary_df, amp_tmp_df, left_on = 'sample_id', right_on = 'sample_name', how = 'outer' )
    summary_df['sample_id'].fillna(summary_df['sample_name'], inplace=True)
    summary_df.drop('sample_name', axis=1, inplace=True)
    columns = columns + ['samtools_num_failed_amplicons']

stats_df = pd.DataFrame(columns=['samtools_sample', 'insert_size_after_trimming'])
stats_files = glob.glob("*.stats.txt")
for file in stats_files :
    print("Getting insert size from " + file)
    if not ( os.stat(file).st_size==0 ) :
        insert_size = 0
        lines = open(file, 'r').readlines()
        for line in lines :
            if "insert size average" in line:
                insert_size = line.split("\t")[2].strip()
                break        
        sample              = str(file).replace('.stats.txt', '')
        tmp_stats_df        = pd.DataFrame({'samtools_sample': [sample], 'insert_size_after_trimming': [insert_size]})
        stats_df            = pd.concat([stats_df, tmp_stats_df], axis=0 )

if not stats_df.empty :
    summary_df              = pd.merge(summary_df, stats_df, left_on = 'sample_id', right_on = 'samtools_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['samtools_sample'])
    summary_df              = summary_df.drop('samtools_sample', axis=1)
    columns                 = columns + ['insert_size_after_trimming']

ivar_variants_df = pd.DataFrame(columns=['ivar_sample', 'ivar_num_variants_identified'])
ivar_variant_files = glob.glob("*.variants.tsv")
for file in ivar_variant_files :
    print("Counting variants in " + file)
    if not ( os.stat(file).st_size==0 ) :
        ivar_variants = 0
        lines = open(file, 'r').readlines()
        for line in lines :
            if "TRUE" in line.split("\t")[13]:
                ivar_variants += 1

        sample              = str(file).replace('.variants.tsv', '')
        tmp_ivar_df         = pd.DataFrame({'ivar_sample': [sample], 'ivar_num_variants_identified': [ivar_variants]})
        ivar_variants_df    = pd.concat([ivar_variants_df, tmp_ivar_df], axis=0 )

if not ivar_variants_df.empty :
    summary_df              = pd.merge(summary_df, ivar_variants_df, left_on = 'sample_id', right_on = 'ivar_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['ivar_sample'])
    summary_df              = summary_df.drop('ivar_sample', axis=1)
    ivar_variants_columns   = ['ivar_num_variants_identified']
    columns                 = columns + ivar_variants_columns

bcftools_variants_df = pd.DataFrame(columns=['bcftools_sample', 'bcftools_variants_identified'])
bcftools_vcf_files = glob.glob("*.vcf")
for file in bcftools_vcf_files :
    print("Counting variants in " + file)
    if not ( os.stat(file).st_size==0 ) :
        bcftools_variants = 0
        lines = open(file, 'r').readlines()
        for line in lines :
            if "#" not in line:
                bcftools_variants += 1

        sample               = str(file).replace('.vcf', '')
        tmp_bcftools_df      = pd.DataFrame({'bcftools_sample': [sample], 'bcftools_variants_identified': [bcftools_variants]})
        bcftools_variants_df = pd.concat([bcftools_variants_df, tmp_bcftools_df], axis=0 )

if not bcftools_variants_df.empty :
    summary_df              = pd.merge(summary_df, bcftools_variants_df, left_on = 'sample_id', right_on = 'bcftools_sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['bcftools_sample'])
    summary_df              = summary_df.drop('bcftools_sample', axis=1)
    bcftools_variants_cols  = ['bcftools_variants_identified']
    columns                 = columns + bcftools_variants_cols

if exists(samtools_coverage_file) : 
    print("Getting results from samtools coverage file " + samtools_coverage_file)
    scov_df = pd.read_table(samtools_coverage_file)
    # depth_after_trimming is now samtools_meandepth_after_trimming
    # 1X_coverage_after_trimming is now samtools_per_1X_coverage_after_trimming
    scov_df['samtools_meandepth_after_trimming'] = scov_df['meandepth']
    scov_df['samtools_per_1X_coverage_after_trimming'] = scov_df['coverage']
    scov_tmp_df = scov_df[['sample','samtools_meandepth_after_trimming','samtools_per_1X_coverage_after_trimming' ]]

    summary_df              = pd.merge(summary_df, scov_df, left_on = 'sample_id', right_on = 'sample', how = 'outer')
    summary_df['sample_id'] = summary_df['sample_id'].fillna(summary_df['sample'])
    summary_df              = summary_df.drop('sample', axis=1)
    columns                 = columns + ['samtools_meandepth_after_trimming', 'samtools_per_1X_coverage_after_trimming']

if exists(vadr_file) :
    print("Getting results from vadr file " + vadr_file)
    vadr_df = pd.read_csv(vadr_file, dtype = str, usecols = ['name', 'p/f', 'model', 'alerts'], index_col= False)
    vadr_df = vadr_df.add_prefix('vadr_')
    vadr_columns = list(vadr_df.columns)
    vadr_columns.remove('vadr_name')
    vadr_columns.remove('vadr_p/f')

    if 'fasta_line' in summary_df.columns.tolist():
        summary_df = pd.merge(summary_df, vadr_df, left_on = 'fasta_line', right_on = 'vadr_name', how = 'outer')
        summary_df['sample_id'].fillna(summary_df['vadr_name'], inplace=True)
        summary_df.drop('vadr_name', axis=1, inplace=True)
        columns = ['vadr_p/f'] + columns + vadr_columns
    else:
        vadr_df['sample_match'] = vadr_df['vadr_name'].str.replace('Consensus_', '', regex =  False).str.split(".").str[0]
        summary_df = pd.merge(summary_df, vadr_df, left_on = 'sample_id', right_on = 'sample_match', how = 'outer')
        summary_df['sample_id'].fillna(summary_df['sample_match'], inplace=True)
        summary_df.drop('vadr_name', axis=1, inplace=True)
        summary_df.drop('sample_match', axis=1, inplace=True)
        columns = ['vadr_p/f'] + columns + vadr_columns 

if exists(nextclade_file) :
    print("Getting results from nextclade file " + nextclade_file)

    use_cols = ['Sample', 'clade', 'qc_overallstatus', 'qc_overallscore']

    first = pd.read_table(nextclade_file, sep = '\t' , dtype = str, nrows=1)
    if 'clade_who' in first.columns:
        use_cols.append('clade_who')
    if 'outbreak' in first.columns:
        use_cols.append('outbreak')
    if 'lineage' in first.columns:
        use_cols.append('lineage')

    nextclade_df = pd.read_table(nextclade_file, sep = '\t' , dtype = str, usecols = use_cols)
    nextclade_df=nextclade_df.add_prefix('nextclade_')
    nextclade_columns = list(nextclade_df.columns)
    nextclade_df['sample_match'] = nextclade_df['nextclade_Sample'].str.replace('Consensus_', '', regex =  False)
    nextclade_columns.remove('nextclade_Sample')
    nextclade_columns.remove('nextclade_clade')

    summary_df = pd.merge(summary_df, nextclade_df, left_on = 'sample_id', right_on = 'sample_match', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['sample_match'], inplace=True)
    summary_df.drop('nextclade_Sample', axis=1, inplace=True)
    summary_df.drop('sample_match', axis = 1, inplace = True )
    columns = ['nextclade_clade'] + columns + nextclade_columns

if exists(pangolin_file) :
    print("Getting results from pangolin file " + pangolin_file)

    pangolin_df = pd.read_table(pangolin_file, dtype = str)
    pangolin_df = pangolin_df.add_prefix('pangolin_')
    pangolin_columns = list(pangolin_df.columns)
    pangolin_df['sample_match'] = pangolin_df['pangolin_Sample'].str.replace('Consensus_', '', regex= False)
    pangolin_columns.remove('pangolin_Sample')
    pangolin_columns.remove('pangolin_lineage')

    summary_df = pd.merge(summary_df, pangolin_df, left_on = 'sample_id', right_on = 'sample_match', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['sample_match'], inplace=True)
    summary_df.drop('pangolin_Sample', axis=1, inplace=True)
    summary_df.drop('sample_match', axis=1, inplace=True)
    columns = ['pangolin_lineage'] + columns + pangolin_columns

if exists(pango_collapse_file) :
    print("Getting results from pango collapse file " + pango_collapse_file)

    pangocollapse_df = pd.read_csv(pango_collapse_file, dtype = str, usecols=['lineage','Lineage_full', 'Lineage_expanded', 'Lineage_family'])
    pangocollapse_df = pangocollapse_df.add_prefix('pangocollapse_')
    pangocollapse_columns = list(pangocollapse_df.columns)

    summary_df = pd.merge(summary_df, pangocollapse_df, left_on = 'pangolin_lineage', right_on = 'pangocollapse_lineage', how = 'left')
    columns = columns + pangocollapse_columns

if exists(freyja_file) :
    print("Getting results from freyja file " + freyja_file)
    freyja_df = pd.read_table(freyja_file, dtype = str, sep="\t")
    freyja_df = freyja_df.add_prefix('freyja_')
    freyja_df['freyja_Unnamed: 0'] = freyja_df['freyja_Unnamed: 0'].str.replace("_variants.tsv", "", regex = False)
    freyja_columns = ['freyja_summarized']
    
    summary_df = pd.merge(summary_df, freyja_df, left_on = 'sample_id', right_on = 'freyja_Unnamed: 0', how = 'outer')
    summary_df['sample_id'].fillna(summary_df['freyja_Unnamed: 0'], inplace=True)
    summary_df.drop('freyja_Unnamed: 0', axis=1, inplace=True)
    columns = columns + freyja_columns

if exists(versions_file) :
    print("Adding versions to summary file from " + versions_file)
    versions_df = pd.read_csv(versions_file, dtype = str)
    version_columns = list(versions_df.columns)

    for version in version_columns:
        software_version = versions_df[version].iloc[0]
        summary_df[version] = software_version
    columns = columns + version_columns

summary_df['sample']    = summary_df['sample_id'].str.split("_").str[0]
summary_df['sample_id'] = summary_df['sample_id'].astype('string')
summary_df = summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True)
summary_df = summary_df.sort_values(by=['sample_id'], ascending=True)
summary_df = summary_df.replace([" ", ",", "\t", "\n"], [" ", " ", " ", " "], regex=True)
summary_df = summary_df.drop_duplicates(keep='first')
summary_df.to_csv('cecret_results.csv', columns = ['sample_id','sample'] + columns, index=False)
summary_df.to_csv('cecret_results.txt', columns = ['sample_id','sample'] + columns, index=False, sep = "\t")