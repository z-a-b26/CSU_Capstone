#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import sklearn as sk
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
import seaborn as sns
from sklearn import preprocessing
from pathlib import Path
import argparse
import sys


# In[2]:


def get_data(filepath):
    """
    str filepath: filepath to the directory where tsv's are located

    Example:
    TSV file: dc1a4e95f352400ca00df3a69ad3ea84_alignment_summary.tsv
    Output: alignment_summary
    To view dataframe from that tsv file, simply type out the output name

    """
    for file in glob.glob("*.tsv"):
        filename = os.path.basename(file)
        filename = filename.split('.tsv')[0]
        undie_score = filename.index('_')
        filename = filename[undie_score + 1:]
        globals()[filename] = pd.read_csv(file, comment  = '#', sep = '\t')


# In[3]:


def median_hist(values, counts):
    expanded_vector = np.repeat(values,counts)
    result = expanded_vector.quantile(.5)
    return result

def create_rna_insert_mediansize(insert_histogram):
    my_unique_samples = set(rna_insert_histogram['name'])

    # Create dict to store samples into
    sample_to_median_mapping = {}

    for unique_sample in my_unique_samples:
        slice_df = rna_insert_histogram[rna_insert_histogram['name'] == unique_sample]
        sample_to_median_mapping[unique_sample] = median_hist(slice_df['Insert Size'], slice_df['Counts'])

    rna_insert_histogram['Median Insert Size'] = rna_insert_histogram.apply(lambda x: sample_to_median_mapping[x['name']], axis = 1)

    rna_insert_mediansize = rna_insert_histogram.loc[:,['name','Median Insert Size']]

    return rna_insert_mediansize


# In[4]:


def msd(expected, observed):
    result = sum((expected - observed)**2)
    return result

def create_rna_coverageMSD(coverage_histogram):
    my_unique_samples_again = set(rnaseq_coverage_histogram['name'])

    sample_to_median_mapping_again = {}
    for unique_sample in my_unique_samples_again:
        slice_df = rnaseq_coverage_histogram[rnaseq_coverage_histogram['name'] == unique_sample]
        sample_to_median_mapping_again[unique_sample] = msd(slice_df['Normalized Coverage'], 1)

    rnaseq_coverage_histogram['CoverageMSD'] = rnaseq_coverage_histogram.apply(lambda x: sample_to_median_mapping_again[x['name']], axis = 1)

    rna_coverageMSD = rnaseq_coverage_histogram.loc[:, ['name', 'CoverageMSD']]

    return rna_coverageMSD


# In[5]:


def create_master(alignment_summ, dup_metrics, summary_metrics, insert_mediansize, coverageMSD):
    """
    df alignment_summ: alignment_summary dataframe
    df dup_metrics: duplication_metrics dataframe
    df summary_metrics: rnaseq_summary_metrics dataframe
    df insert_mediansize: rna_insert_mediansize dataframe
    df coverageMSD: rna_coverageMSD dataframe
    """
    master = alignment_summ.merge(dup_metrics, on = 'name', how = 'left')
    master = master.merge(summary_metrics, on = 'name', how = 'left')
    master = master.merge(insert_mediansize, on = 'name', how = 'left')
    master = master.merge(coverageMSD, on = 'name', how = 'left')
    master.drop_duplicates(inplace = True)
    
     # Remove unnecessary columns from master matrix
    master = master.drop(labels = ['index_x', 'sample_id_x', 'analysis_id_x', 'analysis_id_y', 'row_x', 'row_y', 'CATEGORY', 'index_y', 'index', 'sample_id_y', 'analysis_id', 'sample_id', 'row'], axis = 1)
    master

    # X will be the master matrix without the name column
    X = master.drop('name', 1)

    # y will contain all the names from the master matrix
    y = master['name']

    # create a features array that only contains the values from each metric
    features = []
    for col in master.columns:
        features.append(col)
    # since the name won't account for any variance in PCA, remove it from features array
    features = features[1:]

    # Standardize the data
    # Seperate the features from the master dataframe (only values should exist here)
    master_matrix = master.loc[:, features].values
    master_matrix = StandardScaler().fit_transform(master_matrix)

    return master_matrix, features, X


# In[6]:


def linear_pca(master, features):
    """
    param df master: master_matrix
    param list features: features; array of metric names from previous method
    
    Will return top metrics along with heatmap for each metrics explained variance per principal component
    """
    
    ## PRINCIPAL COMPONENT ANALYSIS
    pca = PCA(n_components = .95) #len(features))
    principalComponents = pca.fit_transform(master)
    pc_df = pd.DataFrame(data = principalComponents)
    
    # bring_back_metrics will contain the explained variance ratios for each principal component
    bring_back_metrics = pd.DataFrame(pca.components_)
    bring_back_metrics = bring_back_metrics.transpose()

    # have feature names included as a column
    bring_back_metrics['features'] = features
    bring_back_metrics.set_index(keys = 'features', inplace = True)
    bring_back_metrics.index.name = None
    
    ## HEATMAP
    ax = sns.heatmap(pca.components_,
                 cmap= 'magma',# color of heatmap
                 yticklabels=[ "PCA"+str(x) for x in range(1,pca.n_components_+1)], # labeling the ticks on y axis
                 xticklabels=features, # labeling ticks on x axis
                 cbar_kws={"orientation": "vertical"}) # legend
    ax.set_aspect("equal")
    
    plt.savefig('pc_heatmap.png')
    
    ## VARIABLE IMPORTANCE
    # get absolute value of PCA scores (negative coorelation just as important as positive)
    importance = pd.DataFrame(abs(pca.components_))
    importance.columns = features

    # create a list where we store the top 3 metrics of each principal component as an element
    top_three_list = []
    for i in range(0,len(importance)):
        top_three_list.append(importance.iloc[i,:].sort_values(ascending = False).nlargest(3))
    
    # Only want the metric names (not their values)
    metric_list = []
    for pc in top_three_list:
        x = pc.to_frame()
        x.reset_index(level = 0, inplace = True)
        x = x.iloc[:,0].to_frame()
        metric_list.append(x)
    
    # Create an array with PC labels
    comp = []
    for i in range(1, pca.n_components_+1):
        pc = f"PC{i}"
        comp.append(pc)
    
    # Create the final dataframe
    top_df = pd.concat(metric_list, axis = 1)
    top_df.columns = comp
    
    these_are_the_top_metrics = top_df.stack().value_counts(ascending = False)
    these_are_the_top_metrics = these_are_the_top_metrics.to_frame()
    
    #show dataframe
    return bring_back_metrics, these_are_the_top_metrics


# In[7]:


def parse_cmdline_params():
    """
    parses command line arguments
    
    """
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parse.add_argument("-fp", "--file_path",
                       help = "input file path to tsv's",
                       type = str,
                       required = True)
    
    opts, unknown = parse.parse_known_args(args=sys.argv[1:])

    return opts


# In[8]:


def workflow(filepath): #how to test if this works
    file_names = get_data(filepath)
    mediansize = create_rna_insert_mediansize(rna_insert_histogram)
    mean_square = create_rna_coverageMSD(rnaseq_coverage_histogram)
    m_matrix = create_master(alignment_summary, duplication_metrics, rnaseq_summary_metrics, mediansize, mean_square)
    top_metrics = linear_pca(m_matrix[0], m_matrix[1])
    outpath = os.path.join(filepath, "RNASEQ_Top_Metrics.csv")
    top_metrics_df = top_metrics[1].to_csv(outpath)
    
    return top_metrics[1]
    


# In[9]:


#if __name__ == "__main__":
   # opts = parse_cmdline_params()
   # workflow(opts.file_path)


# In[10]:


workflow("~/RNA_SEQ_PCA")

