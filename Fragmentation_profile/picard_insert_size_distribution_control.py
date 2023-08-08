#!/usr/bin/env python
# coding: utf-8

"""
***************************************************************************************************************
Authors: -
Credits: -
Modified Date: 1 June 2023, 1500 BST
license: -
Version: 1.0

This script requires a working Python installation, with the following packages installed (or updated):
python-3.9.16
pandas-2.0.0
numpy-1.24.2
scipy-1.10.1
seaborn-0.12.2
matplotlib-3.7.1
Please pip install -U/conda install -c conda-forge the respective packages or update it where you find appropriate
with !pip install --upgrade (package_name)==(version).

***************************************************************************************************************
"""
# Data wrangling
import fnmatch
import os
import re, sys, datetime
import pandas as pd
import numpy as np
import scipy.stats as stats

# Data visualisation
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

def fre_table_transformer(df, freq_col, cases_cols):
    def itcases():
        for i, row in df.iterrows():
            for j in range(int(row[freq_col])):
                yield row[cases_cols]
    return pd.DataFrame(itcases())

# calculates decriptive statistics for frequency distributions
def descriptives_from_agg(values, freqs):
    values = np.array(values)
    freqs = np.array(freqs)
    arg_sorted = np.argsort(values)
    values = values[arg_sorted]
    freqs = freqs[arg_sorted]
    count = freqs.sum()
    fx = values * freqs
    mean = fx.sum() / count
    variance = ((freqs * values**2).sum() / count) - mean**2
    variance = count / (count - 1) * variance  # dof correction for sample variance
    std = np.sqrt(variance)
    minimum = np.min(values)
    maximum = np.max(values)
    cumcount = np.cumsum(freqs)
    Q1 = values[np.searchsorted(cumcount, 0.25*count)]
    Q2 = values[np.searchsorted(cumcount, 0.50*count)]
    Q3 = values[np.searchsorted(cumcount, 0.75*count)]
    idx = ['count', 'mean', 'std', 'min', '25%', 'median', '75%', 'max']
    pd.options.display.float_format = '{:.2f}'.format
    result = pd.DataFrame([count, mean, std, minimum, Q1, Q2, Q3, maximum], index=idx)
    return result

# This function reads the insert size matrix and splice the dataframe to include insert size between 30 and 500 bp.
# It then generates a histogram to visualise the distribution and calculate the summary statistics.
def get_distribution(file_name):
    # Assume every single txt files follows default picard CollectinsertSize output and present freqeuncy distribution
    df = pd.read_csv(file_name,sep="\t", header=8)

    print(os.path.splitext(os.path.basename(file_name))[0])
    summary = descriptives_from_agg(df["insert_size"], df["All_Reads.fr_count"])
    print(summary)

    ax = sns.barplot(data=df, x="insert_size", y="All_Reads.fr_count")
    plt.title("Healthy Control Insert Size Distribution")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

    plt.table(cellText=summary.values.round(2),
              rowLabels=summary.index,
              colLabels=None,
              cellLoc='right', rowLoc='center',
              loc='right', bbox=[.75, .45, .22, .5]).auto_set_column_width(col=list(range(len(summary.columns))))   # left, bottom, width, height

    t = datetime.datetime.now()
    histogram_name = file_name + "_Distribution_%s_%s_%s_%s.png" % (t.year, t.month, t.day, t.hour)
    plt.savefig(histogram_name, dpi=300, bbox_inches="tight")
    #plt.show()
    return

# path to the insert size.txt files generated from samtools
path = "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/healthy_control"
#path = "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/PTCL_combined"
#path = "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/PBMC_combined"

# use os.walk() to recursively walk a directory and fnmatch.filter() to match against a simple expression
matches = []
for root, dirnames, filenames in os.walk(path):
    for filename in sorted(fnmatch.filter(filenames, "*.txt")):
        matches.append(os.path.join(root, filename))
print(matches)
for items in matches:
    get_distribution(items)
