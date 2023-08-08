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
import seaborn as sns

# This function reads the insert size matrix and splice the dataframe to include insert size between 30 and 500 bp.
# It then generates a histogram to visualise the distribution and calculate the summary statistics.
def get_distribution(file_name):
    df = pd.read_csv(file_name, header=None)
    df = df.rename(columns = {0: "insert_size"})
    df = df[(df.insert_size <= 500) & (df.insert_size >= 30)]
    # print(df[(df.insert_size > 500)].count())
    # print(df.loc[df.insert_size.between(100,200), "insert_size"].count())

    ## Compute the summary stats
    describe = df["insert_size"].describe()
    ## Change the name of the 50% index to median
    idx = describe.index.tolist()
    idx[5] = 'median'
    describe.index = idx
    pd.options.display.float_format = "{:.2f}".format
    print(describe)

    # generate histograms for prelimanry analysis
    data = df["insert_size"]

    #sns.displot(data, bins=470, kde=True, color='#FE6100')
    #plt.title("30-500bp Insert Size Distribution")
    #plt.xlabel("Insert Size")
    #plt.ylabel("Count")
    #t = datetime.datetime.now()
    #histogram_name = "Filtered_Insert_Size_Distribution_%s_%s_%s_%s_%s_%s.png" % (t.year, t.month, t.day, t.hour, t.minute, t.second)
    #plt.savefig(histogram_name, dpi=300, bbox_inches="tight")
    #plt.show()

    # Freedman–Diaconis rule to be more scientific in choosing the "right" bin width:
    q25, q75 = np.percentile(data, [25, 75])
    bin_width = 2 * (q75 - q25) * len(data) ** (-1 / 3)
    bins = round((data.max() - data.min()) / bin_width)
    print("Freedman–Diaconis number of bins:", bins)

    plt.hist(data, density= True, bins=470, alpha=0.5,
             color='#785EF0', label='Skewness:{:.2f}'.format(data.skew()))  #"#009e73", "#FE6100", #DC267F, # Alternatively #FFB000 000000 or Bang Wong #009e73 0072b2 or Paul Tol 882255 aa4499

    mn, mx = plt.xlim()
    plt.xlim(mn, mx)
    kde_xs = np.linspace(mn, mx, 470)
    # kernel-density estimate using Gaussian kernels to estimate the probability density function (PDF) of a random variable in a non-parametric way
    kde = stats.gaussian_kde(data)
    plt.plot(kde_xs, kde.pdf(kde_xs), label="probability density function")

    plt.title("30-500bp Insert Size Distribution")
    plt.xlabel("Insert Size")
    plt.ylabel("Probability")
    plt.xlim([0, 500])

    t = datetime.datetime.now()
    plt.legend()
    histogram_name = "Filtered_Insert_Size_Distribution_%s_%s_%s_%s_%s_%s.png" % (t.year, t.month, t.day, t.hour, t.minute, t.second)
    plt.savefig(histogram_name, dpi=300, bbox_inches="tight")
    plt.show()
    return

# path to the insert size.txt files generated from samtools
path = "/home/yyl23/Desktop/BS7130_Independent_Research_Project/cfDNApro/samtools_insert_size_cal"

# use os.walk() to recursively walk a directory and fnmatch.filter() to match against a simple expression
matches = []
for root, dirnames, filenames in os.walk(path):
    for filename in sorted(fnmatch.filter(filenames, "*.txt")):
        matches.append(os.path.join(root, filename))
print(matches)
for items in matches:
    get_distribution(items)
