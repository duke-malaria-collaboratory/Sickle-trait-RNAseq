#!/usr/bin/env python
"""

"""

import os
import time
import pandas as pd
import argparse
from utilities import *
from joblib import Parallel, delayed

"""
This script computes a measure of dynamic dissimilarity between an rna expression profile and a corresponding protein 
expression profile using the normalized Kendall-Tau (KT) distance function between gene lists over a range of phase 
shifts. The distance between the two profiles is taken to be the minimum KT distance achieved by some phase shift. 
The fraction of randomized expression profiles which achieve a smaller minimum KT distance over the same set of phase 
shifts is also reported to estimate the significance of the computed measure of dissimilarity. 

This script is modified from prot_rna_diss.py to enable rna to rna comparisons. By Robert C. Moseley

command: python rna_rna_diss.py -t1 <ts_file_1> -t2 <ts_file_2 -c <num_curves> -s <num_shifts> -o <out_dir>
"""

def rna_rna_diss(i, gene_name, num_shifts):

    # initialize data frame
    diss_df = pd.DataFrame()

    # loop over each gene name
    print("Gene {}: {}".format(i+1, gene_name))
    start_time = time.time()
    # extract all records corresponding to the current gene
    # gene_indices = [gene_name == gene for gene in gene_names_set]


    # identify the rna profile from ts 1
    temp_rna_ts_1 = ts_1_df.loc[gene_name].values
    # normalize data
    temp_rna_ts_1 = normalize_ts(temp_rna_ts_1, kind='scale')

    # identify the rna profile from ts 2
    temp_rna_ts_2 = ts_2_df.loc[gene_name].values
    # normalize data
    temp_rna_ts_2 = normalize_ts(temp_rna_ts_2, kind='scale')

    m = len(temp_rna_ts_1)

    # compute and save dissimilarities
    # loop over shifts of rna time series 2 to test different alignments of phases
    for shift in range(num_shifts):
        diss_df.loc[gene_name, 'shift_%s' % str(shift)] = kt_diss(temp_rna_ts_1[:m - shift],
                                                                  temp_rna_ts_2[shift:])

    # find the minumum and range of kt distances over the phase shifts
    diss_df.loc[gene_name, 'kt_diss_range'] = max(diss_df.loc[gene_name,
                                                            ['shift_%s' % str(shift)
                                                             for shift in range(num_shifts)]]) - \
                                            min(diss_df.loc[gene_name,
                                                            ['shift_%s' % str(shift)
                                                             for shift in range(num_shifts)]])

    diss_df.loc[gene_name, 'kt_diss'] = min(diss_df.loc[gene_name,
                                                      ['shift_%s' % str(shift)
                                                       for shift in range(num_shifts)]])

    # generate random curves from the rna time series 1 and prna time series 2 profiles
    # and compute KT measure over phase shifts
    kt_diss_cdf = gen_takt_diss_cdf(temp_rna_ts_1, temp_rna_ts_2, num_curves, num_shifts,
                                    dict(w=1, p=1, noise1=None, noise2=None))
    diss_df.loc[gene_name, 'kt_diss_pval'] = np.searchsorted(kt_diss_cdf, diss_df.loc[gene_name, 'kt_diss']) \
                                           / num_curves

    print('completed scoring %s in %f seconds' % (gene_name, time.time() - start_time))

    return diss_df

def est_delay(diss_row, num_shifts):

    min_shift_idx = np.argmin(diss_df.loc[diss_row.name, ['shift_%s' % str(shift) for shift in range(num_shifts)]].values)
    min_shift = diss_df.loc['shift_times', 'shift_%d' % min_shift_idx]

    return min_shift

if __name__ == "__main__":

    description = "Computes a measure of dynamic dissimilarity between rna expression profiles from two time series " \
                  "datasets a using the normalized Kendall-Tau (KT) distance function between gene lists over " \
                  "a range of phase shifts."

    parser = argparse.ArgumentParser(prog="Temporal Alignment Kendall-Tau",
                                     description=description)
    # user-defined parameters
    parser.add_argument("-t1", "--ts_file_1",
                        type=str,
                        default=None,
                        required=True,
                        help=".txt time series file 1")

    parser.add_argument("-t2", "--ts_file_2",
                        type=str,
                        default=None,
                        required=True,
                        help=".txt time series file 2 ")

    parser.add_argument("-c", "--num_curves",
                        type=int,
                        default=10000,
                        required=True,
                        help="the number of permutations of each curve to perform to derive an empirical p-value")

    parser.add_argument("-s", "--num_shifts",
                        type=int,
                        default=None,
                        required=True,
                        help="number of time point shifts (including 0) to consider (match length of one cell-cycle period")

    parser.add_argument("-p", "--num_processors",
                        type=int,
                        default=0,
                        required=False,
                        help="Set to > 0 to run in parallel")

    parser.add_argument("-o", "--out_dir",
                        type=str,
                        default=".",
                        required=True,
                        help="directory to save output")

    args = parser.parse_args()

    ts_file_1 = args.ts_file_1
    ts_file_2 = args.ts_file_2
    num_curves = args.num_curves
    num_shifts = args.num_shifts
    out_dir = args.out_dir
    n_jobs = args.num_processors

    filename1 = os.path.splitext(os.path.basename(ts_file_1))[0]
    filename2 = os.path.splitext(os.path.basename(ts_file_2))[0]

    ts_1_df = pd.read_csv(ts_file_1, delimiter='\t', header=0, index_col=0, comment='#')
    ts_2_df = pd.read_csv(ts_file_2, delimiter='\t', header=0, index_col=0, comment='#')

    # get intersect of indices for unique gene names
    gene_names_set = list(set(ts_1_df.index.values) & set(ts_2_df.index.values))
    ts_1_df = ts_1_df.loc[gene_names_set]
    ts_2_df = ts_2_df.loc[gene_names_set]
    print("{} genes being analyzed".format(len(gene_names_set)+1))

    # create the output directory if needed
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    diss_df_list = Parallel(n_jobs=n_jobs)(delayed(rna_rna_diss)(i, gene_name, num_shifts)
                                           for i, gene_name in enumerate(gene_names_set))

    diss_df = pd.concat(diss_df_list)
    shift_time_int = float(ts_1_df.columns.values[1]) - float(ts_1_df.columns.values[0])
    for shift in range(num_shifts):
        diss_df.loc['shift_times', 'shift_%s' % str(shift)] = shift * shift_time_int

    # estimate delay
    diss_df['est_delay'] = diss_df.apply(est_delay, num_shifts=num_shifts, axis=1)
    val = 'shift_times'
    new_idx = [val] + diss_df.index.drop(val).tolist()
    diss_df = diss_df.reindex(new_idx)

    with open(out_dir + '/kt_diss_%s-%s.tsv' % (filename1, filename2), 'w') as outfile:
        outfile.write('# RNA vs. RNA Profile Dissimilarity Measure Over Phase Shifts \n')
        outfile.write('# Time series file 1:  %s \n' % os.path.abspath(ts_file_1))
        outfile.write('# Time series file 2:  %s \n' % os.path.abspath(ts_file_2))
        outfile.write('# No. of Randomized Curves:  %d \n' % num_curves)
        diss_df.to_csv(path_or_buf=outfile, sep='\t')
