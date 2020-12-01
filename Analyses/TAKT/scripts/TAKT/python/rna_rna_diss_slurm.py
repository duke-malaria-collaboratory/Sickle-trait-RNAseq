#!/usr/bin/env python
"""

"""

import os
import time
import pandas as pd
import argparse
import itertools
from utilities import *

"""
This script computes a measure of dynamic dissimilarity between two rna expression profiles using the normalized 
Kendall-Tau (KT) distance function between gene lists over a range of phase shifts. The distance between the two 
profiles is taken to be the minimum KT distance achieved by some phase shift. The fraction of randomized expression 
profiles which achieve a smaller minimum KT distance over the same set of phase shifts is also reported to estimate 
the significance of the computed measure of dissimilarity. 

This script is modified from prot_rna_diss.py to enable rna to rna comparisons. It has also been modified 
to work on compute clusters that use Slurm for job management. By Robert C. Moseley

*************************************************************************************
All dataframes must have a the same genes in them and contain no duplicate gene names.
*************************************************************************************

command: python rna_rna_diss_slurm.py -fd <ts_file_dir> -c <num_curves> -s <num_shifts> -o <out_dir>
"""

if __name__ == "__main__":

    description = "Computes a measure of dynamic dissimilarity between rna expression profiles from two time series " \
                  "datasets a using the normalized Kendall-Tau (KT) distance function between gene lists over " \
                  "a range of phase shifts."

    parser = argparse.ArgumentParser(prog="Temporal Alignment Kendall-Tau",
                                     description=description)
    # user-defined parameters
    parser.add_argument("-fd", "--ts_file_dir",
                        type=str,
                        default=None,
                        required=True,
                        help="path to time series files with extension .txt")

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

    parser.add_argument("-o", "--out_dir",
                        type=str,
                        default=".",
                        required=True,
                        help="directory to save output")

    args = parser.parse_args()

    ts_file_dir = args.ts_file_dir
    num_curves = args.num_curves
    num_shifts = args.num_shifts
    out_dir = args.out_dir

    taskID = int(os.environ['SLURM_ARRAY_TASK_ID']) - 1

    file_df_dict = {}
    ts_file_list = []
    for filename in os.listdir(ts_file_dir):
        if filename.endswith('.txt'):
            file_path = os.path.join(ts_file_dir, filename)
            dataname = '_'.join(filename.split('_')[:2])
            ts_file_list.append(dataname)
            file_df_dict[dataname] = pd.read_csv(file_path, delimiter='\t',
                                                               header=0, index_col=0, comment='#')

    # sample df to getting gene list and shift time
    sample_df = file_df_dict[ts_file_list[0].split('.')[0]]    

    # get gene list from first file in file list
    gene_names_set = sample_df.index.values
    gene_names_set.sort()

    # create the output directory if needed
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # initialize data frame
    diss_df = pd.DataFrame()
    shift_time_int = float(sample_df.columns.values[1]) - float(sample_df.columns.values[0])
    for shift in range(num_shifts):
        diss_df.loc['shift_times', 'shift_%s' % str(shift)] = shift * shift_time_int

    diss_df['df_1'] = ''
    diss_df['df_2'] = ''

    # create list of unique pairwise comparisons of data
    df_comparisons = list(itertools.combinations(file_df_dict.keys(), 2))
    # remove inter-parasite comparisons
    df_comparisons = [tup for tup in df_comparisons if tup[0][-3:] == tup[1][-3:]]

    # grab gene name at index equal to job index
    gene_name = gene_names_set[taskID]
    print("Job index {} : Gene {}".format(taskID, gene_name))

    # loop over pairwise combos
    for i, df_pair in enumerate(df_comparisons):
        print("Gene {}: Dataframes {}, {}".format(gene_name, df_pair[0], df_pair[1]))
        start_time = time.time()

        # identify the rna profile from ts 1
        temp_rna_ts_1 = file_df_dict[df_pair[0]].loc[gene_name].values
        # normalize data
        temp_rna_ts_1 = normalize_ts(temp_rna_ts_1, kind='scale')

        # identify the rna profile from ts 2s
        temp_rna_ts_2 = file_df_dict[df_pair[1]].loc[gene_name].values
        # normalize data
        temp_rna_ts_2 = normalize_ts(temp_rna_ts_2, kind='scale')

        m = len(temp_rna_ts_1)

        diss_gene_name = '{}__{}-{}'.format(gene_name, df_pair[0], df_pair[1])

        # compute and save dissimilarities
        # loop over shifts of rna time series 2 to test different alignments of phases
        for shift in range(num_shifts):
            diss_df.loc[diss_gene_name, 'shift_%s' % str(shift)] = kt_diss(temp_rna_ts_1[:m - shift],
                                                                      temp_rna_ts_2[shift:])

        # find the minumum and range of kt distances over the phase shifts
        diss_df.loc[diss_gene_name, 'kt_diss_range'] = max(diss_df.loc[diss_gene_name,
                                                                ['shift_%s' % str(shift)
                                                                 for shift in range(num_shifts)]]) - \
                                                min(diss_df.loc[diss_gene_name,
                                                                ['shift_%s' % str(shift)
                                                                 for shift in range(num_shifts)]])

        diss_df.loc[diss_gene_name, 'kt_diss'] = min(diss_df.loc[diss_gene_name,
                                                          ['shift_%s' % str(shift)
                                                           for shift in range(num_shifts)]])

        
        # check if num_curves if greater than the fatorial of the number of time points.
        # if True, then exhaustive search is performed
        if num_curves >= np.math.factorial(len(sample_df.columns)):
            # generate random curves from the rna profiles and compute KT measure over phase shifts
            kt_diss_cdf = gen_takt_diss_cdf_ex(temp_rna_ts_1, temp_rna_ts_2, num_shifts,
                                            dict(w=1, p=1, noise1=None, noise2=None))
        else:
            # generate random curves from the rna profiles and compute KT measure over phase shifts
            kt_diss_cdf = gen_takt_diss_cdf(temp_rna_ts_1, temp_rna_ts_2, num_curves, num_shifts,
                                            dict(w=1, p=1, noise1=None, noise2=None))
        diss_df.loc[diss_gene_name, 'kt_diss_pval'] = np.searchsorted(kt_diss_cdf, diss_df.loc[diss_gene_name, 'kt_diss']) \
                                               / num_curves

        # estimate the delay as the shift which achieves the minimum kt distance
        min_shift_idx = np.argmin(diss_df.loc[diss_gene_name, ['shift_%s' % str(shift) for shift in range(num_shifts)]].values)
        min_shift = diss_df.loc['shift_times', 'shift_%d' % min_shift_idx]
        diss_df.loc[diss_gene_name, 'est_delay'] = min_shift

        diss_df.loc[diss_gene_name, 'df_1'] = df_pair[0]
        diss_df.loc[diss_gene_name, 'df_2'] = df_pair[1]

        print('completed scoring {} between {} and {} in {} seconds'.format(gene_name, df_pair[0], df_pair[1],
                                                                            time.time() - start_time))
    
    # drop __df1-df2 from gene name in index
    new_index = diss_df.index.str.split('__')
    new_index = [gene_dfs[0] for gene_dfs in new_index]
    diss_df.index = new_index

    with open(out_dir + '/kt_diss_{}.tsv'.format(gene_name), 'w') as outfile:
        outfile.write('# RNA vs. RNA Profile Dissimilarity Measure Over Phase Shifts \n')
        outfile.write('# Time series files:  {} \n'.format(', '.join(ts_file_list)))
        if num_curves >= np.math.factorial(len(sample_df.columns)):
            outfile.write('# No. of Randomized Curves (Exhaustive):  %d \n' % np.math.factorial(len(sample_df.columns)))
        else:
            outfile.write('# No. of Randomized Curves:  %d \n' % num_curves)
        outfile.write('# Author: Robert C. Moseley \n')
        diss_df.to_csv(path_or_buf=outfile, sep='\t')
