import os
import time
import pandas as pd
from utilities import *

"""
This script computes a measure of dynamic dissimilarity between an rna expression profile and a corresponding protein 
expression profile using the normalized Kendall-Tau (KT) distance function between ranked lists over a range of phase 
shifts. The distance between the two profiles is taken to be the minimum KT distance achieved by some phase shift. 
The fraction of randomized expression profiles which achieve a smaller minimum KT distance over the same set of phase 
shifts is also reported to estimate the significance of the computed measure of dissimilarity. 
"""
# user-defined parameters
out_dir = 'data/prot_rna_diss_output/'  # where to say the output
ts_file = 'data/clb1-6_linearSpline_time-series.tsv'  # the input data file
num_curves = 10000  # the number of permutations of each curve to perform to derive an empirical p-value
num_shifts = 11  # number of time point shifts (including 0) to consider (match length of one cell-cycle period)
# ---------------------------------------------------------------

filename = os.path.splitext(os.path.basename(ts_file))[0]
ts_df = pd.read_csv(ts_file, delimiter='\t', header=0, index_col=0, comment='#')

# identify a list of unique gene names
gene_names = []
for gene_id in ts_df.index.values:
    gene_names.append(get_gene_name(gene_id))
gene_names_set = set(gene_names)

# create the output directory if needed
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# initialize data frame
diss_df = pd.DataFrame()
shift_time_int = float(ts_df.columns.values[1]) - float(ts_df.columns.values[0])
for shift in range(11):
    diss_df.loc['shift_times', 'shift_%s' % str(shift)] = shift * shift_time_int

# loop over each gene name
for gene_name in gene_names_set:
    # extract all records corresponding to the current gene
    gene_indices = [gene_name == gene for gene in gene_names]
    gene_ts_df = ts_df.iloc[gene_indices]

    temp_rna_ts = None
    # identify the rna profile
    for gene_id in gene_ts_df.index.values:
        if gene_id is get_gene_name(gene_id):
            temp_rna_ts = ts_df.loc[gene_id].values
            # normalize data
            temp_rna_ts = normalize_ts(temp_rna_ts, kind='scale')

    # loop over each protein profile and compare to rna profile
    for gene_id in gene_ts_df.index.values:
        start_time = time.time()
        # only try to make comparisons if rna exists
        if temp_rna_ts is not None:
            # loop over protein reps/pep profiles
            if gene_id is not get_gene_name(gene_id):
                temp_prot_ts = ts_df.loc[gene_id].values
                # normalize data
                temp_prot_ts = normalize_ts(temp_prot_ts, kind='scale')
                m = len(temp_prot_ts)

                # compute and save dissimilarities
                # loop over shifts of protein time series to test different alignments of phases
                for shift in range(11):
                    diss_df.loc[gene_id, 'shift_%s' % str(shift)] = kt_diss(temp_rna_ts[:m - shift],
                                                                            temp_prot_ts[shift:])

                # find the minumum and range of kt distances over the phase shifts
                diss_df.loc[gene_id, 'kt_diss_range'] = max(diss_df.loc[gene_id,
                                                                        ['shift_%s' % str(shift)
                                                                         for shift in range(11)]]) - \
                                                        min(diss_df.loc[gene_id,
                                                                        ['shift_%s' % str(shift)
                                                                         for shift in range(11)]])

                diss_df.loc[gene_id, 'kt_diss'] = min(diss_df.loc[gene_id,
                                                                  ['shift_%s' % str(shift)
                                                                   for shift in range(11)]])

                # generate random curves from the rna and prot profiles and compute KT measure over phase shifts
                kt_diss_cdf = gen_takt_diss_cdf(temp_rna_ts, temp_prot_ts, num_curves, num_shifts,
                                                dict(w=1, p=1, noise1=None, noise2=None))
                diss_df.loc[gene_id, 'kt_diss_pval'] = np.searchsorted(kt_diss_cdf, diss_df.loc[gene_id, 'kt_diss']) \
                                                       / num_curves

                # estimate the protein delay as the shift which achieves the minimum kt distance
                min_shift_idx = np.argmin(diss_df.loc[gene_id, ['shift_%s' % str(shift) for shift in range(11)]].values)
                min_shift = diss_df.loc['shift_times', 'shift_%d' % min_shift_idx]
                diss_df.loc[gene_id, 'est_delay'] = min_shift
                diss_df.loc[gene_id, 'min_dist_flag'] = 0
                diss_df.loc[gene_id, 'min_pval_flag'] = 0

                print('completed scoring %s in %f seconds' % (gene_id, time.time() - start_time))

    # extract all dissimilarity records corresponding to the current gene
    gene_diss_df = diss_df.iloc[[gene_name == get_gene_name(gene_id) for gene_id in diss_df.index.values]]

    # only attempt to update gene records for which comparisons between prot and rna were made
    # (in particular because ASE1 for wild-type data only has RNA data)
    if len(gene_diss_df.index.values) > 0:
        # find the protein which achieves the minimum kt dissimilarity and p-value
        min_kt_diss_gene_id = gene_diss_df.iloc[np.argmin(gene_diss_df['kt_diss'].values)].name
        min_kt_pval_gene_id = gene_diss_df.iloc[np.argmin(gene_diss_df['kt_diss_pval'].values)].name

        # flag the protein profile which achieves the minimum kt distance (and min. p-value)
        diss_df.loc[min_kt_diss_gene_id, 'min_dist_flag'] = 1
        diss_df.loc[min_kt_pval_gene_id, 'min_pval_flag'] = 1

with open(out_dir + '/kt_diss_%s.tsv' % filename, 'w') as outfile:
    outfile.write('# Prot. vs. RNA Profile Dissimilarity Measure Over Phase Shifts \n')
    outfile.write('# Time series file:  %s \n' % os.path.abspath(ts_file))
    outfile.write('# No. of Randomized Curves:  %d \n' % num_curves)
    diss_df.to_csv(path_or_buf=outfile, sep='\t')
