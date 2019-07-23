import numpy as np

def get_gene_name(gene_id):
    return gene_id.split('_')[0]


def get_pep_id(gene_id):
    return gene_id.split('_')[1]


def get_rep_id(gene_id):
    return gene_id.split('_')[3][1]


def normalize_ts(ts_array, kind='scale'):
    if kind is 'scale':
        ts_array = ts_array - np.mean(ts_array)
        ts_array = ts_array / np.std(ts_array)

    elif kind is 'max':
        ts_array = ts_array / max(ts_array)

    elif kind is '01':
        ts_array = ts_array - min(ts_array)
        ts_array = ts_array / max(ts_array)

    return ts_array


def kt_diss(L1, L2, w=1, p=1, noise1=None, noise2=None):
    """
    Function for computing the Kendall distance with penalty parameter p as defined in
    http://epubs.siam.org/doi/abs/10.1137/05063088X

    Optional argument/value pairs:
        - w: 0 <= w <= 1 specifies the penalty for discordance due to differences in rank order between the lists
             (type 1) discordance.  (default = 1)
        - p: 0 <= p <= 1 specifies the penalty for discordance due to one list having ambiguous ranking (type 2).
             This function is a metric only if 1/2 < p <= 1, is a near metric if 0 < p <= 1/2, and is a
             "near pseudometric" when p = 0. (default = 1)
        - noise1/noise2: lists of non-negative numbers the same length as L1/L2.  noise1[i] is the uncertainty of
                         L1[i]. (default set to 0 for each i, i.e. no noise)
    """
    m = len(L1)

    # Check that lists have same length
    if m != len(L2):
        length_mismatch()
        return

    # Initialize variables and set defaults in parameter dictionary if not defined by user
    disc_pair_count = 0

    if noise1 is None:
        e1 = np.zeros((m, 1))  # Default noise model for L1
    else:
        e1 = noise1

    if noise2 is None:
        e2 = np.zeros((m, 1))  # Default noise model for L2
    else:
        e2 = noise2

    # Check that noise models contain allowable values and are valid length
    if m != len(e1) or m != len(e2):
        length_mismatch()
        return
    elif min(e1) < 0 or min(e2) < 0:
        print('ERROR: Error models must contain non-negative error bounds.')
        return

    # Sum penalty functions \bar{K}_{i,j}(L1, L2) over all pairs (i, j), i < j
    for i in range(len(L1)):
        for j in range(i + 1, len(L1)):
            if L1[i] - e1[i] >= L1[j] + e1[j] and L2[i] + e2[i] <= L2[j] - e2[j]:
                disc_pair_count += w

            elif L1[i] + e1[i] <= L1[j] - e1[j] and L2[i] - e2[i] >= L2[j] + e2[j]:
                disc_pair_count += w

            elif abs(L1[i] - L1[j]) >= e1[i] + e1[j] and abs(L2[i] - L2[j]) < e2[i] + e2[j]:
                disc_pair_count += p

            elif abs(L1[i] - L1[j]) < e1[i] + e1[j] and abs(L2[i] - L2[j]) >= e2[i] + e2[j]:
                disc_pair_count += p

    return float(disc_pair_count) / float(m*(m-1)/2)


def gen_takt_diss_cdf(rna_curve, prot_curve, num_curves, num_shifts, args_dict):
    # generate num_curves random rna curves from the given rna profile and compute their dissimilarity scores from the
    # corresponding protein profile
    m = len(prot_curve)
    diss_cdf = np.empty(num_curves, dtype=float)
    # generate num_curves permutations of rna_curve
    for i in range(num_curves):
        # compute kt_diss over temporal shifts 0-num_shifts
        temp_shift = np.empty(shape=num_shifts)
        for shift in range(num_shifts):
            temp_shift[shift]= kt_diss(prot_curve[shift:], np.random.permutation(rna_curve[:m - shift]), **args_dict)
        # compute the takt dissimilarity measure
        diss_cdf[i] = np.min(temp_shift)  # minimum dissimilarity over shifts

    diss_cdf.sort()

    return diss_cdf


def length_mismatch():
    print('ERROR: The lengths of required variables must match.')