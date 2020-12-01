import os
import pandas as pd


def sep_patient_dfs(full_df, out_dir, time_points):
    parasite = list(set([col.split('_')[1] for col in full_df.columns]))[0]
    patients = list(set([col.split('_')[0][:-2] for col in full_df.columns]))
    df_dict = {}
    print(parasite)
    for pat in patients:
        pat_df = full_df.filter(regex=pat, axis=1)
        int_cols = [int(col.split('_')[0][-2:]) for col in pat_df.columns]
        int_cols.sort()
        pat_df.columns = int_cols
        pat_df.columns = time_points
        df_dict[pat] = pat_df
        pat_df.to_csv(os.path.join(out_dir, '{}_{}_ts.txt'.format(pat, parasite)), sep='\t')
        print('Patient {}; {} genes'.format(pat, pat_df.shape[0]))
    return df_dict


def ring_stage_dfs(ring_stage_tps, patient_df_dir, output_dir):
    for pat_file in os.listdir(patient_df_dir):
        if pat_file.endswith('.txt'):
            pat_df = pd.read_csv(os.path.join(patient_df_dir, pat_file), sep='\t', index_col=0)
            if 'FUP' in pat_file:
                pat = 'FUP'
            else:
                pat = '3D7'
            # columns truncated on ring stage
            ring_df = pat_df.loc[:, :str(ring_stage_tps[pat])]
            ring_df.to_csv(os.path.join(output_dir, pat_file.split('.')[0] + '_ring.txt'), sep='\t')


def make_slurm_jobs_file(batch_args_dict, input_dir, output_dir, takt_dir, perms, shifts, sp):
    with open(f'takt_{sp}_ring_rcm.sh', 'w') as job_file:
        job_file.write('#!/bin/bash\n')
        for flag, arg in batch_args_dict.items():
            if '--' in flag:
                job_file.write(f'#SBATCH {flag}={arg}\n')
            else:
                job_file.write(f'#SBATCH {flag} {arg}\n')
        job_file.write('\nmodule load Python/3.8.1\n\n')
        job_file.write(f"ring_file_dir='{input_dir}'\n")
        job_file.write(f"ring_out_dir='{output_dir}'\n\n")
        job_file.write('python -u %s -fd ${ring_file_dir} -c %d -s %d -o ${ring_out_dir}\n' % (takt_dir,
                                                                                               perms,
                                                                                               shifts))


def setup_df(file_path):
    def adj_shift(row):
        est_delay = float(row[-1])
        if est_delay > 24:
            return est_delay - 48
        else:
            return est_delay

    def sep_pat_para_cols(df):
        df[['pat1', 'para1']] = df['df_1'].str.split(expand=True, pat='_')
        df[['pat2', 'para2']] = df['df_2'].str.split(expand=True, pat='_')
        return df

    df = pd.read_csv(file_path, skiprows=4, index_col=0, sep='\t')
    df.drop(['shift_times'], inplace=True)
    df['adj_shift'] = df.apply(adj_shift, axis=1)
    df = sep_pat_para_cols(df)

    return df


def single_result_analysis(series, shift, no_shift, para_name):
    if no_shift:
        shift_move = '<'
        if -shift <= series['adj_shift'].values[0] <= shift:
            return 1, shift_move
        else:
            return 0, shift_move

    elif not no_shift:
        shift_move = '>'
        if series['adj_shift'].values[0] < -shift or series['adj_shift'].values[0] > shift:
            return 1, shift_move
        else:
            return 0, shift_move


def check_shift(df, shift, no_shift):
    shift_type = 'empty'
    if no_shift:
        for s in list(set(df['adj_shift'].tolist())):
            if -shift <= s <= shift:
                shift_type = 'none'
    elif not no_shift:
        for s in list(set(df['adj_shift'].tolist())):
            if s < -shift or s > shift:
                shift_type = 'shift'
    return shift_type


def single_para_analysis(para_df, para_name, shift, no_shift):
    if para_df.shape[0] == 1:
        count, _ = single_result_analysis(para_df, shift, no_shift, para_name)
        return count
    else:
        shift_check = check_shift(para_df, shift, no_shift)

        if no_shift and shift_check == 'none':
            para_shift_filtered_df = para_df[(para_df['adj_shift'] >= -shift) & (para_df['adj_shift'] <= shift)]
            return para_shift_filtered_df.shape[0]

        elif not no_shift and shift_check == 'shift':
            para_shift_filtered_df = para_df[(para_df['adj_shift'] < -shift) | (para_df['adj_shift'] > shift)]
            return para_shift_filtered_df.shape[0]

        elif shift_check == 'empty':
            return 0


def single_gene_analysis(df, type1, type2, shift, pval, no_shift=True):
    gene = df.index.tolist()[0]
    pat_filtered_df = df[((df['pat1'].str.contains(type1)) & (df['pat2'].str.contains(type2))) | (
                (df['pat1'].str.contains(type2)) & (df['pat2'].str.contains(type1)))]
    pval_filtered_df = pat_filtered_df[pat_filtered_df['kt_diss_pval'] <= pval]

    if pval_filtered_df.shape[0] == 0:
        return {'FUP': 0, '3D7': 0}
    else:
        total_fup = pval_filtered_df[
            (pval_filtered_df['para1'].str.contains('FUP')) & (pval_filtered_df['para2'].str.contains('FUP'))]
        total_threed7 = pval_filtered_df[
            (pval_filtered_df['para1'].str.contains('3D7')) & (pval_filtered_df['para2'].str.contains('3D7'))]

        fup_count = single_para_analysis(total_fup, 'FUP', shift, no_shift)
        threed7_count = single_para_analysis(total_threed7, '3D7', shift, no_shift)
        return {'FUP': fup_count, '3D7': threed7_count}

#             if verbose:
#                 print('{} with shifts {} {} hrs between {} and {}:\tFUP: {} 3D7: {} (Possible FUP and 3D7: {} {}) (p-val:{})'.format(
#                     gene, shift_move, shift, type1, type2, fup_df.shape[0], threed7_df.shape[0], total_fup, total_threed7, pval))
#             return {'FUP': fup_df.shape[0], '3D7': threed7_df.shape[0]}
