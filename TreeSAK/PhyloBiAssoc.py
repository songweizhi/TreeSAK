import os
import argparse
import pandas as pd
import multiprocessing as mp
from statsmodels.stats.multitest import multipletests


PhyloBiAssoc_usage = '''
============================= PhyloBiAssoc example commands =============================

BioSAK PhyloBiAssoc -i demo.tre -d demo.txt -o op_dir -t 10 -f

# Note, header for the first two columns in -d has to be "ID" and "cate"!!!

# It will perform:
# 1) binaryPGLMM test if phylosig p-value <= 0.05 (significant phylogenetic signal)
# 2) chi-squared test if phylosig p-value > 0.05  (no phylogenetic signal)
# 3) do nothing if phylosig returns NaN (might due to the same value across all genomes)

# https://www.rdocumentation.org/packages/ape/versions/5.7-1/topics/binaryPGLMM

=========================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def subset_df(file_in, rows_to_keep, cols_to_keep, sep_symbol, row_name_pos, column_name_pos, file_out):

    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)

    if len(rows_to_keep) == 0:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep]
    else:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[rows_to_keep, :]
        else:
            subset_df = df.loc[rows_to_keep, cols_to_keep]

    subset_df.to_csv(file_out, sep=sep_symbol)


def PhyloBiAssoc(args):

    tree_file           = args['i']
    data_file           = args['d']
    op_dir              = args['o']
    num_threads         = args['t']
    force_create_op_dir = args['f']

    pwd_current_script  = os.path.realpath(__file__)
    current_script_path = '/'.join(pwd_current_script.split('/')[:-1])
    PhyloBiAssoc_R      = '%s/PhyloBiAssoc.R' % current_script_path

    cmd_txt                     = '%s/cmds.txt'                         % op_dir
    df_subset_dir               = '%s/df_subset'                        % op_dir
    stats_op_dir                = '%s/stats_results'                    % op_dir
    combined_stats_txt          = '%s/stats_results_all.txt'            % op_dir
    combined_stats_txt_sig      = '%s/stats_results_0.05.txt'           % op_dir
    combined_stats_txt_adjusted = '%s/stats_results_0.05_adjusted.txt'  % op_dir

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output directory exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % df_subset_dir)
    os.system('mkdir %s' % stats_op_dir)

    # read in dataframe
    df = pd.read_csv(data_file, sep='\t', header=0, index_col=0)
    col_header_list = list(df.columns.values)

    subset_dict = dict()
    for each_col in col_header_list[1:]:
        subset_dict[each_col] = ['cate', each_col]

    # subset dataframe
    cmd_txt_handle = open(cmd_txt, 'w')
    stats_cmd_list = []
    op_stats_txt_set = set()
    for each_subset in subset_dict:
        cols_to_keep   = subset_dict[each_subset]
        df_subset_file = '%s/%s.tab' % (df_subset_dir, each_subset)
        stats_out_txt  = '%s/%s.txt' % (stats_op_dir, each_subset)
        subset_df(data_file, set(), cols_to_keep, '\t', 0, 0, df_subset_file)
        stats_cmd = 'Rscript %s -t %s -d %s > %s' % (PhyloBiAssoc_R, tree_file, df_subset_file, stats_out_txt)
        cmd_txt_handle.write(stats_cmd + '\n')
        stats_cmd_list.append(stats_cmd)
        op_stats_txt_set.add(stats_out_txt)
    cmd_txt_handle.close()

    print('Processing %s objects with %s cores' % (len(stats_cmd_list), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, stats_cmd_list)
    pool.close()
    pool.join()
    #os.system('cp /Users/songweizhi/Documents/Research/Sponge/12_PhyloBiAssoc_wd/PhyloBiAssoc_wd_backup/stats_results/*.txt %s/' % stats_op_dir)

    # combine stats results
    sig_list_id = []
    sig_list_value = []

    combined_stats_txt_handle = open(combined_stats_txt, 'w')
    combined_stats_txt_handle.write('ID	phylosig	binaryPGLMM	chisq.test	coefficient	significance\n')
    combined_stats_txt_sig_handle = open(combined_stats_txt_sig, 'w')
    for each_file in sorted(list(op_stats_txt_set)):
        f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        for each_line in open(each_file):
            if not each_line.startswith('ID\tphylosig\tbinaryPGLMM\tchisq.test'):
                each_line_split = each_line.strip().split('\t')
                significance = each_line_split[5]
                combined_stats_txt_handle.write('%s\t%s\n' % (f_base, '\t'.join(each_line_split[1:])))
                if significance == 'y':
                    combined_stats_txt_sig_handle.write(f_base + '\n')
                    sig_bi = each_line_split[2]
                    sig_chi = each_line_split[3]
                    current_sig = ''
                    if sig_bi == 'na':
                        current_sig = sig_chi
                    elif sig_chi == 'na':
                        current_sig = sig_bi
                    sig_list_id.append(f_base)
                    sig_list_value.append(float(current_sig))
    combined_stats_txt_handle.close()
    combined_stats_txt_sig_handle.close()

    # perform Bonferroni correction
    sig_list_value_adjusted = list(multipletests(sig_list_value, alpha=0.1, method='bonferroni')[1])

    # write out adjusted p values
    combined_stats_txt_adjusted_handle = open(combined_stats_txt_adjusted, 'w')
    combined_stats_txt_adjusted_handle.write('ID\tadjusted_p_value\n')
    for (id, adjusted_p) in zip(sig_list_id, sig_list_value_adjusted):
        if adjusted_p <= 0.05:
            combined_stats_txt_adjusted_handle.write('%s\t%s\n' % (id, adjusted_p))
    combined_stats_txt_adjusted_handle.close()

    # Final report
    print('Results exported to: \n%s\n%s' % (combined_stats_txt, combined_stats_txt_adjusted))
    print('Done')


if __name__ == "__main__":

    PhyloBiAssoc_parser = argparse.ArgumentParser(usage=PhyloBiAssoc_usage)
    PhyloBiAssoc_parser.add_argument('-i', required=True,                       help='tree file')
    PhyloBiAssoc_parser.add_argument('-d', required=True,                       help='data file')
    PhyloBiAssoc_parser.add_argument('-o', required=True,                       help='output directory')
    PhyloBiAssoc_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads, default: 1')
    PhyloBiAssoc_parser.add_argument('-f', required=False, action="store_true", help='force overwrite')
    args = vars(PhyloBiAssoc_parser.parse_args())
    PhyloBiAssoc(args)
