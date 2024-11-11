import os
import argparse


pruneMSA_usage = '''
==================== pruneMSA example commands ====================

# Dependencies: perl and alignment_pruner.pl

TreeSAK pruneMSA -i input_msa.fasta -c 10
TreeSAK pruneMSA -i input_msa.fasta -c 5,10,20,30,40

Note:
1. This is a wrapper for alignment_pruner.pl (--chi2_prune mode) 
2. For details: https://doi.org/10.1038/s41467-020-17408-w

===================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def pruneMSA(args):

    msa_in              = args['i']
    conserved_cutoffs   = args['c']

    _, msa_path, msa_base, msa_ext = sep_path_basename_ext(msa_in)

    current_file_path   = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    alignment_pruner_pl = '%s/alignment_pruner.pl'  % current_file_path
    cutoff_list         = conserved_cutoffs.split(',')

    op_file_list = []
    for each_cutoff in cutoff_list:
        cutoff_formatted = str(float(each_cutoff)/100)
        current_msa_out     = '%s/%s_chi2p%s.%s'                        % (msa_path, msa_base, each_cutoff, msa_ext)
        perl_cmd            = 'perl %s --file %s --chi2_prune f%s > %s' % (alignment_pruner_pl,   msa_in, cutoff_formatted, current_msa_out)
        perl_cmd_to_print   = 'perl %s --file %s --chi2_prune f%s > %s' % ('alignment_pruner.pl', msa_in, cutoff_formatted, current_msa_out)
        op_file_list.append(current_msa_out)
        print(perl_cmd_to_print)
        os.system(perl_cmd)

    # report
    print('Pruned MSA exported to:')
    print('\n'.join(op_file_list))


if __name__ == '__main__':

    pruneMSA_parser = argparse.ArgumentParser()
    pruneMSA_parser.add_argument('-i', required=True, help='input MSA file')
    pruneMSA_parser.add_argument('-c', required=True, help='conservation cutoffs, comma separated')
    args = vars(pruneMSA_parser.parse_args())
    pruneMSA(args)
