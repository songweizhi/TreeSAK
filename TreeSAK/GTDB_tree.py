import os
import argparse


GTDB_tree_usage = '''
======================== GTDB_tree example command ========================

export GTDBTK_DATA_PATH=/scratch/PI/boqianpy/Database/gtdb_r220/release220
TreeSAK GTDB_tree -p Demo -i gnm_folder -x fa -t 12

# This is a wrapper for the following commands
gtdbtk identify --genome_dir gnm_folder -x fa --out_dir op_dir --cpus 12
gtdbtk align --identify_dir Demo_op_dir --out_dir op_dir --cpus 12
gtdbtk infer --msa_file Demo_op_dir/align/gtdbtk.bac120.user_msa.fasta.gz --out_dir op_dir --cpus 12 --prefix Demo_bac120
gtdbtk infer --msa_file Demo_op_dir/align/gtdbtk.ar53.user_msa.fasta.gz --out_dir op_dir --cpus 12 --prefix Demo_ar53

===========================================================================
'''

def GTDB_tree(args):

    input_gnm_dir       = args['i']
    output_prefix       = args['p']
    file_extension      = args['x']
    num_threads         = args['t']

    output_dir          = '%s_GTDB_tree'                                                            % output_prefix
    msa_bac120_gz       = '%s/align/gtdbtk.bac120.user_msa.fasta.gz'                                % output_dir
    msa_bac120          = '%s/align/gtdbtk.bac120.user_msa.fasta'                                   % output_dir
    msa_ar53_gz         = '%s/align/gtdbtk.ar53.user_msa.fasta.gz'                                  % output_dir
    msa_ar53            = '%s/align/gtdbtk.ar53.user_msa.fasta'                                     % output_dir

    cmd_identify        = 'gtdbtk identify --genome_dir %s -x %s --out_dir %s --cpus %s'            % (input_gnm_dir, file_extension, output_dir, num_threads)
    cmd_align           = 'gtdbtk align --identify_dir %s --out_dir %s --cpus %s'                   % (output_dir, output_dir, num_threads)
    cmd_gunzip_bac120   = 'gunzip %s'                                                               % msa_bac120_gz
    cmd_gunzip_ar53     = 'gunzip %s'                                                               % msa_ar53_gz
    cmd_infer_bac120    = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix %s_bac120'    % (msa_bac120, output_dir, num_threads, output_prefix)
    cmd_infer_ar53      = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix %s_ar53'      % (msa_ar53, output_dir, num_threads, output_prefix)

    print(cmd_identify)
    os.system(cmd_identify)
    print(cmd_align)
    os.system(cmd_align)

    if os.path.isfile(msa_bac120_gz):
        print(cmd_gunzip_bac120)
        os.system(cmd_gunzip_bac120)
        print(cmd_infer_bac120)
        os.system(cmd_infer_bac120)

    if os.path.isfile(msa_ar53_gz):
        print(cmd_gunzip_ar53)
        os.system(cmd_gunzip_ar53)
        print(cmd_infer_ar53)
        os.system(cmd_infer_ar53)

    inferred_bac120_tree = '%s/%s_bac120.unrooted.tree' % (output_dir, output_prefix)
    inferred_ar53_tree   = '%s/%s_ar53.unrooted.tree'   % (output_dir, output_prefix)

    if os.path.isfile(inferred_bac120_tree):
        print('Inferred bacterial tree:\t%s' % inferred_bac120_tree)
    if os.path.isfile(inferred_ar53_tree):
        print('Inferred archaeal tree:\t%s' % inferred_ar53_tree)

    print('Done!')


if __name__ == '__main__':

    GTDB_tree_parser = argparse.ArgumentParser(usage=GTDB_tree_usage)
    GTDB_tree_parser.add_argument('-p',               required=True,                            help='output prefix')
    GTDB_tree_parser.add_argument('-i',               required=True,                            help='genome folder')
    GTDB_tree_parser.add_argument('-x',               required=True,                            help='genome file extension')
    GTDB_tree_parser.add_argument('-t',               required=False, type=int, default=1,      help='number of threads')
    args = vars(GTDB_tree_parser.parse_args())
    GTDB_tree(args)
