import os
import argparse
from distutils.spawn import find_executable


GeneTree_usage = '''
================ GeneTree example commands ================

TreeSAK GeneTree -i amoA.faa -o amoA_tree -t 36 -f -bmge

===========================================================
'''


def check_dependencies(program_list):

    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not found, program exited!' % ','.join(not_detected_programs))
        exit()


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def GeneTree(args):

    seq_file                    = args['i']
    num_threads                 = args['t']
    op_dir                      = args['o']
    force_create_op_dir         = args['f']
    trim_with_bmge              = args['bmge']
    bmge_trim_model             = args['bmge_m']
    bmge_entropy_score_cutoff   = args['bmge_esc']

    # check dependencies
    check_dependencies(['java', 'mafft-einsi'])

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar      = '%s/BMGE.jar' % current_file_path

    # determine the version of iqtree available on the system
    if find_executable('iqtree2'):
        iqtree_exe = 'iqtree2'
    elif find_executable('iqtree'):
        iqtree_exe = 'iqtree'
    else:
        print('iqtree not detected, program exited!')
        exit()

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    ######################################## define output file name ########################################

    sep_file_path, sep_file_base, sep_file_ext = sep_path_basename_ext(seq_file)
    get_gene_tree_cmds_txt  = '%s/cmds.txt'     % op_dir
    msa_file                = '%s/%s.aln'       % (op_dir, sep_file_base)
    msa_file_bmge           = '%s/%s.bmge.aln'  % (op_dir, sep_file_base)

    msa_for_iqtree = msa_file
    if trim_with_bmge is True:
        msa_for_iqtree = msa_file_bmge

    #########################################################################################################

    # prepare commands
    mafft_cmd   = 'mafft-einsi --thread %s --quiet %s > %s'                 % (num_threads, seq_file, msa_file)
    trim_cmd    = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'             % (pwd_bmge_jar, msa_file, bmge_trim_model, bmge_entropy_score_cutoff, msa_file_bmge)
    iqtree_cmd  = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s/%s'    % (iqtree_exe, num_threads, msa_for_iqtree, op_dir, sep_file_base)

    # write out commands
    with open(get_gene_tree_cmds_txt, 'w') as f:
        if trim_with_bmge is True:
            f.write('%s\n%s\n%s\n' % (mafft_cmd, trim_cmd, iqtree_cmd))
        else:
            f.write('%s\n%s\n' % (mafft_cmd, iqtree_cmd))

    # run mafft
    print(mafft_cmd)
    os.system(mafft_cmd)

    # run BMGE
    print(trim_cmd)
    os.system(trim_cmd)

    # run iqtree
    print(iqtree_cmd)
    os.system(iqtree_cmd)

    print('Done!')


if __name__ == '__main__':

    GeneTree_parser = argparse.ArgumentParser()
    GeneTree_parser.add_argument('-i',          required=False, default=None,           help='sequence file')
    GeneTree_parser.add_argument('-o',          required=True,                          help='output dir, i.e., OMA working directory')
    GeneTree_parser.add_argument('-t',          required=False, type=int, default=3,    help='number of threads specified in job script, default: 3')
    GeneTree_parser.add_argument('-f',          required=False, action="store_true",    help='force overwrite')
    GeneTree_parser.add_argument('-bmge',       required=False, action="store_true",    help='trim MSA with BMGE, default no trimming')
    GeneTree_parser.add_argument('-bmge_m',     required=False, default='BLOSUM30',     help='BMGE trim model, default: BLOSUM30')
    GeneTree_parser.add_argument('-bmge_esc',   required=False, default='0.55',         help='BMGE entropy score cutoff, default: 0.55')
    args = vars(GeneTree_parser.parse_args())
    GeneTree(args)
