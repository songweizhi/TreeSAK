import os
import argparse
from distutils.spawn import find_executable


PMSF_usage = '''
==================== PMSF example commands ====================

# Dependency: iqtree2

TreeSAK PMSF -i in.aln -o get_PMSF_tree_wd -t 12

# This is a wrapper for:
iqtree2 -T 12 -B 1000 --alrt 1000 --quiet --seqtype AA -s in.aln --prefix guide_tree -m LG+F+G 
iqtree2 -T 12 -B 1000 --alrt 1000 --quiet --seqtype AA -s in.aln --prefix PMSF -m LG+C60+F+G -ft guide_tree.treefile

# more information: http://www.iqtree.org/doc/Complex-Models

# Reference: The evolutionary origin of host association in the Rickettsiales
Maximum likelihood phylogenetic reconstructions were done under the PMSF approximation 
(with 100 non-parametric bootstraps; guidetree under LG+G+F) of the LG+C60+F+Γ4 model 
(selected by ModelFinder) for both supermatrix alignments with IQTREE v1.6.5.

===============================================================
'''


def PMSF(args):

    msa_in                  = args['i']
    iqtree_model_guide_tree = args['gm']
    iqtree_model            = args['m']
    op_dir                  = args['o']
    tree_prefix             = args['p']
    force_overwrite         = args['f']
    num_of_threads          = args['t']

    guide_tree_wd  = '%s/guide_tree'           % op_dir
    pwd_guide_tree = '%s/guide_tree.treefile'  % guide_tree_wd
    pwd_cmd_txt    = '%s/cmds.txt'             % op_dir

    iqtree_exe = ''
    if find_executable('iqtree2'):
        iqtree_exe = 'iqtree2'
    elif find_executable('iqtree'):
        iqtree_exe = 'iqtree'
    else:
        print('iqtree not detected, program exited!')
        exit()

    # check input file
    if os.path.isfile(msa_in) is False:
        print('MSA file not found, program exited!')
        exit()

    # create output dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % guide_tree_wd)

    guidetree_cmd = '%s -s %s --prefix %s/guide_tree --seqtype AA -m %s -T %s -B 1000 --alrt 1000 --quiet' % (iqtree_exe, msa_in, guide_tree_wd, iqtree_model_guide_tree, num_of_threads)
    iqtree_cmd    = '%s -s %s --prefix %s/%s --seqtype AA -m %s -T %s -B 1000 --alrt 1000 --quiet -ft %s'  % (iqtree_exe, msa_in, op_dir, tree_prefix, iqtree_model, num_of_threads, pwd_guide_tree)

    # write out commands
    pwd_cmd_txt_handle = open(pwd_cmd_txt, 'w')
    pwd_cmd_txt_handle.write(guidetree_cmd + '\n')
    pwd_cmd_txt_handle.write(iqtree_cmd + '\n')
    pwd_cmd_txt_handle.close()

    # get guide tree
    print('Building guide tree')
    print(guidetree_cmd)
    os.system(guidetree_cmd)

    # get PMSF tree
    print('Building PMSF tree with model %s' % iqtree_model)
    print(iqtree_cmd)
    os.system(iqtree_cmd)

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=True,                          help='input MSA file')
    parser.add_argument('-gm',  required=False, default='LG+F+G',       help='iqtree model for guide tree, default: LG+F+G')
    parser.add_argument('-m',   required=False, default='LG+C60+F+G',   help='iqtree model, default: LG+C60+F+G')
    parser.add_argument('-o',   required=True,                          help='output plot')
    parser.add_argument('-p',   required=False, default='PMSF',         help='tree prefix, default: PMSF')
    parser.add_argument('-t',   required=False, type=int, default=1,    help='num of threads')
    parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    PMSF(args)

