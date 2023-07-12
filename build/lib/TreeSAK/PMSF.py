import os
import argparse


PMSF_usage = '''
==================== PMSF example commands ====================

# Dependency: iqtree

TreeSAK PMSF -i concatenated.phy -o get_PMSF_tree_wd -t 12

# more information: http://www.iqtree.org/doc/Complex-Models

===============================================================
'''


def PMSF(args):

    msa_in                  = args['s']
    iqtree_model_guide_tree = args['gm']
    iqtree_model            = args['m']
    op_dir                  = args['o']
    tree_prefix             = args['p']
    force_overwrite         = args['f']
    num_of_threads          = args['t']

    guide_tree_wd  = '%s/guide_tree'           % op_dir
    pwd_guide_tree = '%s/guide_tree.treefile'  % guide_tree_wd
    pwd_cmd_txt    = '%s/cmds.txt'             % op_dir

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

    # get guide tree
    guidetree_cmd = 'iqtree -s %s --prefix %s/guide_tree --seqtype AA -m %s -T %s -B 1000 --alrt 1000 --quiet' % (msa_in, guide_tree_wd, iqtree_model_guide_tree, num_of_threads)
    os.system(guidetree_cmd)

    # get PMSF tree
    iqtree_cmd = 'iqtree -s %s --prefix %s/%s --seqtype AA -m %s -T %s -B 1000 --alrt 1000 --quiet -ft %s' % (msa_in, op_dir, tree_prefix, iqtree_model, num_of_threads, pwd_guide_tree)
    os.system(iqtree_cmd)

    # write out commands
    pwd_cmd_txt_handle = open(pwd_cmd_txt, 'w')
    pwd_cmd_txt_handle.write(guidetree_cmd + '\n')
    pwd_cmd_txt_handle.write(iqtree_cmd + '\n')
    pwd_cmd_txt_handle.close()
    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',   required=True,                          help='MSA file')
    parser.add_argument('-gm',  required=False, default='LG+F+G',       help='iqtree model for guide tree, default: LG+F+G')
    parser.add_argument('-m',   required=False, default='LG+C60+F+G',   help='iqtree model, default: LG+C60+F+G')
    parser.add_argument('-o',   required=True,                          help='output plot')
    parser.add_argument('-p',   required=False, default='PMSF',         help='tree prefix, default: PMSF')
    parser.add_argument('-t',   required=False, type=int, default=1,    help='num of threads')
    parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    PMSF(args)

'''

cd /home-user/wzsong/DateArTree_GTDB_r214_2
python3 PMSF.py -s gtdbtk.ar53.r214.phy -o gtdbtk.ar53.r214.PMSF_wd -p gtdbtk.ar53.r214 -t 12 

'''