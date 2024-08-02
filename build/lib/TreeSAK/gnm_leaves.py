import os
import argparse
from ete3 import Tree


gnm_leaves_usage = '''
========== gnm_leaves example commands ==========

TreeSAK gnm_leaves -i input.tree -o output.tree

=================================================
'''


def gnm_leaves(args):

    tree_file_in  = args['i']
    tree_file_out = args['o']
    tree_format   = args['fmt']

    if os.path.isfile(tree_file_in) is False:
        print('Tree file not found, program exited!')
        exit()

    t = Tree(tree_file_in, format=tree_format)

    for leaf in t:
        leaf_name = leaf.name
        leaf_name_new = '_'.join(leaf_name.split('_')[:-1])
        leaf.name = leaf_name_new
    t.write(format=tree_format, outfile=tree_file_out)

    print('Done!')


if __name__ == '__main__':

    gnm_leaves_parser = argparse.ArgumentParser()
    gnm_leaves_parser.add_argument('-i',    required=True,                       help='input tree')
    gnm_leaves_parser.add_argument('-o',    required=True,                       help='output tree')
    gnm_leaves_parser.add_argument('-fmt',  required=False, default=1, type=int, help='tree format, default: 1')
    args = vars(gnm_leaves_parser.parse_args())
    gnm_leaves(args)
