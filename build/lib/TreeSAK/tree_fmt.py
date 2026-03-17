import os
import argparse
from ete3 import Tree


tree_fmt_usage = '''
===================== tree_fmt example commands =====================

TreeSAK tree_fmt -i in.tree -fi 0 -o subset.tree -fo 5

# https://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html
0	flexible with support values
1	flexible with internal node names
2	all branches + leaf names + internal supports
3	all branches + all names
4	leaf branches + leaf names
5	internal and leaf branches + leaf names
6	internal branches + leaf names
7	leaf branches + all names
8	all names
9	leaf names
100	topology only

=====================================================================
'''


def tree_fmt(args):

    tree_file_in    = args['i']
    input_tree_fmt  = args['fi']
    tree_file_out   = args['o']
    output_tree_fmt = args['fo']

    if os.path.isfile(tree_file_in) is False:
        print("Input tree file does not exist, program exited!")
        exit()

    input_tree = Tree(tree_file_in, quoted_node_names=False, format=input_tree_fmt)
    input_tree.write(outfile=tree_file_out, format=output_tree_fmt)
    print('Tree exported to: %s' % tree_file_out)


if __name__ == '__main__':

    tree_fmt_parser = argparse.ArgumentParser()
    tree_fmt_parser.add_argument('-i',      required=True,                       help='input tree')
    tree_fmt_parser.add_argument('-fi',     required=False, default=1, type=int, help='input tree format, default: 1')
    tree_fmt_parser.add_argument('-o',      required=True,                       help='output tree')
    tree_fmt_parser.add_argument('-fo',     required=False, default=1, type=int, help='output tree format, default: 1')
    args = vars(tree_fmt_parser.parse_args())
    tree_fmt(args)
