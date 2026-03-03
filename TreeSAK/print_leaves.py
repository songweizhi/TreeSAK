import argparse
from ete3 import Tree


print_leaves_usage = '''
=========== print_leaves example commands ===========

TreeSAK print_leaves -i in.tree -o leaves.txt
TreeSAK print_leaves -i in.tree -o leaves.txt -fmt 2

=====================================================
'''


def print_leaves(args):

    tree_file_in = args['i']
    op_txt       = args['o']
    tree_fmt     = args['fmt']

    leaf_list = []
    for leaf in Tree(tree_file_in, format=tree_fmt):
        leaf_name = leaf.name
        leaf_list.append(leaf_name)

    with open(op_txt, 'w') as f:
        f.write('\n'.join(sorted(leaf_list)))


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,                           help='input tree file')
    parser.add_argument('-fmt',    required=False, default=1, type=int,     help='tree format, default is 1')
    parser.add_argument('-o',      required=True,                           help='output file')
    args = vars(parser.parse_args())
    print_leaves(args)
