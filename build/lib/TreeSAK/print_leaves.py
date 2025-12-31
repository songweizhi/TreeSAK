import argparse
from ete3 import Tree


print_leaves_usage = '''
======= print_leaves example commands =======

TreeSAK print_leaves -i in.tree

=============================================
'''


def print_leaves(args):

    tree_file_in = args['i']
    op_txt = args['o']

    leaf_list = []
    for leaf in Tree(tree_file_in, format=1):
        leaf_name = leaf.name
        leaf_list.append(leaf_name)

    with open(op_txt, 'w') as f:
        f.write('\n'.join(sorted(leaf_list)))


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,   help='input tree file')
    parser.add_argument('-o',      required=True,   help='output text file')
    args = vars(parser.parse_args())
    print_leaves(args)
