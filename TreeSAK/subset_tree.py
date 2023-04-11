import argparse
from ete3 import Tree


subset_tree_usage = '''
================ subset_tree example commands ================

TreeSAK subset_tree -i in.tree -k leaves.txt -o subset.tree

==============================================================
'''


def subset_tree(args):

    tree_file_in    = args['i']
    to_keep_txt     = args['k']
    tree_file_out   = args['o']

    genomes_to_keep = [i.strip() for i in open(to_keep_txt)]
    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(genomes_to_keep, preserve_branch_length=True)
    subset_tree.write(outfile=tree_file_out)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,   help='input tree file')
    parser.add_argument('-k',      required=True,   help='leaves to keep, one id per line')
    parser.add_argument('-o',      required=True,   help='output tree file')
    args = vars(parser.parse_args())
    subset_tree(args)
