import argparse
from ete3 import Tree


subset_usage = '''
====================== subset example commands ======================

TreeSAK subset_tree -fmt 1 -i in.tree -k leaves.txt -o subset.tree
TreeSAK subset_tree -fmt 1 -i in.tree -r leaves.txt -o subset.tree

=====================================================================
'''


def subset(args):

    tree_file_in    = args['i']
    tree_file_out   = args['o']
    to_keep_txt     = args['k']
    to_remove_txt   = args['r']
    tree_fmt        = args['fmt']


    genomes_to_keep = []
    if (to_keep_txt is None) and (to_remove_txt is None):
        print('Please specify either -k or -r, program exited!')
        exit()

    elif (to_keep_txt is not None) and (to_remove_txt is not None):
        print('Please do NOT specify -k and -r at the same time, program exited!')
        exit()

    elif (to_keep_txt is not None) and (to_remove_txt is None):
        genomes_to_keep = [i.strip() for i in open(to_keep_txt)]

    elif (to_keep_txt is None) and (to_remove_txt is not None):

        genomes_to_remove = [i.strip() for i in open(to_remove_txt)]

        leaf_list = []
        for leaf in Tree(tree_file_in, quoted_node_names=True, format=tree_fmt):
            leaf_name = leaf.name
            leaf_list.append(leaf_name)

        for each_leaf in leaf_list:
            if each_leaf not in genomes_to_remove:
                genomes_to_keep.append(each_leaf)

        if len(leaf_list) == len(genomes_to_keep):
            print('No leaf to remove, program exited!')
            exit()

    input_tree = Tree(tree_file_in, quoted_node_names=True, format=tree_fmt)
    subset_tree = input_tree.copy()
    subset_tree.prune(genomes_to_keep, preserve_branch_length=True)
    subset_tree.write(outfile=tree_file_out)

    print('Subset tree exported to: %s' % tree_file_out)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,                       help='input tree file')
    parser.add_argument('-o',      required=True,                       help='output tree file')
    parser.add_argument('-k',      required=False, default=None,        help='leaves to keep')
    parser.add_argument('-r',      required=False, default=None,        help='leaves to remove')
    parser.add_argument('-fmt',    required=False, default=1, type=int, help='tree format, default: 1')
    args = vars(parser.parse_args())
    subset(args)
