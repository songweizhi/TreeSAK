import argparse
from ete3 import Tree


root_tree_usage = '''
====================== root_tree example commands ======================

BioSAK root_tree -i input.tree -og outgroup_genomes.txt -o rooted.tree

========================================================================
'''


# root_with_out_group
def root_tree(args):

    tree_file           = args['i']
    out_group_txt       = args['og']
    tree_file_rooted    = args['o']

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=1)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted)


if __name__ == '__main__':

    root_tree_parser = argparse.ArgumentParser()
    root_tree_parser.add_argument('-i',     required=True,                       help='input tree')
    root_tree_parser.add_argument('-og',    required=True,                       help='rename file')
    root_tree_parser.add_argument('-o',     required=True,                       help='output tree')
    args = vars(root_tree_parser.parse_args())
    root_tree(args)
