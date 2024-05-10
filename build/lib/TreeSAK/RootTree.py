import argparse
from ete3 import Tree


RootTree_usage = '''
====================== RootTree example commands ======================

BioSAK RootTree -i input.tree -og outgroup_genomes.txt -o rooted.tree

=======================================================================
'''

def RootTree(args):

    tree_file           = args['i']
    out_group_txt       = args['og']
    tree_file_rooted    = args['o']
    tree_fmt            = args['fmt']

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=tree_fmt)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted, format=tree_fmt)


if __name__ == '__main__':

    RootTree_parser = argparse.ArgumentParser()
    RootTree_parser.add_argument('-i',      required=True,                          help='input tree')
    RootTree_parser.add_argument('-og',     required=True,                          help='out group leaves')
    RootTree_parser.add_argument('-o',      required=True,                          help='output tree')
    RootTree_parser.add_argument('-fmt',    required=False, default=1, type=int,    help='tree format, default: 1')
    args = vars(RootTree_parser.parse_args())
    RootTree(args)
