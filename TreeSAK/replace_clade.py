import os
import argparse
from ete3 import Tree


replace_clade_usage = '''
======================= replace_clade example commands =======================

TreeSAK replace_clade -m main.tree -s sub.tree -l leaves.txt -o output.tree

==============================================================================
'''


def replace_clade(args):

    main_tree_file  = args['m']
    sub_tree_file   = args['s']
    tree_out        = args['o']
    leaf_txt        = args['l']

    # read in sub tree
    sub_tre = Tree(sub_tree_file, quoted_node_names=True, format=1)

    # read in leaves
    leaf_list = []
    for each_leaf in open(leaf_txt):
        leaf_list.append(each_leaf.strip())

    # read in main tree
    main_tre = Tree(main_tree_file, quoted_node_names=True, format=1)

    # remove clades
    lca = main_tre.get_common_ancestor(leaf_list)

    # if len(lca.get_leaf_names()) != len(subtree_leaf_name_list):
    #     print('LCA of subtree leaves in main tree contain extra leaves, program exited!')
    #     exit()

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(sub_tre)

    # write out updated tree
    main_tre.write(outfile=tree_out, format=8)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-m',   required=True,  help='main file')
    parser.add_argument('-s',   required=True,  help='sub file')
    parser.add_argument('-l',   required=True,  help='main tree leaves')
    parser.add_argument('-o',   required=True,  help='output tree')
    args = vars(parser.parse_args())
    replace_clade(args)
