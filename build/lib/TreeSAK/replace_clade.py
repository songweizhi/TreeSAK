import os
import argparse
from ete3 import Tree


replace_clade_usage = '''
====================== replace_clade example commands ======================

TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaves.txt
TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaf1
TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaf1,leaf2

============================================================================
'''


def replace_clade(args):

    main_tree_file  = args['m']
    main_tree_fmt   = args['mf']
    sub_tree_file   = args['s']
    sub_tree_fmt    = args['sf']
    leaf_txt        = args['l']
    tree_out        = args['o']
    tree_out_fmt    = args['of']

    if os.path.isfile(main_tree_file) is False:
        print('%s does not exist, program exited!' % main_tree_file)
        exit()

    if os.path.isfile(sub_tree_file) is False:
        print('%s does not exist, program exited!' % sub_tree_file)
        exit()

    # read in subtree
    sub_tre = Tree(sub_tree_file, quoted_node_names=True, format=sub_tree_fmt)

    # read in leaves
    leaf_list = []
    if os.path.isfile(leaf_txt) is True:
        for each_leaf in open(leaf_txt):
            leaf_list.append(each_leaf.strip())
    else:
        leaf_list = leaf_txt.split(',')

    # read in main tree
    main_tre = Tree(main_tree_file, quoted_node_names=True, format=main_tree_fmt)

    # remove clades
    if len(leaf_list) == 1:
        lca = main_tre.search_nodes(name=leaf_list[0])[0]
    else:
        lca = main_tre.get_common_ancestor(leaf_list)

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(sub_tre)

    # write out updated tree
    main_tre.write(outfile=tree_out, format=tree_out_fmt)


if __name__ == '__main__':

    replace_clade_parser = argparse.ArgumentParser()
    replace_clade_parser.add_argument('-m',   required=True,                        help='main tree file')
    replace_clade_parser.add_argument('-mf',  required=False, default=1, type=int,  help='main tree format, default: 1')
    replace_clade_parser.add_argument('-s',   required=True,                        help='subtree file')
    replace_clade_parser.add_argument('-sf',  required=False, default=1, type=int,  help='subtree format, default: 1')
    replace_clade_parser.add_argument('-l',   required=True,                        help='leaves on main tree to be replaced')
    replace_clade_parser.add_argument('-o',   required=True,                        help='output tree')
    replace_clade_parser.add_argument('-of',  required=False, default=9, type=int,  help='output tree format, default is 9')
    args = vars(replace_clade_parser.parse_args())
    replace_clade(args)
