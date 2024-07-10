import argparse
from ete3 import Tree


replace_clade_usage = '''
===================== replace_clade example commands =====================

TreeSAK replace_clade -m main.tree -s sub.tree -l leaves.txt -o out.tree

==========================================================================
'''


def replace_clade(args):

    main_tree_file  = args['m']
    sub_tree_file   = args['s']
    leaf_txt        = args['l']
    tree_out        = args['o']
    tree_out_fmt    = args['of']

    # read in subtree
    sub_tre = Tree(sub_tree_file, quoted_node_names=True, format=1)

    # read in leaves
    leaf_list = []
    for each_leaf in open(leaf_txt):
        leaf_list.append(each_leaf.strip())

    # read in main tree
    main_tre = Tree(main_tree_file, quoted_node_names=True, format=1)

    # remove clades
    lca = main_tre.get_common_ancestor(leaf_list)

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(sub_tre)

    # write out updated tree
    main_tre.write(outfile=tree_out, format=tree_out_fmt)




if __name__ == '__main__':

    replace_clade_parser = argparse.ArgumentParser()
    replace_clade_parser.add_argument('-m',   required=True,                        help='main tree file')
    replace_clade_parser.add_argument('-s',   required=True,                        help='subtree file')
    replace_clade_parser.add_argument('-l',   required=True,                        help='leaves on main tree to be replaced')
    replace_clade_parser.add_argument('-o',   required=True,                        help='output tree')
    replace_clade_parser.add_argument('-of',  required=False, default=9, type=int,  help='output tree format, default is 9')
    args = vars(replace_clade_parser.parse_args())
    replace_clade(args)
