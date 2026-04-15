import os
import argparse
from ete3 import Tree


replace_clade_usage = '''
====================== replace_clade example commands ======================

TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaves.txt
TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaf1
TreeSAK replace_clade -m main.tree -s sub.tree -o out.tree -l leaf1,leaf2

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

python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/replace_clade.py


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
    batch_file      = args['b']

    if os.path.isfile(main_tree_file) is False:
        print('%s does not exist, program exited!' % main_tree_file)
        exit()

    # read in main tree
    main_tre = Tree(main_tree_file, quoted_node_names=True, format=main_tree_fmt)

    if batch_file is None:

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

        # remove clades
        if len(leaf_list) == 1:
            lca = main_tre.search_nodes(name=leaf_list[0])[0]
        else:
            lca = main_tre.get_common_ancestor(leaf_list)

        lca_p = lca.up
        lca_p.remove_child(lca)
        lca_p.add_child(sub_tre)

    else:
        if os.path.isfile(batch_file) is False:
            print('%s does not exist, program exited!' % batch_file)
            exit()

        subtrees_not_found = []
        replace_dict = dict()
        for each_line in open(batch_file):
            each_line_split = each_line.strip().split('\t')
            leaves_to_replace = each_line_split[0]
            currents_subset_tree = each_line_split[1]
            replace_dict[leaves_to_replace] = currents_subset_tree
            if os.path.isfile(currents_subset_tree) is False:
                subtrees_not_found.append(currents_subset_tree)

        if len(subtrees_not_found) > 0:
            print('The following subtrees do not exist in the batch file, program exited!')
            print('\n'.join(subtrees_not_found))
            exit()

        for each_replace in replace_dict:

            current_subtree_file    = replace_dict[each_replace]
            current_leaf_list       = each_replace.split(',')
            current_sub_tre         = Tree(current_subtree_file, quoted_node_names=True, format=sub_tree_fmt)

            # remove clades
            if len(current_leaf_list) == 1:
                lca = main_tre.search_nodes(name=current_leaf_list[0])[0]
            else:
                lca = main_tre.get_common_ancestor(current_leaf_list)

            lca_p = lca.up
            lca_p.remove_child(lca)
            lca_p.add_child(current_sub_tre)

    # write out updated tree
    main_tre.write(outfile=tree_out, format=tree_out_fmt)


if __name__ == '__main__':

    replace_clade_parser = argparse.ArgumentParser()
    replace_clade_parser.add_argument('-m',   required=True,                        help='main tree file')
    replace_clade_parser.add_argument('-mf',  required=False, default=1, type=int,  help='main tree format, default: 1')
    replace_clade_parser.add_argument('-s',   required=True,                        help='subtree file')
    replace_clade_parser.add_argument('-b',   required=False, default=None,         help='replace multiple clades in batch manner')
    replace_clade_parser.add_argument('-sf',  required=False, default=1, type=int,  help='subtree format, default: 1')
    replace_clade_parser.add_argument('-l',   required=True,                        help='leaves on main tree to be replaced')
    replace_clade_parser.add_argument('-o',   required=True,                        help='output tree')
    replace_clade_parser.add_argument('-of',  required=False, default=9, type=int,  help='output tree format, default is 9')
    args = vars(replace_clade_parser.parse_args())
    replace_clade(args)
