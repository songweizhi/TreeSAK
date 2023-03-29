from ete3 import Tree


def replace_clades(main_tree, sub_tree, tree_out):

    # read in sub tree
    tre_sub = Tree(sub_tree, format=1)

    # get all leaves in sub tree
    subtree_leaf_name_list = tre_sub.get_leaf_names()

    # read in main tree
    tre_main = Tree(main_tree)

    # remove clades
    lca = tre_main.get_common_ancestor(subtree_leaf_name_list)

    if len(lca.get_leaf_names()) != len(subtree_leaf_name_list):
        print('LCA of subtree leaves in main tree contain extra leaves, program exited!')
        exit()

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(tre_sub)

    # write out updated tree
    tre_main.write(outfile=tree_out, format=8)

