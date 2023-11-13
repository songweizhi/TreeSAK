import os
import argparse
from ete3 import Tree


LcaToLeaves_usage = '''
==================== LcaToLeaves example commands ====================

BioSAK LcaToLeaves -s species_tree.ale.stree -n 123
BioSAK LcaToLeaves -s species_tree.ale.stree -n internal_nodes.txt

======================================================================
'''


def lca_to_two_leaves(species_tree_from_ale, internal_node_id):

    # read in ale species tree
    stree_ale = Tree(species_tree_from_ale, format=1)

    # get all leaves of the internal node
    internal_node = stree_ale.search_nodes(name=internal_node_id)[0]
    internal_node_leaf_object = internal_node.get_leaves()
    internal_node_leaf_set = set()
    for each_leaf in internal_node_leaf_object:
        internal_node_leaf_set.add(each_leaf.name)

    # get the two leaves needed
    targeted_two_leaves = []
    leaves_found = False
    for leaf_1 in internal_node_leaf_set:
        for leaf_2 in internal_node_leaf_set:
            if leaf_1 != leaf_2:
                if leaves_found is False:
                    current_lca_id = stree_ale.get_common_ancestor(leaf_1, leaf_2).name
                    if current_lca_id == internal_node_id:
                        targeted_two_leaves.append(leaf_1)
                        targeted_two_leaves.append(leaf_2)
                        leaves_found = True

    return targeted_two_leaves[0], targeted_two_leaves[1]


def LcaToLeaves(args):

    species_tree_from_ale = args['s']
    internal_node         = args['n']

    internal_node_set = set()
    if os.path.isfile(internal_node) is False:
        internal_node_set.add(internal_node)
    else:
        for each_node in open(internal_node):
            internal_node_set.add(each_node.strip())

    for each_internal_node in internal_node_set:
        leaf_1, leaf_2 = lca_to_two_leaves(species_tree_from_ale, each_internal_node)
        print('%s\t%s\t%s' % (each_internal_node, leaf_1, leaf_2))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s',   required=True,   help='the .stree file from ALE')
    parser.add_argument('-n',   required=True,   help='internal node(s)')
    args = vars(parser.parse_args())
    LcaToLeaves(args)
