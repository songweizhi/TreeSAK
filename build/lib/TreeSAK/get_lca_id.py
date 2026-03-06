import os
import argparse
from ete3 import Tree


get_lca_id_usage = '''
=================== get_lca_id example commands ===================

TreeSAK get_lca_id -i mcmc.tree -l GCF000022205_1,GCF000172995_2
TreeSAK get_lca_id -i mcmc.tree -l leaves.txt

===================================================================
'''


def get_lca_of_leaves(tree_subject, leaf_1_name, leaf_2_name):

    leaf1_node = tree_subject.search_nodes(name=leaf_1_name)
    leaf2_node = tree_subject.search_nodes(name=leaf_2_name)

    leaf1 = leaf1_node[0]
    leaf2 = leaf2_node[0]

    lca_node = leaf1.get_common_ancestor(leaf2)
    lca_node_name = lca_node.name

    return lca_node_name


def get_lca_id(args):

    tree_file  = args['i']
    leaf_txt   = args['l']

    if os.path.isfile(tree_file) is False:
        print('Error: input tree file does not exist')
        exit()

    tree_subject = Tree(tree_file, quoted_node_names=True, format=1)

    if os.path.isfile(leaf_txt) is False:
        leaf_split = leaf_txt.split(',')
        lca_node_name = get_lca_of_leaves(tree_subject, leaf_split[0], leaf_split[1])
        print('%s\t%s' % (leaf_txt.strip(), lca_node_name))
    else:
        for each_line in open(leaf_txt):
            leaf_split = each_line.strip().split(',')
            lca_node_name = get_lca_of_leaves(tree_subject, leaf_split[0], leaf_split[1])
            print('%s\t%s' % (each_line.strip(), lca_node_name))


if __name__ == '__main__':

    get_lca_id_parser = argparse.ArgumentParser()
    get_lca_id_parser.add_argument('-i', required=True,  help='input tree file')
    get_lca_id_parser.add_argument('-l', required=True,  help='leaves file')
    args = vars(get_lca_id_parser.parse_args())
    get_lca_id(args)
