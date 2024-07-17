import argparse
from ete3 import Tree


mcmc2tree_usage = '''
============ mcmc2tree example commands ============

TreeSAK mcmc2tree -i mcmctree.out -o renamed.tree

====================================================
'''


def mcmc2tree(args):

    mamctree_out = args['i']
    tree_file    = args['o']

    # get tree string from mcmctree_out
    tree_str = ''
    tree_line = 0
    current_line = 1
    for each_line in open(mamctree_out):
        if 'Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels' in each_line:
            tree_line = current_line + 1
        if tree_line == current_line:
            tree_str = each_line.strip()
        current_line += 1

    tree_str_no_space = tree_str.replace(' ', '')

    # rename tree nodes
    t = Tree(tree_str_no_space, format=1)
    for each_node in t.traverse():
        if each_node.is_leaf():
            node_name_new = '_'.join(each_node.name.split('_')[1:])
        else:
            node_name_new = 't_n%s' % each_node.name
        each_node.name = node_name_new

    tree_str_renamed = t.write(format=8)

    tree_file_handle = open(tree_file, 'w')
    tree_file_handle.write(tree_str_renamed + '\n')
    tree_file_handle.close()


if __name__ == '__main__':

    mcmc2tree_parser = argparse.ArgumentParser()
    mcmc2tree_parser.add_argument('-i',  required=True,  help='the .out file from mcmctree')
    mcmc2tree_parser.add_argument('-o',  required=True,  help='output tree file')
    args = vars(mcmc2tree_parser.parse_args())
    mcmc2tree(args)
