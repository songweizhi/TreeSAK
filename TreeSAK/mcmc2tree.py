import os
import argparse
from ete3 import Tree


mcmc2tree_usage = '''
=============== mcmc2tree example commands ===============

TreeSAK mcmc2tree -i mcmctree.out

# This module will extract two trees from the input file
1. one tree with mcmctree added internal node ids
2. another tree with node ages

=========================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def mcmc2tree(args):

    mamctree_out            = args['i']

    # define output file name
    _, _, file_in_base, _ = sep_path_basename_ext(mamctree_out)

    op_tree_with_internal_node  = '%s_with_mcmctree_internal_node.treefile' % file_in_base
    op_tree_with_node_age       = '%s_with_node_age.treefile'               % file_in_base

    if os.path.isfile(mamctree_out) is False:
        print('%s not found, program exited!' % mamctree_out)

    # get tree string from mcmctree_out
    tree_str_list = []
    add_str_to_list = False
    for each_line in open(mamctree_out):
        each_line = each_line.strip()
        if 'Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels' in each_line:
            add_str_to_list = True
        if 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width' in each_line:
            add_str_to_list = False

        if add_str_to_list is True:
            if (not each_line.startswith('Species tree for FigTree.')) and (len(each_line) > 0):
                tree_str_list.append(each_line.replace(' ', ''))

    tree_str_with_internal_node_id  = tree_str_list[0]
    tree_str_with_node_age          = tree_str_list[-1]

    # rename internal nodes
    t = Tree(tree_str_with_internal_node_id, format=1)
    for each_node in t.traverse():
        if each_node.is_leaf():
            node_name_new = '_'.join(each_node.name.split('_')[1:])
        else:
            node_name_new = 't_n%s' % each_node.name
        each_node.name = node_name_new

    # write out tree_str_with_node_age
    op_tree_with_node_age_handle = open(op_tree_with_node_age, 'w')
    op_tree_with_node_age_handle.write(tree_str_with_node_age + '\n')
    op_tree_with_node_age_handle.close()

    # rename tree nodes
    t = Tree(tree_str_with_internal_node_id, format=1)
    for each_node in t.traverse():
        if each_node.is_leaf():
            node_name_new = '_'.join(each_node.name.split('_')[1:])
        else:
            node_name_new = 't_n%s' % each_node.name
        each_node.name = node_name_new

    op_tree_with_internal_node_handle = open(op_tree_with_internal_node, 'w')
    op_tree_with_internal_node_handle.write(t.write(format=8) + '\n')
    op_tree_with_internal_node_handle.close()


if __name__ == '__main__':

    mcmc2tree_parser = argparse.ArgumentParser()
    mcmc2tree_parser.add_argument('-i',  required=True,  help='the .out file from mcmctree')
    args = vars(mcmc2tree_parser.parse_args())
    mcmc2tree(args)
