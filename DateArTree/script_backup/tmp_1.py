import os
import glob
import argparse
import pandas as pd
from ete3 import Tree
import plotly.express as px


def mamctree_out_to_tree_str(mamctree_out):

    # get tree string from mamctree_out
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

    return tree_str_renamed


def get_internal_node_to_plot(node_txt, mo_file):

    tree_str = ''
    if os.path.isfile(mo_file):
        tree_str = mamctree_out_to_tree_str(mo_file)

    # get nodes to plot
    node_set = set()
    node_rename_dict = dict()
    if os.path.isfile(node_txt) is True:
        for each in open(node_txt):
            each_split = each.strip().split('\t')
            node_str = each_split[0]

            # get internal_node_to_plot
            internal_node_to_plot = ''
            if ',' not in node_str:
                internal_node_to_plot = each_split[0]
            else:
                leaf_list = node_str.split(',')
                if tree_str == '':
                    print('MCMCTree out file not found, program exited!')
                    exit()
                current_lca = Tree(tree_str, format=1).get_common_ancestor(leaf_list)
                internal_node_to_plot = current_lca.name

            # add internal_node_to_plot to  node_set
            if internal_node_to_plot != '':
                node_set.add(internal_node_to_plot)

            # read in name to show in plot
            if len(each_split) == 2:
                if each_split[1] != '':
                    node_rename_dict[internal_node_to_plot] = each_split[1]
    else:
        node_set = node_txt.split(',')

    return node_set, node_rename_dict

# mamctree_out = '/Users/songweizhi/Desktop/777/PA_75_DeltaLL_50_clock3_nsample500000_out.txt'
# mamctree_formatted_tree_str = mamctree_out_to_tree_str(mamctree_out)
# print('mamctree_formatted_tree_str\t%s' % mamctree_formatted_tree_str)

df_txt = '/Users/songweizhi/Desktop/777/multi_run_multi_nodes.pdf.txt'
df = pd.read_table(df_txt, sep=',')

setting_list = df['Setting'].unique()
setting_list_sorted = sorted(setting_list)

print(setting_list_sorted)