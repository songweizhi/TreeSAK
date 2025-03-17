import os
import argparse
from ete3 import Tree


mcmctree_vs_reltime_usage = '''
====================== mcmctree_vs_reltime example command ======================

TreeSAK mcmctree_vs_reltime -m mcmc.out -r reltime.txt -n nodes.txt -o ages.pdf

# Example data
https://github.com/songweizhi/TreeSAK/tree/master/DemoData/mcmctree_vs_reltime

=================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_lca(reltime_txt, leaf_1_name, leaf_2_name):

    leaf_set = set()
    child_to_parent_dict = dict()
    id_to_name_dict = dict()
    name_to_id_dict = dict()
    for each_line in open(reltime_txt):
        if not each_line.startswith('NodeLabel'):
            each_line_split = each_line.strip().split('\t')
            each_line_split = [i.strip() for i in each_line_split]
            if len(each_line_split) > 1:
                node_name = each_line_split[0].replace(' ', '_')
                node_id = each_line_split[1]
                des1 = each_line_split[2]
                des2 = each_line_split[3]
                id_to_name_dict[node_id] = node_name
                name_to_id_dict[node_name] = node_id
                child_to_parent_dict[des1] = node_id
                child_to_parent_dict[des2] = node_id
                if (des1 == '-') and (des2 == '-'):
                    leaf_set.add(node_id)

    leaf_to_lineage_dict = dict()
    for leaf in sorted([i for i in leaf_set]):
        original_leaf = leaf
        lineage_list = [leaf]
        while leaf in child_to_parent_dict:
            leaf_p = child_to_parent_dict[leaf]
            lineage_list.append(leaf_p)
            leaf = leaf_p
        leaf_to_lineage_dict[original_leaf] = lineage_list

    leaf_1_id     = name_to_id_dict[leaf_1_name]
    leaf_2_id     = name_to_id_dict[leaf_2_name]
    leaf_1_linage = leaf_to_lineage_dict[leaf_1_id]
    leaf_2_linage = leaf_to_lineage_dict[leaf_2_id]

    lca = ''
    for each_p in leaf_1_linage[::-1]:
        if each_p in leaf_2_linage:
            lca = each_p

    return lca


def parse_reltime(reltime_txt, interested_nodes_txt, op_txt):

    lca_to_leaves_dict = dict()
    interested_node_desc_dict = dict()
    for interested_node in open(interested_nodes_txt):
        interested_node_split = interested_node.strip().split('\t')
        paired_leaves = interested_node_split[0]
        interested_node_desc = paired_leaves
        if len(interested_node_split) > 1:
            interested_node_desc = interested_node_split[1]
        interested_node_desc_dict[paired_leaves] = interested_node_desc
        leaf_1 = paired_leaves.split(',')[0]
        leaf_2 = paired_leaves.split(',')[1]
        lca_id = get_lca(reltime_txt, leaf_1, leaf_2)
        lca_to_leaves_dict[lca_id] = paired_leaves.strip()

    op_txt_handle = open(op_txt, 'w')
    line_num_index = 0
    for each_line in open(reltime_txt):
        each_line_split = each_line.strip().split('\t')
        each_line_split = [i.strip() for i in each_line_split]
        if line_num_index == 0:
            op_txt_handle.write('ColorBy\tNode\tMean\tLow\tHigh\n')
        else:
            if len(each_line_split) > 1:
                node_id = each_line_split[1]
                if node_id in lca_to_leaves_dict:
                    node_id = each_line_split[1]
                    div_time = each_line_split[7]
                    ci_lower = each_line_split[8]
                    ci_upper = each_line_split[9]
                    corresponding_leaves = lca_to_leaves_dict[node_id]
                    interested_node_desc = interested_node_desc_dict[corresponding_leaves]
                    op_txt_handle.write('RelTime\t%s\t%s\t%s\t%s\n' % (interested_node_desc, div_time, ci_lower, ci_upper))
        line_num_index += 1
    op_txt_handle.close()


def mcmctree_out_to_tree_str(mamctree_out):

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
        tree_str = mcmctree_out_to_tree_str(mo_file)

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

    return node_set, node_rename_dict, tree_str


def read_in_posterior_mean(mcmctree_out):

    # read in Posterior mean
    node_to_mean_hpd95_dict = dict()
    current_line = 1
    posterior_mean_header_line = 0
    for each_line in open(mcmctree_out):
        if 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width' in each_line:
            posterior_mean_header_line = current_line

        if (posterior_mean_header_line != 0) and (current_line > posterior_mean_header_line):
            each_line_split = each_line.strip().split(' ')

            each_line_split_no_empty = []
            for each_element in each_line_split:
                if each_element not in ['', '(']:
                    each_element_value = each_element.replace('(', '').replace(')', '').replace(',', '')
                    each_line_split_no_empty.append(each_element_value)
            if len(each_line_split_no_empty) == 9:
                node_id           = each_line_split_no_empty[0]
                value_mean        = each_line_split_no_empty[1]
                value_hpd95_small = each_line_split_no_empty[4]
                value_hpd95_big   = each_line_split_no_empty[5]
                node_to_mean_hpd95_dict[node_id] = [value_mean, value_hpd95_small, value_hpd95_big]
        current_line += 1

    return node_to_mean_hpd95_dict


def parse_mcmc_out(mcmc_out_file, node_txt, dm_out):

    dm_out_handle = open(dm_out, 'a')
    #dm_out_handle.write('Test\tShape\tVar\tMean\tLow\tHigh\n')
    node_set, node_rename_dict, tree_str = get_internal_node_to_plot(node_txt, mcmc_out_file)
    node_to_mean_95_hpd_dict = read_in_posterior_mean(mcmc_out_file)
    for each_node in node_set:
        node_name_to_write = node_rename_dict.get(each_node, each_node)
        mean_95_hpd_list = node_to_mean_95_hpd_dict.get(each_node)
        dm_out_handle.write('MCMCTree\t%s\t%s\n' % (node_name_to_write, '\t'.join(mean_95_hpd_list)))
    dm_out_handle.close()


def mcmctree_vs_reltime(args):

    mcmc_out_file           = args['m']
    reltime_txt             = args['r']
    interested_nodes_txt    = args['n']
    pdf_out                 = args['o']

    dm_out_combined = '.'.join(pdf_out.split('.')[:-1]) + '.txt'

    # define mcmctree_vs_reltime_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    mcmctree_vs_reltime_R  = '%s/mcmctree_vs_reltime.R' % current_file_path

    parse_reltime(reltime_txt, interested_nodes_txt, dm_out_combined)
    parse_mcmc_out(mcmc_out_file, interested_nodes_txt, dm_out_combined)

    plot_cmd   = 'Rscript %s -i %s -x %s -y %s -o %s' % (mcmctree_vs_reltime_R, dm_out_combined, 8, 5, pdf_out)
    os.system(plot_cmd)
    print('Plot exported to: %s' % pdf_out)


if __name__ == '__main__':

    mcmctree_vs_reltime_parser = argparse.ArgumentParser()
    mcmctree_vs_reltime_parser.add_argument('-m',   required=True,  help='.out file from MCMCTree')
    mcmctree_vs_reltime_parser.add_argument('-r',   required=True,  help='output from elTime')
    mcmctree_vs_reltime_parser.add_argument('-n',   required=True,  help='interested nodes txt file')
    mcmctree_vs_reltime_parser.add_argument('-o',   required=True,  help='output pdf')
    args = vars(mcmctree_vs_reltime_parser.parse_args())
    mcmctree_vs_reltime(args)
