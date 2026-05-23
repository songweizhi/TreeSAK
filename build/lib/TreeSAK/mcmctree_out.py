import os
import glob
import argparse
from ete3 import Tree


mcmctree_out_usage = '''
===================== mcmctree_out example commands =====================

TreeSAK mcmctree_out -i topo1_out.txt -n important_nodes.txt -o op.txt 
TreeSAK mcmctree_out -i dating_result -n important_nodes.txt -o op.txt

# The input file need to be in the format of *_out.txt

# Format of node file (tab separated)

t_n141	Oxygen Age Constraint, Thermoproteales
t_n201	Emergence of the D1 clade
Povist_bin29,SMArc018	LCA_D7
SMArc016,GCA964507825_1	LCA_D8

=========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def get_node_age_dict(out_txt):
    node_age_dict = dict()
    keep_the_next_value_line = False
    for each_line in open(out_txt):
        if 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width' in each_line:
            keep_the_next_value_line = True
        else:
            if keep_the_next_value_line is True:
                each_line_split = each_line.strip().replace('(', '').replace(')', '').split()
                if len(each_line_split) > 0:
                    node_id = each_line_split[0]
                    node_mean = each_line_split[1].replace(',', '')
                    hpd_ci_95_l = each_line_split[4].replace(',', '')
                    hpd_ci_95_r = each_line_split[5].replace(',', '')
                    value_str = '%s\t%s-%s' % (node_mean, hpd_ci_95_l, hpd_ci_95_r)
                    node_age_dict[node_id] = value_str
    return node_age_dict


def mcmc2tree(mamctree_out):

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

    # rename internal nodes
    t = Tree(tree_str_with_internal_node_id, format=1)
    for each_node in t.traverse():
        if each_node.is_leaf():
            node_name_new = '_'.join(each_node.name.split('_')[1:])
        else:
            node_name_new = 't_n%s' % each_node.name
        each_node.name = node_name_new

    # rename tree nodes
    t = Tree(tree_str_with_internal_node_id, format=1)
    for each_node in t.traverse():
        if each_node.is_leaf():
            node_name_new = '_'.join(each_node.name.split('_')[1:])
        else:
            node_name_new = 't_n%s' % each_node.name
        each_node.name = node_name_new

    tree_str = str(t.write(format=8))

    return tree_str


def get_lca_of_leaves(tree_file, leaf_1_name, leaf_2_name):

    tree_subject = Tree(tree_file, quoted_node_names=True, format=1)
    leaf1_node = tree_subject.search_nodes(name=leaf_1_name)
    leaf2_node = tree_subject.search_nodes(name=leaf_2_name)

    leaf1 = leaf1_node[0]
    leaf2 = leaf2_node[0]

    lca_node = leaf1.get_common_ancestor(leaf2)
    lca_node_name = lca_node.name

    return lca_node_name


def mcmctree_out(args):

    mcmc_out_dir         = args['i']
    interested_nodes_txt = args['n']
    op_txt               = args['o']

    # read in interested_nodes_txt
    interested_node_dict = dict()
    for each_line in open(interested_nodes_txt):
        each_line_split = each_line.strip().split('\t')
        if len(each_line_split) > 1:
            node_id   = each_line_split[0]
            node_desc = each_line_split[1]
            interested_node_dict[node_id] = node_desc

    tree_dict = dict()
    node_age_dod = dict()
    if os.path.isfile(mcmc_out_dir) is True:
        _, _, f_base, _ = sep_path_basename_ext(mcmc_out_dir)
        node_age_dict  = get_node_age_dict(mcmc_out_dir)
        node_age_dod[f_base] = node_age_dict
        tree_str = mcmc2tree(mcmc_out_dir)
        tree_dict[f_base] = tree_str
    elif os.path.isdir(mcmc_out_dir) is True:
        mcmcout_file_re = '%s/*_out.txt' % mcmc_out_dir
        mcmcout_file_list = glob.glob(mcmcout_file_re)

        if len(mcmcout_file_list) == 0:
            print('No file found in %s, program exited!' % mcmc_out_dir)
            exit()

        for each_mcmcout_file in mcmcout_file_list:
            _, _, f_base, _ = sep_path_basename_ext(each_mcmcout_file)
            current_node_age_dict = get_node_age_dict(each_mcmcout_file)
            node_age_dod[f_base] = current_node_age_dict
            tree_str = mcmc2tree(each_mcmcout_file)
            tree_dict[f_base] = tree_str
    else:
        print('No input file found, program exited!')
        exit()

    node_list_sorted = sorted(list(interested_node_dict.keys()))
    node_list_sorted_with_desc = [('%s__%s' % (i, interested_node_dict[i])) for i in node_list_sorted]

    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('Node\tPosterior_mean\tHPD_CI_95\tDescription\n')
    if len(node_age_dod) ==1:
        for each_setting in node_age_dod:
            current_setting_tree_str = tree_dict[each_setting]
            current_value_dict       = node_age_dod[each_setting]
            for each_node_id in node_list_sorted:
                if ',' in each_node_id:
                    each_node_id_split = each_node_id.split(',')
                    lca = get_lca_of_leaves(current_setting_tree_str, each_node_id_split[0], each_node_id_split[1])
                    current_node_value = current_value_dict[lca]
                else:
                    current_node_value = current_value_dict[each_node_id]
                op_txt_handle.write('%s\t%s\t%s\n' % (each_node_id, current_node_value, interested_node_dict.get(each_node_id, '')))

    elif len(node_age_dod) >=2:
        op_txt_handle.write('\t' + '\t'.join(node_list_sorted_with_desc) + '\n')
        for each_setting in node_age_dod:
            current_setting_tree_str = tree_dict[each_setting]
            current_value_dict = node_age_dod[each_setting]
            current_value_list = [each_setting]
            current_mean_list  = [each_setting]
            for each_node_id in node_list_sorted:
                if ',' in each_node_id:
                    each_node_id_split = each_node_id.split(',')
                    lca = get_lca_of_leaves(current_setting_tree_str, each_node_id_split[0], each_node_id_split[1])
                    current_node_value = current_value_dict[lca]
                else:
                    current_node_value = current_value_dict[each_node_id]
                current_node_mean = current_node_value.split('(')[0]
                current_value_list.append(current_node_value)
                current_mean_list.append(current_node_mean)
            op_txt_handle.write('\t'.join(current_value_list) + '\n')
    op_txt_handle.close()


if __name__ == '__main__':

    mcmctree_out_parser = argparse.ArgumentParser()
    mcmctree_out_parser.add_argument('-i', required=True,  help='MCMCTree out file/dir')
    mcmctree_out_parser.add_argument('-n', required=True,  help='interested nodes')
    mcmctree_out_parser.add_argument('-o', required=True,  help='output table')
    args = vars(mcmctree_out_parser.parse_args())
    mcmctree_out(args)
