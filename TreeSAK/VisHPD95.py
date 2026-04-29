import os
import argparse
from ete3 import Tree


VisHPD95_usage = '''
============================ VisHPD95 example command ============================

TreeSAK VisHPD95 -i samples.txt -n nodes.txt -o HPD95.pdf

# The input file has three columns (tab separated)
# Samples will be colored by the 1st column and shaped by the 2nd column
DeltaLL_50	IR  dating_results/M1_PA_75_DeltaLL_50_clock2_nsample500000_out.txt
DeltaLL_50	AR  dating_results/M1_PA_75_DeltaLL_50_clock3_nsample500000_out.txt
DeltaLL_75	IR  dating_results/M1_PA_75_DeltaLL_75_clock2_nsample500000_out.txt
DeltaLL_75	AR  dating_results/M1_PA_75_DeltaLL_75_clock3_nsample500000_out.txt
JTT_default	RelTime	dating_results/topo1_RelTime_JTT_default.txt
JTT_Gamma4	RelTime	dating_results/topo1_RelTime_JTT_Gamma4.txt

# Example data
https://github.com/songweizhi/TreeSAK/tree/master/DemoData/VisHPD95

# format of nodes.txt
GCA038869675_1,GCA013002265_1	LCA_AOA
GCF900177045_1,GCA013002265_1	LCA_f__Nitrosopumilaceae 
GCA965218725_1,GCA026706865_1	D1

==================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


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
            if len(each_line_split_no_empty) >= 7:
                node_id           = each_line_split_no_empty[0]
                value_mean        = each_line_split_no_empty[1]
                value_hpd95_small = each_line_split_no_empty[4]
                value_hpd95_big   = each_line_split_no_empty[5]
                node_to_mean_hpd95_dict[node_id] = [value_mean, value_hpd95_small, value_hpd95_big]
        current_line += 1

    return node_to_mean_hpd95_dict


def get_lca_from_reltime_op(reltime_txt, leaf_1_name, leaf_2_name):

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


def parse_reltime_op(reltime_txt, interested_nodes_txt):

    scale_factor = 1

    lca_to_leaves_dict = dict()
    interested_node_desc_dict = dict()
    if interested_nodes_txt is not None:
        if os.path.isfile(interested_nodes_txt) is False:
            print('%s not found, program exited!' % interested_nodes_txt)
            exit()
        for interested_node in open(interested_nodes_txt):
            interested_node_split = interested_node.strip().split('\t')
            paired_leaves = interested_node_split[0]
            interested_node_desc = paired_leaves
            if len(interested_node_split) > 1:
                interested_node_desc = interested_node_split[1]
            interested_node_desc_dict[paired_leaves] = interested_node_desc
            leaf_1 = paired_leaves.split(',')[0]
            leaf_2 = paired_leaves.split(',')[1]
            lca_id = get_lca_from_reltime_op(reltime_txt, leaf_1, leaf_2)
            lca_to_leaves_dict[lca_id] = paired_leaves.strip()

    print(lca_to_leaves_dict)
    print(interested_node_desc_dict)
    print()

    value_list = []
    line_num_index = 0
    for each_line in open(reltime_txt):
        each_line_split = each_line.strip().split('\t')
        each_line_split = [i.strip() for i in each_line_split]
        if line_num_index != 0:
            if len(each_line_split) > 1:
                if len(lca_to_leaves_dict) == 0:
                    value_list.append('\t'.join(each_line_split))
                else:
                    node_id = each_line_split[1]
                    if node_id in lca_to_leaves_dict:
                        corresponding_leaves = lca_to_leaves_dict[node_id]
                        interested_node_desc = interested_node_desc_dict[corresponding_leaves]
                        value_list.append('%s\t%s\t%s' % (corresponding_leaves, interested_node_desc, '\t'.join(each_line_split)))
        line_num_index += 1

    to_write_dict = dict()
    line_num_index = 0
    for each_line in value_list:
        if line_num_index > 0:
            each_line_split = each_line.strip().split('\t')
            node_label = each_line_split[1]
            div_time = each_line_split[7]
            ci_lower = each_line_split[8]
            ci_upper = each_line_split[9]
            if len(lca_to_leaves_dict) != 0:
                div_time = each_line_split[9]
                ci_lower = each_line_split[10]
                ci_upper = each_line_split[11]
            if not ((div_time == '-') and (ci_lower == '-') and (ci_upper == '-')):
                div_time = float("{0:.2f}".format(float(div_time) * scale_factor))
                ci_lower = float("{0:.2f}".format(float(ci_lower) * scale_factor))
                ci_upper = float("{0:.2f}".format(float(ci_upper) * scale_factor))
                to_write_dict[node_label] = '%s\t%s-%s\t%s' % (div_time, ci_upper, ci_lower, node_label)
                to_write_dict[node_label] = [div_time, ci_lower, ci_upper]
        line_num_index += 1

    return to_write_dict


def VisHPD95(args):

    input_txt   = args['i']
    node_txt    = args['n']
    plot_out    = args['o']
    plot_width  = args['x']
    plot_height = args['y']

    _, op_path, op_base, _ = sep_path_basename_ext(plot_out)

    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    VisHPD95_R        = '%s/VisHPD95.R'     % current_file_path
    VisHPD95_R_v2     = '%s/VisHPD95_v2.R'  % current_file_path
    dm_out            = '%s/%s_table.txt'   % (op_path, op_base)

    if os.path.isfile(input_txt) is False:
        print('%s not found, program exited!' % input_txt)
        exit()

    if os.path.isfile(node_txt) is False:
        print('%s not found, program exited!' % node_txt)
        exit()

    node_desc_list = []
    for i in open(node_txt):
        i_split = i.strip().split('\t')
        node_desc_list.append(i_split[-1])
    order_str = ','.join(node_desc_list)

    dating_result_file_list = []
    label_dict              = dict()
    color_dict              = dict()
    shape_dict              = dict()
    missing_file_set        = set()
    missing_tree_file_set   = set()
    for each_file in open(input_txt):
        each_file_split = each_file.strip().split('\t')
        if len(each_file_split) == 4:
            dating_result_file  = each_file_split[-1].strip()
            dating_result_file_list.append(dating_result_file)
            if os.path.isfile(dating_result_file) is False:
                missing_file_set.add(dating_result_file)
            label_dict[dating_result_file] = each_file_split[0]
            color_dict[dating_result_file] = each_file_split[1]
            shape_dict[dating_result_file] = each_file_split[2]
        else:
            print('There is something wrong with the format of %s' % input_txt)
            print('run "TreeSAK VisHPD95 -h" for help, Program exited!')
            exit()

    if len(missing_file_set) > 0:
        print('The following files are missing, program exited!')
        print('\n'.join(sorted(list(missing_file_set))))
        exit()

    if len(missing_tree_file_set) > 0:
        print('The following tree files are missing, program exited!')
        print('\n'.join(sorted(list(missing_tree_file_set))))
        exit()

    ########## read in dating results ##########
    value_dict = dict()
    for dating_result_file in dating_result_file_list:
        dating_result_file_name, _, _, _ = sep_path_basename_ext(dating_result_file)
        label_col_to_write  = label_dict.get(dating_result_file, dating_result_file)
        color_col_to_write  = color_dict.get(dating_result_file, dating_result_file)
        shape_col_to_write  = shape_dict.get(dating_result_file, dating_result_file)
        dating_approach = ''
        with open(dating_result_file) as f:
            first_line = f.readline()
            if first_line.startswith('MCMCTREE'):
                dating_approach = 'mcmctree'
            elif first_line.startswith('NodeLabel'):
                dating_approach = 'reltime'
            else:
                print('Unrecognisable dating approach for %s, program exited!' % dating_result_file_name)
                exit()

        if dating_approach == 'mcmctree':
            node_set, node_rename_dict, tree_str = get_internal_node_to_plot(node_txt, dating_result_file)
            node_to_mean_95_hpd_dict = read_in_posterior_mean(dating_result_file)
            for each_node in node_set:
                node_name_to_write = node_rename_dict.get(each_node, each_node)
                mean_95_hpd_list = node_to_mean_95_hpd_dict.get(each_node)
                value_str ='%s\t%s\t%s\t%s\t%s' % (label_col_to_write, color_col_to_write, shape_col_to_write, node_name_to_write, '\t'.join(mean_95_hpd_list))
                if node_name_to_write not in value_dict:
                    value_dict[node_name_to_write] = []
                value_dict[node_name_to_write].append(value_str)
        if dating_approach == 'reltime':
            reltime_value_dict = parse_reltime_op(dating_result_file, node_txt)
            for each_node in reltime_value_dict:
                node_age_value_list = reltime_value_dict[each_node]
                value_str = '%s\t%s\t%s\t%s\t%s' % (label_col_to_write, color_col_to_write, shape_col_to_write, each_node, '\t'.join([str(i) for i in node_age_value_list]))
                if each_node not in value_dict:
                    value_dict[each_node] = []
                value_dict[each_node].append(value_str)

    dm_out_handle = open(dm_out, 'w')
    dm_out_handle.write('Label\tColor\tShape\tNode\tMean\tLow\tHigh\tPost_X\n')
    x_loci = 3
    v_line_list = ['0']
    break_range_dict = dict()
    for each_node in node_desc_list:
        if each_node not in break_range_dict:
            break_range_dict[each_node] = []
        str_list = value_dict[each_node]
        for each_str in str_list:
            break_range_dict[each_node].append(x_loci)
            dm_out_handle.write('%s\t%s\n' % (each_str, x_loci))
            x_loci +=1
        v_line_list.append(str(x_loci + 3))
        x_loci += 6
    dm_out_handle.close()
    v_line_str = ','.join(v_line_list)
    break_point_list = []
    for each_range in break_range_dict:
        range_value_list = sorted(break_range_dict[each_range])
        start_point = range_value_list[0]
        end_point   = range_value_list[-1]
        break_point = start_point + ((end_point - start_point)/2)
        break_point_list.append(break_point)
    break_point_list_sorted = sorted(break_point_list)
    break_point_str = ','.join([str(i) for i in break_point_list_sorted])

    plot_cmd    = 'Rscript %s -i %s -x %s -y %s -o %s -l "%s" -b %s -v %s' % (VisHPD95_R,    dm_out, plot_width, plot_height, plot_out, order_str, break_point_str, v_line_str)
    plot_cmd_v2 = 'Rscript %s -i %s -x %s -y %s -o %s -l "%s" -b %s -v %s' % (VisHPD95_R_v2, dm_out, plot_width, plot_height, plot_out, order_str, break_point_str, v_line_str)
    os.system(plot_cmd_v2)

    print('Plot exported to: %s' % plot_out)
    print('Done!')


if __name__ == '__main__':

    VisHPD95_parser = argparse.ArgumentParser()
    VisHPD95_parser.add_argument('-i',      required=True,                          help='input file')
    VisHPD95_parser.add_argument('-n',      required=True,                          help='nodes to plot')
    VisHPD95_parser.add_argument('-x',      required=False, default=10,type=int,    help='plot width, default: 10')
    VisHPD95_parser.add_argument('-y',      required=False, default=6,type=int,     help='plot height, default: 6')
    VisHPD95_parser.add_argument('-o',      required=True,                          help='output plot')
    args = vars(VisHPD95_parser.parse_args())
    VisHPD95(args)
