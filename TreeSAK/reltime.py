import os
import argparse


reltime_usage = '''
==================== reltime example commands ====================

TreeSAK reltime -i RelTime.txt -o op.txt -scale 1000
TreeSAK reltime -i RelTime.txt -o op.txt -n interested_nodes.txt

# Format of interested_nodes.txt
COS36386_bin19,GCA026706865_1	LCA_shallow-water_D1
GCA965218725_1,COS36386_bin19	LCA_D1
GCF025699255_1,GCA965218725_1	LCA_D1_diverging

==================================================================
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


def reltime(args):

    reltime_txt          = args['i']
    interested_nodes_txt = args['n']
    op_txt               = args['o']
    scale_factor         = args['scale']

    _, f_path, f_base, _ = sep_path_basename_ext(op_txt)
    op_txt_all_info      = '%s/%s_all_info.txt' % (f_path,f_base)

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
            lca_id = get_lca(reltime_txt, leaf_1, leaf_2)
            lca_to_leaves_dict[lca_id] = paired_leaves.strip()

    op_txt_all_info_handle = open(op_txt_all_info, 'w')
    line_num_index = 0
    for each_line in open(reltime_txt):
        each_line_split = each_line.strip().split('\t')
        each_line_split = [i.strip() for i in each_line_split]
        if line_num_index == 0:
            op_txt_all_info_handle.write('Leaves\tDescription\t%s\n' % ('\t'.join(each_line_split)))
        else:
            if len(each_line_split) > 1:
                if len(lca_to_leaves_dict) == 0:
                    op_txt_all_info_handle.write('%s\n' % ('\t'.join(each_line_split)))
                else:
                    node_id = each_line_split[1]
                    if node_id in lca_to_leaves_dict:
                        corresponding_leaves = lca_to_leaves_dict[node_id]
                        interested_node_desc = interested_node_desc_dict[corresponding_leaves]
                        op_txt_all_info_handle.write('%s\t%s\t%s\n' % (corresponding_leaves, interested_node_desc, '\t'.join(each_line_split)))
        line_num_index += 1
    op_txt_all_info_handle.close()

    to_write_dict = dict()
    line_num_index = 0
    for each_line in open(op_txt_all_info):
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
        line_num_index += 1

    # write out
    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('DivTime\tCI95(Upper-Lower)\tNode\n')
    for each_node in sorted(list(to_write_dict.keys())):
        op_txt_handle.write(to_write_dict[each_node] + '\n')
    op_txt_handle.close()

    os.system('rm %s' % op_txt_all_info)


if __name__ == '__main__':

    reltime_parser = argparse.ArgumentParser()
    reltime_parser.add_argument('-i',       required=True,                          help='reltime output file')
    reltime_parser.add_argument('-n',       required=False, default=None,           help='interested node txt')
    reltime_parser.add_argument('-scale',   required=False, default=1, type=int,    help='scale the value, default is 1 (no scaling)')
    reltime_parser.add_argument('-o',       required=True,                          help='output txt file')
    args = vars(reltime_parser.parse_args())
    reltime(args)
