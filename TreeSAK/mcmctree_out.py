import os
import glob
import argparse


mcmctree_out_usage = '''
===================== mcmctree_out example commands =====================

TreeSAK mcmctree_out -i topo1_out.txt -n important_nodes.txt -o op.txt 
TreeSAK mcmctree_out -i dating_result -n important_nodes.txt -o op.txt

# The input file need to be in the format of *_out.txt

# Format of node file
t_n141	Oxygen Age Constraint, Thermoproteales
t_n201	Emergence of the D1 clade

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
                    value_str = '%s(%s %s)' % (node_mean, hpd_ci_95_l, hpd_ci_95_r)
                    node_age_dict[node_id] = value_str
    return node_age_dict


def mcmctree_out(args):

    mcmc_out_dir            = args['i']
    interested_nodes_txt    = args['n']
    op_txt                  = args['o']

    _, op_f_path, op_f_base, op_f_ext = sep_path_basename_ext(op_txt)
    op_txt_with_hpd_ci_95 = '%s/%s_HPD_CI_95.%s' % (op_f_path, op_f_base, op_f_ext)

    # read in interested_nodes_txt
    interested_node_dict = dict()
    for each_line in open(interested_nodes_txt):
        each_line_split = each_line.strip().split('\t')
        node_id   = each_line_split[0]
        node_desc = each_line_split[1]
        interested_node_dict[node_id] = node_desc

    node_age_dod = dict()
    if os.path.isfile(mcmc_out_dir) is True:
        _, _, f_base, _ = sep_path_basename_ext(mcmc_out_dir)
        node_age_dict  = get_node_age_dict(mcmc_out_dir)
        node_age_dod[f_base] = node_age_dict
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
    else:
        print('No input file found, program exited!')
        exit()

    node_list_sorted = sorted(list(interested_node_dict.keys()))
    node_list_sorted_with_desc = [('%s__%s' % (i, interested_node_dict[i])) for i in node_list_sorted]

    op_txt_handle = open(op_txt_with_hpd_ci_95, 'w')
    op_txt_mean_handle = open(op_txt, 'w')
    if len(node_age_dod) ==1:
        for each_setting in node_age_dod:
            current_value_dict = node_age_dod[each_setting]
            for each_node_id in node_list_sorted:
                current_node_value = current_value_dict[each_node_id]
                current_node_mean = current_node_value.split('(')[0]
                op_txt_handle.write('%s\t%s\t%s\n'      % (each_node_id, current_node_value, interested_node_dict.get(each_node_id, '')))
                op_txt_mean_handle.write('%s\t%s\t%s\n' % (each_node_id, current_node_mean, interested_node_dict.get(each_node_id, '')))

    elif len(node_age_dod) >=2:
        op_txt_handle.write('\t' + '\t'.join(node_list_sorted_with_desc) + '\n')
        op_txt_mean_handle.write('\t' + '\t'.join(node_list_sorted_with_desc) + '\n')
        for each_setting in node_age_dod:
            current_value_dict = node_age_dod[each_setting]
            current_value_list = [each_setting]
            current_mean_list = [each_setting]
            for each_node_id in node_list_sorted:
                current_node_value = current_value_dict[each_node_id]
                current_node_mean = current_node_value.split('(')[0]
                current_value_list.append(current_node_value)
                current_mean_list.append(current_node_mean)
            op_txt_handle.write('\t'.join(current_value_list) + '\n')
            op_txt_mean_handle.write('\t'.join(current_mean_list) + '\n')
    op_txt_handle.close()
    op_txt_mean_handle.close()


if __name__ == '__main__':

    mcmctree_out_parser = argparse.ArgumentParser()
    mcmctree_out_parser.add_argument('-i', required=True,  help='MCMCTree out file/dir')
    mcmctree_out_parser.add_argument('-n', required=True,  help='interested nodes')
    mcmctree_out_parser.add_argument('-o', required=True,  help='output table')
    args = vars(mcmctree_out_parser.parse_args())
    mcmctree_out(args)
