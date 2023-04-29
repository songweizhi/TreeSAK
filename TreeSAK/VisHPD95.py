import os
import glob
import argparse
from ete3 import Tree
from TreeSAK.TreeSAK_config import config_dict


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


VisHPD95_usage = '''
============================ VisHPD95 example command ============================

TreeSAK VisHPD95 -i mcmc_out -o HPD95.pdf -n nodes.txt -label label.txt
TreeSAK VisHPD95 -i mcmc_out -o HPD95.pdf -n nodes.txt -label label.txt -x 9 -y 6

# Example data
https://github.com/songweizhi/TreeSAK/tree/master/example_data/VisHPD95

==================================================================================
'''


def VisHPD95(args):

    mcmc_in     = args['i']
    node_txt    = args['n']
    label_txt   = args['label']
    plot_out    = args['o']
    plot_width  = args['x']
    plot_height = args['y']
    VisHPD95_R  = args['VisHPD95_R']

    dm_out      = '%s.txt' % plot_out

    # check MCMCTree output file/dir
    if os.path.isfile(mcmc_in) is True:
        mcmc_out_file_list = [mcmc_in]
    else:
        mcmc_out_file_re = '%s/*_out.txt' % (mcmc_in)
        mcmc_out_file_list = glob.glob(mcmc_out_file_re)

    if len(mcmc_out_file_list) == 0:
        print('MCMCTree out file not found, program exited!')
        exit()

    # read in y-axis label file
    label_dict = dict()
    color_dict = dict()
    shape_dict = dict()
    if label_txt is not None:
        for each_sample in open(label_txt):
            each_sample_split = each_sample.strip().split('\t')
            if len(each_sample_split) == 3:
                label_dict[each_sample_split[0]] = each_sample_split[1]
                color_dict[each_sample_split[0]] = each_sample_split[1]
                shape_dict[each_sample_split[0]] = each_sample_split[2]
            else:
                print('Format error: %s' % label_txt)
                exit()

    dm_out_handle = open(dm_out, 'w')
    dm_out_handle.write('Test\tShape\tVar\tMean\tLow\tHigh\n')
    for mcmc_out_file in mcmc_out_file_list:
        mcmc_out_file_no_path = mcmc_out_file
        if '/' in mcmc_out_file_no_path:
            mcmc_out_file_no_path = mcmc_out_file_no_path.split('/')[-1]

        color_col_to_write = color_dict.get(mcmc_out_file_no_path, mcmc_out_file_no_path)
        shape_col_to_write = shape_dict.get(mcmc_out_file_no_path, mcmc_out_file_no_path)
        node_set, node_rename_dict, tree_str = get_internal_node_to_plot(node_txt, mcmc_out_file)
        node_to_mean_95_hpd_dict = read_in_posterior_mean(mcmc_out_file)

        for each_node in node_set:
            node_name_to_write = node_rename_dict.get(each_node, each_node)
            mean_95_hpd_list = node_to_mean_95_hpd_dict.get(each_node)
            dm_out_handle.write('%s\t%s\t%s\t%s\n' % (color_col_to_write, shape_col_to_write, node_name_to_write, '\t'.join(mean_95_hpd_list)))
    dm_out_handle.close()

    plot_cmd   = 'Rscript %s -i %s -x %s -y %s -o %s' % (VisHPD95_R, dm_out, plot_width, plot_height, plot_out)
    os.system(plot_cmd)
    print('Plot exported to: %s' % plot_out)


if __name__ == '__main__':

    # arguments for rename_seq_parser
    VisHPD95_parser = argparse.ArgumentParser()
    VisHPD95_parser.add_argument('-i',      required=True,                      help='mcmc.txt file or folder')
    VisHPD95_parser.add_argument('-n',      required=True,                      help='Nodes to plot')
    VisHPD95_parser.add_argument('-label',  required=False, default=None,       help='labels on y axis')
    VisHPD95_parser.add_argument('-x',      required=False, default=8,type=int, help='plot width, default: 8')
    VisHPD95_parser.add_argument('-y',      required=False, default=5,type=int, help='plot height, default: 5')
    VisHPD95_parser.add_argument('-o',      required=True,                      help='Output plot')
    args = vars(VisHPD95_parser.parse_args())
    args['VisHPD95_R'] = config_dict['VisHPD95_R']
    VisHPD95(args)

'''

cd /Users/songweizhi/Desktop/777
python3 ~/PycharmProjects/TreeSAK/TreeSAK/VisHPD95.py -i M1_mcmc_txt -o M1_HPD95.pdf -n nodes_five.txt -label y_label_out.txt

'''
