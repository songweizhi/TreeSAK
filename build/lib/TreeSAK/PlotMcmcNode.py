import os
import glob
import argparse
import pandas as pd
from ete3 import Tree
import plotly.express as px


PlotMcmcNode_usage = '''
============================ PlotMcmcNode example commands ============================

TreeSAK PlotMcmcNode -i McmcTree_op_files -n n179 -o Clock2_n179.pdf
TreeSAK PlotMcmcNode -i McmcTree_op_files -n n161,n186 -o Clock3_n161_n186.pdf
TreeSAK PlotMcmcNode -i McmcTree_op_files -n nodes.txt -o multi_runs_multi_nodes.pdf

# File name of the mcmc.txt and the corresponding mcmc out file need to follow 
# the rule as specified below:
[prefix]_mcmc.txt
[prefix]_out.txt

# file format (-n, tab separated)
# leave the 2nd column blank for nodes without renaming
node1	Bacteria
node2
node3,node9	Archaea

# Y-axis label file format (-l, tab separated)
PA_75_DeltaLL_50_clock3_mcmc.txt	DeltaLL_50
PA_75_DeltaLL_75_clock3_mcmc.txt	DeltaLL_75

=======================================================================================
'''

def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def mcmctree_out_to_tree_str(mamctree_out):

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

    return tree_str_renamed


def plot_distribution(df_txt, output_plot):

    df = pd.read_table(df_txt, sep=',')
    run_id_list = df['Setting'].unique()
    node_id_list = df['Node'].unique()

    # sort dataframe by run id
    df = df.sort_values(by='Setting', ascending=False)

    plot_width  = 900
    plot_height = len(run_id_list)*100
    if plot_height < 360:
        plot_height = 360

    fig = px.violin(df, x="Value", y="Setting", color="Node", points=False, orientation="h", width=plot_width, height=plot_height)
    if len(node_id_list) == 1:
        fig.update_traces(side="positive", fillcolor='lightblue', width=1.6, opacity=0.75)
    else:
        fig.update_traces(side="positive", fillcolor='rgba(0,0,0,0)', width=1.6)

    fig.update_traces(showlegend=True)
    fig.layout.template = "simple_white"
    # fig.layout.width = 700
    # fig.layout.height = 750
    # fig.update_xaxes(range=[40, 0])
    # fig.update_layout(margin_t=10, title_text='Demo', title_x=0.5)
    fig.write_image(output_plot)


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
                    print('*out.txt file not found, program exited!')
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


def PlotMcmcNode(args):

    mcmc_in             = args['i']
    node_txt            = args['n']
    output_plot         = args['o']
    specified_out_file  = args['of']
    y_label_txt         = args['l']

    # check MCMCTree output file/dir
    if os.path.isfile(mcmc_in) is True:
        mcmc_file_list = [mcmc_in]
    else:
        mcmc_file_re = '%s/*_mcmc.txt' % (mcmc_in)
        mcmc_file_list = glob.glob(mcmc_file_re)

    if len(mcmc_file_list) == 0:
        print('*mcmc.txt file not found, program exited!')
        exit()

    if specified_out_file is None:
        missed_out_file_list = []
        for each_mcmc_file in mcmc_file_list:
            pwd_out_file = each_mcmc_file.replace('_mcmc.txt', '_out.txt')
            if os.path.isfile(pwd_out_file) is False:
                missed_out_file_list.append(pwd_out_file)
        if len(missed_out_file_list) > 0:
            print('The following *out.txt files are missing, program exited!')
            print('\n'.join(sorted(missed_out_file_list)))
            exit()

    # read in y-axis label file
    y_label_dict = dict()
    if y_label_txt is not None:
        for each_sample in open(y_label_txt):
            each_sample_split = each_sample.strip().split('\t')
            if len(each_sample_split) == 2:
                y_label_dict[each_sample_split[0]] = each_sample_split[1]
            else:
                print('Format error: %s' % y_label_txt)
                exit()

    _, f_base, _ = sep_path_basename_ext(output_plot)
    found_matched_node = False
    op_tree_tmp  = '%s_tree.txt'  % f_base
    op_df_tmp    = '%s_data.txt'  % f_base
    op_label_tmp = '%s_label.txt' % f_base

    op_label_tmp_handle = open(op_label_tmp, 'w')
    op_tree_tmp_handle = open(op_tree_tmp, 'w')
    op_df_tmp_handle = open(op_df_tmp, 'w')
    op_df_tmp_handle.write('Value,Node,Setting\n')
    for mcmc_file in mcmc_file_list:

        mcmc_file_no_path = mcmc_file
        if '/' in mcmc_file_no_path:
            mcmc_file_no_path = mcmc_file_no_path.split('/')[-1]

        if specified_out_file is None:
            pwd_current_run_mcmc_out = mcmc_file.replace('_mcmc.txt', '_out.txt')
        else:
            pwd_current_run_mcmc_out = specified_out_file
        node_set, node_rename_dict, tree_str = get_internal_node_to_plot(node_txt, pwd_current_run_mcmc_out)
        op_tree_tmp_handle.write('%s\t%s\n' % (mcmc_file_no_path.replace('_mcmc.txt', ''), tree_str))
        label_to_write = y_label_dict.get(mcmc_file_no_path, mcmc_file_no_path)
        mcmc_df = pd.read_table(mcmc_file, index_col=0)
        for each_col in mcmc_df:
            if each_col in node_set:
                node_name_to_write = node_rename_dict.get(each_col, each_col)
                found_matched_node = True
                value_list = mcmc_df[each_col].values
                for each_value in value_list:
                    op_df_tmp_handle.write('%s,%s,%s\n' % (each_value, node_name_to_write, label_to_write))

                op_label_tmp_handle.write('%s\t%s\t%s\n' % (label_to_write, each_col, node_name_to_write))
    op_df_tmp_handle.close()
    op_label_tmp_handle.close()
    op_tree_tmp_handle.close()

    if found_matched_node is False:
        print('Provided node(s) not found, program exited!')
        exit()

    # plot distribution
    plot_distribution(op_df_tmp, output_plot)

    print('Plot exported to %s, done!' % output_plot)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,                   help='folder holds the *mcmc.txt and *out.txt files')
    parser.add_argument('-of',     required=False, default=None,    help='the *out.txt file')
    parser.add_argument('-n',      required=True,                   help='Nodes to plot')
    parser.add_argument('-l',      required=False, default=None,    help='labels on y axis')
    parser.add_argument('-o',      required=True,                   help='Output plot')
    args = vars(parser.parse_args())
    PlotMcmcNode(args)
