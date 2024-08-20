import os
import glob
import argparse
import pandas as pd
import plotly.express as px


PlotMcmcNode_usage = '''
========================== PlotMcmcNode example commands ==========================

TreeSAK PlotMcmcNode -i McmcTree_op_files -n topo123_nodes.txt -o topo123_age.pdf

# file format (-n, tab separated, no header)
# column 1: file name, no extension
# column 2: node id
# column 3: file name on the figure
# column 4: node name on the figure
topo1_clock3_mcmc	t_n171	Topo1	Symbiosis_event_1
topo1_clock3_mcmc	t_n151	Topo1	Symbiosis_event_2
topo1_clock3_mcmc	t_n131	Topo1	Symbiosis_event_3
topo2_clock3_mcmc	t_n171	Topo2	Symbiosis_event_1
topo2_clock3_mcmc	t_n171	Topo3	Symbiosis_event_1

===================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


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


def PlotMcmcNode(args):

    mcmc_in         = args['i']
    mcmc_file_ext   = args['x']
    node_txt        = args['n']
    output_plot     = args['o']

    # check MCMCTree output file/dir
    if os.path.isfile(mcmc_in) is True:
        mcmc_file_list = [mcmc_in]
    else:
        mcmc_file_re = '%s/*.%s' % (mcmc_in, mcmc_file_ext)
        mcmc_file_list = glob.glob(mcmc_file_re)

    if len(mcmc_file_list) == 0:
        print('mcmc file not found, program exited!')
        exit()

    if output_plot[-4:] not in ['.PDF', '.pdf']:
        output_plot = output_plot + '.pdf'
    op_df_tmp = output_plot + '.txt'

    mcmc_file_full_path_dict = dict()
    for each_mcmc_file in mcmc_file_list:
        mcmc_path, mcmc_base, mcmc_ext = sep_path_basename_ext(each_mcmc_file)
        mcmc_file_full_path_dict[mcmc_base] = each_mcmc_file

    # read in nodes to plot
    file_to_node_dict = dict()
    file_rename_dict = dict()
    node_rename_dict = dict()
    not_found_file_set = set()
    for each_node in open(node_txt):
        each_node_split = each_node.strip().split()
        file_name = each_node_split[0]
        if '/' in file_name:
            file_name = file_name.split('/')[-1]

        if file_name not in mcmc_file_full_path_dict:
            not_found_file_set.add(file_name)

        node_id = each_node_split[1]
        node_id_with_file_name = '%s_____%s' % (file_name, node_id)

        if file_name not in file_to_node_dict:
            file_to_node_dict[file_name] = set()
        file_to_node_dict[file_name].add(node_id)

        if len(each_node_split) >= 3:
            file_short_name = each_node_split[2]
            file_rename_dict[file_name] = file_short_name

        if len(each_node_split) == 4:
            node_desc = each_node_split[3]
            node_rename_dict[node_id_with_file_name] = node_desc

    # report not found files and exit
    if len(not_found_file_set) > 0:
        print('MCMC file not found, program exited!')
        exit()

    found_matched_node = False
    op_df_tmp_handle = open(op_df_tmp, 'w')
    op_df_tmp_handle.write('Value,Node,Setting\n')
    for mcmc_file in file_to_node_dict:
        pwd_mcmc_file = mcmc_file_full_path_dict[mcmc_file]
        current_node_set = file_to_node_dict.get(mcmc_file, set())
        file_name_to_plot = file_rename_dict.get(mcmc_file, mcmc_file)
        mcmc_df = pd.read_table(pwd_mcmc_file, index_col=0)
        for each_col in mcmc_df:
            if each_col in current_node_set:
                node_desc_to_plot = node_rename_dict.get(('%s_____%s' % (mcmc_file, each_col)), each_col)
                found_matched_node = True
                value_list = mcmc_df[each_col].values
                for each_value in value_list:
                    op_df_tmp_handle.write('%s,%s,%s\n' % (each_value, node_desc_to_plot, file_name_to_plot))
    op_df_tmp_handle.close()

    if found_matched_node is False:
        print('Provided node(s) not found, program exited!')
        exit()

    # plot distribution
    plot_distribution(op_df_tmp, output_plot)

    # remove tmp files
    os.system('rm %s' % op_df_tmp)

    # final report
    print('Plot exported to %s, done!' % output_plot)


if __name__ == '__main__':

    PlotMcmcNode_parser = argparse.ArgumentParser()
    PlotMcmcNode_parser.add_argument('-i',  required=True,                  help='folder holds the *mcmc.txt files')
    PlotMcmcNode_parser.add_argument('-x',  required=False, default='txt',  help='file extension, default txt')
    PlotMcmcNode_parser.add_argument('-n',  required=True,                  help='nodes to plot')
    PlotMcmcNode_parser.add_argument('-o',  required=True,                  help='output plot')
    args = vars(PlotMcmcNode_parser.parse_args())
    PlotMcmcNode(args)
