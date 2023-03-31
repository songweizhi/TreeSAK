import os
import glob
import argparse
import pandas as pd
import plotly.express as px


PlotMcmcNode_usage = '''
======================= PlotMcmcNode example commands =======================

TreeSAK PlotMcmcNode -i Clock2.txt -n t_n179 -o Clock2_t_n179.pdf
TreeSAK PlotMcmcNode -i MCMC_op -x txt -n nodes.txt -o MCMC_op_nodes.pdf

=============================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def plot_distribution(df_txt, output_plot):

    df = pd.read_table(df_txt, sep=',')
    fig = px.violin(df, y="Setting", x="Value", color="Node", points=False, orientation="h")
    fig.update_traces(side="positive", fillcolor='rgba(0,0,0,0)', width=1.8)
    fig.update_traces(showlegend=True)
    fig.layout.template = "simple_white"
    # fig.layout.width = 700
    # fig.layout.height = 750
    # fig.update_xaxes(range=[40, 0])
    # fig.update_layout(margin_t=10, title_text='Demo', title_x=0.5)
    fig.write_image(output_plot)


def PlotMcmcNode(args):

    mcmc_dir    = args['i']
    mcmc_txt    = args['x']
    node_txt    = args['n']
    output_plot = args['o']

    # get nodes to plot
    node_set = set()
    if os.path.isfile(node_txt) is True:
        for each in open(node_txt):
            node_set.add(each.strip())
    else:
        node_set.add(node_txt)

    # check MCMCTree output file/dir
    if os.path.isfile(mcmc_dir) is True:
        mcmc_file_list = [mcmc_dir]
    else:
        mcmc_file_re = '%s/*.%s' % (mcmc_dir, mcmc_txt)
        mcmc_file_list = glob.glob(mcmc_file_re)

    if len(mcmc_file_list) == 0:
        print('MCMCTree output file not found, program exited!')
        exit()

    found_matched_node = False
    op_df_tmp = output_plot + '.txt'
    op_df_tmphandle = open(op_df_tmp, 'w')
    op_df_tmphandle.write('Value,Node,Setting\n')
    for mcmc_txt in mcmc_file_list:
        file_path, file_basename, file_extension = sep_path_basename_ext(mcmc_txt)
        mcmc_df = pd.read_table(mcmc_txt, index_col=0)
        for each_col in mcmc_df:
            if each_col in node_set:
                found_matched_node = True
                value_list = mcmc_df[each_col].values
                for each_value in value_list:
                    op_df_tmphandle.write('%s,%s,%s\n' % (each_value, each_col, file_basename))
    op_df_tmphandle.close()

    if found_matched_node is False:
        print('Provided node(s) not found, program exited!')
        exit()

    plot_distribution(op_df_tmp, output_plot)

    print('Plot exported to %s, done!' % output_plot)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,                   help='MCMCTree output file or folder')
    parser.add_argument('-x',      required=False, default='txt',   help='Extension of MCMCTree output file')
    parser.add_argument('-n',      required=True,                   help='Nodes to plot')
    parser.add_argument('-o',      required=True,                   help='Output plot')
    args = vars(parser.parse_args())
    PlotMcmcNode(args)
