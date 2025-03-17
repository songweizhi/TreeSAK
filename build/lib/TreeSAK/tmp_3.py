import os
import argparse
import arviz as az
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


CompareMCMC_usage = '''
====================================== CompareMCMC example commands ======================================

TreeSAK CompareMCMC -mx IR_mcmc.txt -my AR_mcmc.txt -lx IR -ly AR -o convergence_plot.png -max 40 -fs 12

cd /Users/songweizhi/Desktop
TreeSAK CompareMCMC -mx /Users/songweizhi/Desktop/Sponge_r220/6_dating/MCMCTree/dating_outputs/topo2p10_clock3_nsample250000_run1_mcmc.txt -my /Users/songweizhi/Desktop/Sponge_r220/6_dating/MCMCTree/dating_outputs/topo2p10_clock3_nsample250000_run2_mcmc.txt -lx IR -ly AR -o convergence_plot.png -max 40 -fs 12

==========================================================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def CompareMCMC():

    # mcmc_txt_x      = args['mx']
    # mcmc_txt_y      = args['my']
    # label_x         = args['lx']
    # label_y         = args['ly']
    # pwd_figure      = args['o']
    # max_axis_value  = args['max']
    # label_fs        = args['fs']

    label_fs = 16
    file_x = '/Users/songweizhi/Desktop/x.txt'
    file_y = '/Users/songweizhi/Desktop/y.txt'
    pwd_figure = '/Users/songweizhi/Desktop/Figures.pdf'
    min_value = 0
    max_value = 1
    max_axis_value = 1
    label_x = []
    num_list_x = []
    x_err_l = []
    x_err_r = []
    line_num_index = 0
    for line in open(file_x):
        line_split = line.strip().split('\t')
        if line_num_index > 0:
            label_x.append(line_split[0])
            num_list_x.append(float(line_split[1]))
            x_err_l.append(float(line_split[1]) - float(line_split[2]))
            x_err_r.append(float(line_split[3]) - float(line_split[1]))
        line_num_index += 1


    label_y = []
    num_list_y = []
    y_err_l = []
    y_err_u = []
    line_num_index = 0
    for line in open(file_y):
        line_split = line.strip().split('\t')
        if line_num_index > 0:
            label_y.append(line_split[0])
            num_list_y.append(float(line_split[1]))
            y_err_l.append(float(line_split[1]) - float(line_split[2]))
            y_err_u.append(float(line_split[3]) - float(line_split[1]))
        line_num_index += 1






    figure(figsize=(6, 6), dpi=300)
    plt.plot([min_value, max_value], [min_value, max_value], color='black', linestyle='dashed', linewidth=1, alpha=0.5)
    plt.scatter(num_list_x, num_list_y, s=0)
    plt.errorbar(num_list_x, num_list_y, xerr=[x_err_l, x_err_r], yerr=[y_err_l, y_err_u],
                 ls='none', ecolor='skyblue', elinewidth=1, alpha=0.5)

    if max_axis_value is not None:
        plt.xlim([0, max_axis_value])
        plt.ylim([0, max_axis_value])

    # Set the font size of xticks and yticks
    plt.xticks(fontsize=label_fs)
    plt.yticks(fontsize=label_fs)
    plt.xlabel(label_x, fontsize=label_fs)
    plt.ylabel(label_y, fontsize=label_fs)

    # write out
    plt.tight_layout()
    plt.savefig(pwd_figure)
    plt.close()

    print('Plot exported to %s, done!' % pwd_figure)


CompareMCMC()

# if __name__ == '__main__':
#
#     # initialize the options parser
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-mx',      required=True,                          help='mcmc.txt for x axis')
#     parser.add_argument('-my',      required=True,                          help='mcmc.txt for y axis')
#     parser.add_argument('-lx',      required=False, default=None,           help='label for x axis')
#     parser.add_argument('-ly',      required=False, default=None,           help='label for y axis')
#     parser.add_argument('-max',     required=False, default=None, type=int, help='maximum axis value')
#     parser.add_argument('-fs',      required=False, default=16, type=int,   help='label font size, default: 16')
#     parser.add_argument('-o',       required=True,                          help='output plot')
#     args = vars(parser.parse_args())
#     CompareMCMC(args)
