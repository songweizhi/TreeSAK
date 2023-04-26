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

==========================================================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def CompareMCMC(args):

    mcmc_txt_x      = args['mx']
    mcmc_txt_y      = args['my']
    label_x         = args['lx']
    label_y         = args['ly']
    pwd_figure      = args['o']
    max_axis_value  = args['max']
    label_fs        = args['fs']

    x_path, x_basename, x_ext = sep_path_basename_ext(mcmc_txt_x)
    y_path, y_basename, y_ext = sep_path_basename_ext(mcmc_txt_y)

    if label_x is None:
        label_x = x_basename
    if label_y is None:
        label_y = y_basename

    # read in dataframe
    df_x = pd.read_table(mcmc_txt_x, index_col=0)
    df_y = pd.read_table(mcmc_txt_y, index_col=0)

    # get Mean value for each column
    df_x_col_to_mean_dict = {col_name: mean for col_name, mean in df_x.mean().iteritems()}
    df_y_col_to_mean_dict = {col_name: mean for col_name, mean in df_y.mean().iteritems()}

    # get CI95 for each column
    df_x_col_to_ci_dict = {col_name: az.hdi(col.values, hdi_prob=0.95) for col_name, col in df_x.iteritems()}
    df_y_col_to_ci_dict = {col_name: az.hdi(col.values, hdi_prob=0.95) for col_name, col in df_y.iteritems()}

    num_list_x = []
    num_list_y = []
    err_range_x = []
    err_range_y = []
    for col_name, col in df_x.iteritems():
        if col_name not in ['mu', 'sigma2', 'lnL']:
            num_list_x.append(df_x_col_to_mean_dict[col_name])
            num_list_y.append(df_y_col_to_mean_dict[col_name])
            err_range_x.append(df_x_col_to_ci_dict[col_name])
            err_range_y.append(df_y_col_to_ci_dict[col_name])

    x_err_l = []
    x_err_r = []
    y_err_l = []
    y_err_u = []
    max_value = 0
    min_value = 100000000000000
    n = 0
    while n < len(num_list_x):
        x_value = num_list_x[n]
        y_value = num_list_y[n]
        x_range = err_range_x[n]
        y_range = err_range_y[n]
        x_l_dist = abs(x_value - x_range[0])
        x_r_dist = abs(x_range[1] - x_value)
        y_l_dist = abs(y_value - y_range[0])
        y_u_dist = abs(y_range[1] - y_value)
        x_err_l.append(x_l_dist)
        x_err_r.append(x_r_dist)
        y_err_l.append(y_l_dist)
        y_err_u.append(y_u_dist)

        current_max = max(x_value, y_value, x_range[0], x_range[1], y_range[0], y_range[1])
        current_min = min(x_value, y_value, x_range[0], x_range[1], y_range[0], y_range[1])

        if current_max > max_value:
            max_value = current_max
        if current_min < min_value:
            min_value = current_min
        n += 1

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


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-mx',      required=True,                          help='mcmc.txt for x axis')
    parser.add_argument('-my',      required=True,                          help='mcmc.txt for y axis')
    parser.add_argument('-lx',      required=False, default=None,           help='label for x axis')
    parser.add_argument('-ly',      required=False, default=None,           help='label for y axis')
    parser.add_argument('-max',     required=False, default=None, type=int, help='maximum axis value')
    parser.add_argument('-fs',      required=False, default=16, type=int,   help='label font size, default: 16')
    parser.add_argument('-o',       required=True,                          help='output plot')
    args = vars(parser.parse_args())
    CompareMCMC(args)
