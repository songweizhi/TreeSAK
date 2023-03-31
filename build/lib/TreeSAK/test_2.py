import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def PlotMCMC(num_list_x, num_list_y, err_range_x, err_range_y, pwd_figure):

    x_err_l = []
    x_err_r = []
    y_err_l = []
    y_err_u = []
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
        n += 1

    plt.scatter(num_list_x, num_list_y, s=0)
    plt.errorbar(num_list_x, num_list_y, xerr=[x_err_l, x_err_r], yerr=[y_err_l, y_err_u], ls='none', ecolor='black', elinewidth=1)
    plt.tight_layout()
    plt.savefig(pwd_figure)
    plt.close()


num_list_x = [1, 2, 3, 4, 5]
num_list_y = [1, 2, 3, 4, 5]
err_range_x = [[0.8, 1.1], [1.5, 2.1], [2.9, 3.7], [3.2, 4.1], [4.5, 5.5]]
err_range_y = [[0.7, 1.3], [1.9, 2.5], [2.9, 3.7], [3.1, 4.4], [4.5, 5.5]]
pwd_figure = '/Users/songweizhi/Desktop/aaa.png'

PlotMCMC(num_list_x, num_list_y, err_range_x, err_range_y, pwd_figure)