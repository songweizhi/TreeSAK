import os
import argparse
from TreeSAK.TreeSAK_config import config_dict


VisHPD95_usage = '''
============= VisHPD95 example command =============

TreeSAK VisHPD95 -i dm.txt -o plot.pdf
TreeSAK VisHPD95 -i dm.txt -o plot.pdf -x 9 -y 6

# format of dm.txt
Test Var  Mean  Low High
Topo_1   Bacteria  1.03    0.95    1.10
Topo_1   Archaea    0.93    0.86    1.01
Topo_2   Bacteria  1.09    1.02    1.16
Topo_2   Archaea    1.04    0.96    1.13

====================================================
'''


def VisHPD95(args):

    data_in              = args['i']
    plot_width           = args['x']
    plot_height          = args['y']
    plot_out             = args['o']
    plot_grouped_HPD95_R = args['plot_grouped_HPD95_R']

    plot_cmd = 'Rscript %s -i %s -x %s -y %s -o %s' % (plot_grouped_HPD95_R, data_in, plot_width, plot_height, plot_out)
    os.system(plot_cmd)
    print('Plot exported to: %s' % plot_out)


if __name__ == '__main__':

    # arguments for rename_seq_parser
    VisHPD95_parser = argparse.ArgumentParser()
    VisHPD95_parser.add_argument('-i',  required=True,                      help='input data matrix')
    VisHPD95_parser.add_argument('-x',  required=False, default=8,type=int, help='plot width, default: 8')
    VisHPD95_parser.add_argument('-y',  required=False, default=5,type=int, help='plot height, default: 5')
    VisHPD95_parser.add_argument('-o',  required=True,                      help='output plot')
    args = vars(VisHPD95_parser.parse_args())
    args['plot_grouped_HPD95_R'] = config_dict['plot_grouped_HPD95_R']
    VisHPD95(args)
