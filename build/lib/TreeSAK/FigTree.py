import os
import argparse
from ete3 import Tree


FigTree_usage = '''
====================== FigTree example commands ======================

TreeSAK FigTree -h

======================================================================
'''


def FigTree(args):

    input_txt_file  = args['i']
    op_txt          = args['o']

    if os.path.isfile(input_txt_file) is False:
        print('Metadata file not found, program exited!')
        exit()


if __name__ == '__main__':

    FigTree_parser = argparse.ArgumentParser(usage=FigTree_usage)
    FigTree_parser.add_argument('-i',     required=True,                         help='input metadata')
    FigTree_parser.add_argument('-tree',  required=False, default=None,          help='gene id, in tree file')
    FigTree_parser.add_argument('-txt',   required=False, default=None,          help='gene id, in txt file')
    FigTree_parser.add_argument('-o',     required=True,                         help='output metadata')
    FigTree_parser.add_argument('-na',    required=False, action='store_true',   help='include leaves with na values')
    args = vars(FigTree_parser.parse_args())
    FigTree(args)
