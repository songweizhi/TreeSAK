import os
import argparse
from ete3 import Tree


rm_leaf_usage = '''
=============== rm_leaf example commands ===============

TreeSAK rm_leaf -i in.tree -o out.tree -l leaves.txt
TreeSAK rm_leaf -i in.tree -o out.tree -l leaf1
TreeSAK rm_leaf -i in.tree -o out.tree -l leaf1,leaf2

# Tree formats
0	flexible with support values
1	flexible with internal node names
2	all branches + leaf names + internal supports
3	all branches + all names
4	leaf branches + leaf names
5	internal and leaf branches + leaf names
6	internal branches + leaf names
7	leaf branches + all names
8	all names
9	leaf names
100	topology only

========================================================
'''


def rm_leaf(args):

    input_tree_file   = args['i']
    input_tree_fmt    = args['if']
    output_tree_file  = args['o']
    output_tree_fmt   = args['of']
    leaf_txt          = args['l']

    if os.path.isfile(input_tree_file) is False:
        print('%s does not exist, program exited!' % input_tree_file)
        exit()

    # read in leaves to remove
    leaves_to_remove_list = []
    if os.path.isfile(leaf_txt) is True:
        for each_leaf in open(leaf_txt):
            leaves_to_remove_list.append(each_leaf.strip())
    else:
        leaves_to_remove_list = leaf_txt.split(',')

    if len(leaves_to_remove_list) == 0:
        print('Please provide at least one leaf to remove, program exited!')
        exit()

    # read in tree and remobe leaf/leaves
    main_tre = Tree(input_tree_file, quoted_node_names=True, format=input_tree_fmt)
    main_tre_leaves = [i.name for i in main_tre.get_leaves()]

    # get leaves to keep
    leaves_to_keep = []
    for each_leaf in main_tre_leaves:
        if each_leaf not in leaves_to_remove_list:
            leaves_to_keep.append(each_leaf)

    if len(main_tre_leaves) == len(leaves_to_keep):
        print('No leaves to remove, program exited!')
        exit()

    main_tre.prune(leaves_to_keep)

    # write out updated tree
    main_tre.write(outfile=output_tree_file, format=output_tree_fmt)


if __name__ == '__main__':

    rm_leaf_parser = argparse.ArgumentParser()
    rm_leaf_parser.add_argument('-i',   required=True,                        help='input tree file')
    rm_leaf_parser.add_argument('-if',  required=False, default=0, type=int,  help='input tree format, default is 0')
    rm_leaf_parser.add_argument('-l',   required=True,                        help='leaf/leaves to remove')
    rm_leaf_parser.add_argument('-o',   required=True,                        help='output tree file')
    rm_leaf_parser.add_argument('-of',  required=False, default=0, type=int,  help='output tree format, default is 0')
    args = vars(rm_leaf_parser.parse_args())
    rm_leaf(args)
