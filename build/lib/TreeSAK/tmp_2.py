import os
import glob
import argparse
from ete3 import Tree
import multiprocessing as mp


gnm_tree_leaf_rename_txt = '/Users/songweizhi/Documents/Research/Sponge/11_ALE_wd/ALE2_op_dir/genome_tree_leaf_rename.txt'
genome_tree_file_rooted  = '/Users/songweizhi/Documents/Research/Sponge/11_ALE_wd/OMA_cov85_213_top25_BMGE.rooted.treefile'


# prepare genome tree for running ALE
gnm_tree_leaf_rename_txt_handle = open(gnm_tree_leaf_rename_txt, 'w')
gnm_tree_in = Tree(genome_tree_file_rooted, format=1)
rename_dict = dict()
for leaf in gnm_tree_in:
    leaf_name = leaf.name
    leaf_name_new = leaf_name.replace('_', '')
    gnm_tree_leaf_rename_txt_handle.write('%s\t%s\n' % (leaf_name_new, leaf.name))
    leaf.name = leaf_name_new
    rename_dict[leaf_name] = leaf_name_new
gnm_tree_leaf_rename_txt_handle.close()

