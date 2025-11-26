import os
import glob
import argparse
from ete3 import Tree


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


file_dir    = '/Users/songweizhi/Desktop/00'
file_ext    = 'treefile'
tree_format = 1
opdir       = '/Users/songweizhi/Desktop/00_renamed'


file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)

for tree_file in file_list:
    f_name, f_path, f_base, f_ext = sep_path_basename_ext(tree_file)
    tree_out = '%s/%s' % (opdir, f_name)
    t = Tree(tree_file, format=tree_format)
    for leaf in t:
        leaf_name_new = '_'.join(leaf.name.split('_')[:-1])
        leaf.name = leaf_name_new
    t.write(format=tree_format, outfile=tree_out)
