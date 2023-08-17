import glob
from ete3 import Tree


gnm_tree_file_re = '/Users/songweizhi/Desktop/in_files/*_genome_tree_for_ALE.treefile'


gnm_tree_file_list = glob.glob(gnm_tree_file_re)

for each_tree in gnm_tree_file_list:
    t = Tree(each_tree, format=1)
    leaf_set = set()
    for leaf in t.iter_leaves():
        leaf_set.add(leaf.name)
    print('%s\t%s\t%s' % (each_tree, len(leaf_set), leaf_set))

