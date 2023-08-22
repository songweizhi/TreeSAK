import glob
from ete3 import Tree


tree_file = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs/ALE3_op_dir_c0.75/SpeciesTreeRef.newick'


tree_in = Tree(tree_file, format=1)


# n = 0
# for leaf in tree_in:
#     print(leaf.name)
#     n += 1
# print(n)

leaf_node_num = 0
nonleaf_node_num = 0
nonleaf_node_list = []
for node in tree_in.traverse():
    if node.is_leaf():
        leaf_node_num += 1
    else:
        nonleaf_node_num += 1
        print(node.name)
        nonleaf_node_list.append(node.name)

print(leaf_node_num)
print(nonleaf_node_num)
print(sorted(nonleaf_node_list))


