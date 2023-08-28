# import glob
# from ete3 import Tree
#
#
# tree_file = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs/ALE3_op_dir_c0.75/SpeciesTreeRef.newick'
#
#
# tree_in = Tree(tree_file, format=1)
#
#
# # n = 0
# # for leaf in tree_in:
# #     print(leaf.name)
# #     n += 1
# # print(n)
#
# leaf_node_num = 0
# nonleaf_node_num = 0
# nonleaf_node_list = []
# for node in tree_in.traverse():
#     if node.is_leaf():
#         leaf_node_num += 1
#     else:
#         nonleaf_node_num += 1
#         print(node.name)
#         nonleaf_node_list.append(node.name)
#
# print(leaf_node_num)
# print(nonleaf_node_num)
# print(sorted(nonleaf_node_list))




aaa_txt = '/Users/songweizhi/Desktop/aaa.txt'


p_set = set()
c_set = set()
o_set = set()
f_set = set()
g_set = set()
for each_line in open(aaa_txt):
    each_line_split = each_line.strip().split('\t')
    print(each_line_split)
    p_set.add(each_line_split[0])
    c_set.add(each_line_split[1])
    o_set.add(each_line_split[3])
    f_set.add(each_line_split[4])
    g_set.add(each_line_split[5])

print(len(p_set))
print(len(c_set))
print(len(o_set))
print(len(f_set))
print(len(g_set))







