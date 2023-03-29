from ete3 import Tree


def root_with_out_group(tree_file, out_group_txt, tree_file_rooted):

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=1)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted)


def replace_clades(main_tree, sub_tree, tree_out):

    # read in sub tree
    tre_sub = Tree(sub_tree, format=1)

    # get all leaves in sub tree
    subtree_leaf_name_list = tre_sub.get_leaf_names()

    # read in main tree
    tre_main = Tree(main_tree)

    # remove clades
    lca = tre_main.get_common_ancestor(subtree_leaf_name_list)

    if len(lca.get_leaf_names()) != len(subtree_leaf_name_list):
        print('LCA of subtree leaves in main tree contain extra leaves, program exited!')
        exit()

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(tre_sub)

    # write out updated tree
    tre_main.write(outfile=tree_out, format=8)


tree_file                           = '/Users/songweizhi/Desktop/777/PA_75_DeltaLL_50_raw.treefile'
out_group_txt                       = '/Users/songweizhi/Desktop/777/out_group.txt'
tree_file_rooted                    = '/Users/songweizhi/Desktop/777/PA_75_DeltaLL_50_rooted.treefile'
eu_tree                             = '/Users/songweizhi/Desktop/777/27.nwk'
rooted_tree_with_time_constraints   = '/Users/songweizhi/Desktop/777/PA_75_DeltaLL_50_rooted_with_time_constraints.treefile'


root_with_out_group(tree_file, out_group_txt, tree_file_rooted)
replace_clades(tree_file_rooted, eu_tree, rooted_tree_with_time_constraints)


root_with_out_group_from_tianhua = '''
from ete3 import Tree
tpath = ''
nodes = ['F7','B7_3','A7']
tre = Tree(tpath)
tre2 = tre.copy()
lca = tre.get_common_ancestor(nodes)
lca_leaves = lca.get_leaf_names()
# intersect = set(lca_leaves).intersection(set(nodes))
ratio = len(nodes)/len(lca_leaves)
if ratio > 0.5:
    tre2.set_outgroup(lca)  # inplace
    tre2.write(outfile='')
else:
    print(ratio)
'''


replace_clade_from_tianhua ='''
## replace clade with a new tree    
tre2 = Tree(tpath2,format=3)  
nodes1 = []
tre = Tree(tpath)
lca = tre.get_common_ancestor(nodes1)
lca_p = lca.up
lca_p.remove_child(lca)
if len(tre2.children)==1:
    lca_p.add_child(tre2.children[0])
else:
    lca_p.add_child(tre2)
tre.write(outfile='',format='')
'''