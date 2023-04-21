from ete3 import Tree


def combine_trees(t1_with_len, t2_with_name, op_tree_with_both):

    # assume t1 has brancn length
    # assume t2 has internal node name

    t1 = Tree(t1_with_len, format=0)
    t2 = Tree(t2_with_name, format=1)

    t1_leaves_to_node_dict = dict()
    for t1_node in t1.traverse():
        leaf_str = '__'.join(sorted(list(t1_node.get_leaf_names())))
        t1_leaves_to_node_dict[leaf_str] = t1_node

    t2_leaves_to_node_dict = dict()
    for t2_node in t2.traverse():
        leaf_str = '__'.join(sorted(list(t2_node.get_leaf_names())))
        t2_leaves_to_node_dict[leaf_str] = t2_node

    t1_node_to_t2_node_dict = dict()
    for index, t1_node in t1_leaves_to_node_dict.items():
        t2_node = t2_leaves_to_node_dict[index]
        t1_node_to_t2_node_dict[t1_node] = t2_node

    merged_tree = t1.copy()
    for node, t1_node in zip(merged_tree.traverse(), t1.traverse()):
        node.name = t1_node_to_t2_node_dict[t1_node].name
    merged_tree.write(outfile=op_tree_with_both, format=3)


def prefix_internal_nodes(tree_in, prefix_str, tree_out):
    t = Tree(tree_in, format=3)
    t_renamed = t.copy()
    for node in t_renamed.traverse():
        if not node.is_leaf():
            node_name_prefixed = '%s%s' % (prefix_str, node.name)
            node.name = node_name_prefixed
        t_renamed.write(outfile=tree_out, format=3)


# t1_file_with_length             = '/Users/songweizhi/Desktop/test/with_branck_length.treefile'
# t2_file_with_internal_node_name = '/Users/songweizhi/Desktop/test/with_internal_node.tree'
# op_tree_with_both               = '/Users/songweizhi/Desktop/test/with_both.tree'
# op_tree_with_both_renamed       = '/Users/songweizhi/Desktop/test/with_both_renamed.tree'
# interal_node_prefix             = 'IN'

t1_file_with_length             = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_genome_tree_for_ALE.treefile'
t2_file_with_internal_node_name = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_for_ALE.ufboot_ALE_renamed_genome_tree.tree'
op_tree_with_both               = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_with_both.tree'
op_tree_with_both_renamed       = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_with_both_renamed.tree'
interal_node_prefix             = 'IN'


# combine_trees
combine_trees(t1_file_with_length, t2_file_with_internal_node_name, op_tree_with_both)

# prefix_internal_nodes of combined tree
prefix_internal_nodes(op_tree_with_both, interal_node_prefix, op_tree_with_both_renamed)

