from ete3 import Tree


tree_file = "(A:1,(B:1,(C:1,D:1)CD:0.5)BCD:0.5)ABCD;"
node_a = 'BCD'
node_b = 'D'


def check_a_is_ancestor_of_b(tree_file, node_a, node_b):

    a_is_ancestor_of_b = False
    for node in Tree(tree_file, format=1).traverse():
        node_name = node.name
        if node_name == node_b:
            node_ancestor_list = [i.name for i in node.get_ancestors()]
            if node_a in node_ancestor_list:
                a_is_ancestor_of_b = True

    return a_is_ancestor_of_b


def check_a_is_child_of_b(tree_file, node_a, node_b):

    a_is_child_of_b = False
    for node in Tree(tree_file, format=1).traverse():
        node_name = node.name
        if node_name == node_b:
            node_children_list = [i.name for i in node.get_descendants()]
            if node_a in node_children_list:
                a_is_child_of_b = True

    return a_is_child_of_b

a_is_ancestor_of_b = check_a_is_ancestor_of_b(tree_file, 'ABCD', 'CD')



a_is_child_of_b = check_a_is_child_of_b(tree_file, 'BCD', 'ABCD')

print(a_is_ancestor_of_b)
print(a_is_child_of_b)

# # get the nodes of interest
# node_a = t.search_nodes(name="A")[0]
# node_b = t.search_nodes(name="B")[0]
# node_c = t.search_nodes(name="C")[0]
#
# # check if node_a is node_b's parent
# if node_a.is_ancestor(node_b):
#     print("Node A is a parent of Node B.")
#
# # check if node_c is node_b's child
# if node_c.is_descendant(node_b):
#     print("Node C is a child of Node B.")



