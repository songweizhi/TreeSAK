import random
import dendropy
import argparse
from ete3 import Tree


RootTree_usage = '''
====================== RootTree example commands ======================

TreeSAK RootTree -i input.tree -og outgroup_genomes.txt -o rooted.tree

=======================================================================
'''


def root_with_outgroup(input_tree, out_group_list, add_root_branch, tree_file_rooted):

    """
    Reroot the tree using the given outgroup.
    modified based on: https://github.com/Ecogenomics/GTDBTk/blob/master/gtdbtk/reroot_tree.py

    input_tree:  File containing Newick tree to rerooted.
    output_tree: Name of file for rerooted tree.
    outgroup:    Labels of taxa in outgroup.
    """

    tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)

    outgroup_in_tree = set()
    ingroup_leaves = set()
    for n in tree.leaf_node_iter():
        if n.taxon.label in out_group_list:
            outgroup_in_tree.add(n.taxon)
        else:
            ingroup_leaves.add(n)

    # Since finding the MRCA is a rooted tree operation, the tree is first rerooted on an ingroup taxa. This
    # ensures the MRCA of the outgroup can be identified so long as the outgroup is monophyletic. If the
    # outgroup is polyphyletic trying to root on it is ill-defined. To try and pick a "good" root for
    # polyphyletic outgroups, random ingroup taxa are selected until two of them give the same size
    # lineage. This will, likely, be the smallest bipartition possible for the given outgroup though
    # this is not guaranteed.

    mrca = tree.mrca(taxa=outgroup_in_tree)
    mrca_leaves = len(mrca.leaf_nodes())
    while True:
        rnd_ingroup = random.sample(list(ingroup_leaves), 1)[0]
        tree.reroot_at_edge(rnd_ingroup.edge, length1=0.5 * rnd_ingroup.edge_length, length2=0.5 * rnd_ingroup.edge_length)
        mrca = tree.mrca(taxa=outgroup_in_tree)
        if len(mrca.leaf_nodes()) == mrca_leaves:
            break
        mrca_leaves = len(mrca.leaf_nodes())

    if mrca.edge_length is not None:
        tree.reroot_at_edge(mrca.edge, length1=0.5 * mrca.edge_length, length2=0.5 * mrca.edge_length)

        # tree.write_to_path(tree_file_rooted, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        tree_out_string = tree.as_string(schema='newick', suppress_rooting=True, unquoted_underscores=True)
        tree_out_string = tree_out_string.replace("'", "")

        # add the root bar
        if add_root_branch is True:
            tree_out_string = '(' + tree_out_string
            tree_out_string = tree_out_string.replace(');', '):0.02);')

        # write out tree string
        tree_file_rooted_handle = open(tree_file_rooted, 'w')
        tree_file_rooted_handle.write(tree_out_string)
        tree_file_rooted_handle.close()


def RootTree(args):

    tree_file           = args['i']
    out_group_txt       = args['og']
    tree_file_rooted    = args['o']
    tree_fmt            = args['fmt']
    add_root_branch     = args['add_root']

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    # tre = Tree(tree_file, format=tree_fmt)
    # out_group_lca = tre.get_common_ancestor(out_group_set)
    # tre.set_outgroup(out_group_lca)
    # tre.write(outfile=tree_file_rooted, format=tree_fmt)

    root_with_outgroup(tree_file, out_group_set, add_root_branch, tree_file_rooted)


if __name__ == '__main__':

    RootTree_parser = argparse.ArgumentParser()
    RootTree_parser.add_argument('-i',          required=True,                          help='input tree')
    RootTree_parser.add_argument('-og',         required=True,                          help='out group leaves')
    RootTree_parser.add_argument('-o',          required=True,                          help='output tree')
    RootTree_parser.add_argument('-add_root',   required=False, action='store_true',    help='add the root branch')
    RootTree_parser.add_argument('-fmt',        required=False, default=1, type=int,    help='tree format, default: 1')
    args = vars(RootTree_parser.parse_args())
    RootTree(args)
