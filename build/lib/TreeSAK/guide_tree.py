import os
import argparse
from ete3 import Tree


guide_tree_usage = '''
======================== guide_tree example command ========================

TreeSAK guide_tree -i in.tree -g gnm_group.txt -o out.tree

============================================================================
'''


def guide_tree(args):

    tree_in         = args['i']
    tree_out        = args['o']
    gnm_group_txt   = args['g']

    group_to_member_dict = dict()
    for each in open(gnm_group_txt):
        each_split = each.strip().split('\t')
        gnm_id  = each_split[0]
        gnm_grp = each_split[1]
        if gnm_grp not in group_to_member_dict:
            group_to_member_dict[gnm_grp] = set()
        group_to_member_dict[gnm_grp].add(gnm_id)

    t = Tree(tree_in, format=0)
    for leaf in t:
        leaf_name = leaf.name
        if leaf_name in group_to_member_dict:
            leaf_member = group_to_member_dict[leaf_name]
            leaf_member_str = ','.join(leaf_member)
            if len(leaf_member) >= 2:
                leaf_member_str = '(' + leaf_member_str + ')'
            if len(leaf_member) == 1:
                leaf.name = leaf_member_str
            else:
                leaf_p = leaf.up
                leaf_p.add_child(Tree((leaf_member_str + ';'), format=0))
                leaf_p.remove_child(leaf)
    t.write(outfile=tree_out, format=9)


if __name__ == '__main__':

    guide_tree_parser = argparse.ArgumentParser(usage=guide_tree_usage)
    guide_tree_parser.add_argument('-i',    required=True,  help='input tree')
    guide_tree_parser.add_argument('-g',    required=True,  help='genome group')
    guide_tree_parser.add_argument('-o',    required=True,  help='output tree')
    args = vars(guide_tree_parser.parse_args())
    guide_tree(args)

