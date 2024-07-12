import argparse
from ete3 import Tree


mcmcTC_usage = '''
===================== mcmcTC example commands =====================

TreeSAK mcmcTC -i in.tree -o out.tree -tc time_constraints.txt

# Format of constraint file (tab separated columns)
IMG2264867070_yang,GCF900696045_1_yang	3.46-4.38	Archaeal root
GCF000015225_1_yang,GCF000007225_1_yang	-2.32	Oxygen Age Constraint, Thermoproteales
GCF000213215_1_yang,GCA000024305_1_yang	-2.32	Oxygen Age Constraint, Sulfolobales
GCF000152265_2_yang,GCF000195915_1_yang	-2.32	Oxygen Age Constraint, Thermoplasma
GCF000376445_1_yang,GCF000172995_2_yang	-1.579	Chitin Age Constraint, Halobacteriales
GCF000195935_2_yang,GCF000151105_2_yang	-1.579	Chitin Age Constraint, Thermococcales
GCA000802205_2_yang,GCA000200715_1	0.75-1.49	HGT from Viridiplantae to Thaumarchaeota

===================================================================
'''


def mcmcTC(args):

    tree_file_in        = args['i']
    time_constraint_txt = args['tc']
    tree_file_out       = args['o']

    constraint_set = set()
    constraint_dict = dict()
    not_recognizable_time_constraint_set = set()
    for each_constraint in open(time_constraint_txt):
        each_constraint_split = each_constraint.strip().split('\t')
        leaf_ids = each_constraint_split[0]
        provided_age = each_constraint_split[1]

        str_to_add = ''
        if provided_age.startswith('-'):
            str_to_add = '<%s' % provided_age[1:]
        elif provided_age.endswith('-'):
            str_to_add = '>%s' % provided_age[:-1]
        elif '-' in provided_age:
            provided_age_split = provided_age.split('-')
            str_to_add = '>%s<%s' % (provided_age_split[0], provided_age_split[1])
        else:
            not_recognizable_time_constraint_set.add(provided_age)

        constraint_dict[leaf_ids] = str_to_add
        constraint_set.add(str_to_add)

    if len(not_recognizable_time_constraint_set) > 0:
        print('Format of the following time constraints are not recognizable, program exited')
        print(','.join(not_recognizable_time_constraint_set))
        exit()

    # read in tree
    tree_in = Tree(tree_file_in, quoted_node_names=True, format=1)

    # add time constraints as node name
    for each_node in constraint_dict:
        node_age = constraint_dict[each_node]
        node_split = each_node.split(',')
        current_lca = tree_in.get_common_ancestor(node_split)
        current_lca.add_features(custom_label=node_age)
        current_lca.name = node_age

    tree_out_str = tree_in.write(format=1)

    # remove branch length of 1
    tree_out_str = tree_out_str.replace(':1', '')

    # quote constraint strings
    for each_constraint in constraint_set:
        tree_out_str = tree_out_str.replace(each_constraint, ("'%s'" % each_constraint))

    # remove the most outside parenthesis
    tree_out_str = tree_out_str[1:].replace(');', ';')

    # write tree to file
    with open(tree_file_out, 'w') as tree_file_out_handle:
        tree_file_out_handle.write('%s\t1\n' % len(tree_in.get_leaves()))
        tree_file_out_handle.write(tree_out_str)


if __name__ == '__main__':

    mcmcTC_parser = argparse.ArgumentParser()
    mcmcTC_parser.add_argument('-i',   required=True,  help='input tree')
    mcmcTC_parser.add_argument('-o',   required=True,  help='output tree')
    mcmcTC_parser.add_argument('-tc',  required=True,  help='time constraint file')
    args = vars(mcmcTC_parser.parse_args())
    mcmcTC(args)
