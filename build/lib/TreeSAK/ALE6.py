import os
import argparse
from Bio import SeqIO
from ete3 import Tree


ALE6_usage = '''
====================================== ALE6 example commands ======================================

# This module is developed to faa ancestral genomes based on ALE outputs
TreeSAK ALE6 -i1 ALE1_op_dir -i3 ALE3_op_dir_70 -s species_tree.rooted.treefile -o ALE6_op_dir_70 -n 380
TreeSAK ALE6 -i1 ALE1_op_dir -i3 ALE3_op_dir_70 -s species_tree.rooted.treefile -o ALE6_op_dir_70 -n 294,309,380,404
TreeSAK ALE6 -i1 ALE1_op_dir -i3 ALE3_op_dir_70 -s species_tree.rooted.treefile -o ALE6_op_dir_70 -n interested_nodes.txt

===================================================================================================
'''


def get_internal_node_leaves(ale_species_tree_file, internal_node_id):

    ale_species_tree = Tree(ale_species_tree_file, format=1)
    internal_node = ale_species_tree.search_nodes(name=internal_node_id)[0]
    internal_node_leaf_object = internal_node.get_leaves()
    internal_node_leaf_set = set()
    for each_leaf in internal_node_leaf_object:
        internal_node_leaf_set.add(each_leaf.name)

    return internal_node_leaf_set


def ALE6(args):

    ale1_op_dir               = args['i1']
    ale3_op_dir               = args['i3']
    op_dir                    = args['o']
    genome_tree_file_rooted   = args['s']
    force_create_op_dir       = args['f']
    interested_internal_nodes = args['n']

    GeneContent_txt = '%s/GeneContent.txt'          % ale3_op_dir
    SpeciesTreeRef  = '%s/SpeciesTreeRef.newick'    % ale3_op_dir

    ########################################## get the id of nodes to process ##########################################

    overall_internal_node_set = set()
    line_num_index = 0
    for each_line in open(GeneContent_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index > 1:
            node_id = line_split[0]
            if '(' not in node_id:
                overall_internal_node_set.add(node_id)

    internal_nodes_to_process = set()
    if interested_internal_nodes is None:
        internal_nodes_to_process = overall_internal_node_set
    else:
        if os.path.isfile(interested_internal_nodes) is False:
            if ',' in interested_internal_nodes:
                internal_nodes_to_process = interested_internal_nodes.split(',')
            else:
                internal_nodes_to_process.add(interested_internal_nodes)
        else:
            for each_node in open(interested_internal_nodes):
                internal_nodes_to_process.add(each_node.strip())

    ####################################################################################################################

    gnm_name_dict_new_to_old = dict()
    for leaf in Tree(genome_tree_file_rooted, format=1):
        leaf_name = leaf.name
        leaf_name_new = leaf_name.replace('_', '')
        leaf.name = leaf_name_new
        gnm_name_dict_new_to_old[leaf_name_new] = leaf_name

    # create output directory
    if force_create_op_dir is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
    os.system('mkdir %s' % op_dir)

    branch_to_leaf_dict = dict()
    branch_to_content_dict = dict()
    col_header_list = []
    line_num_index = 0
    for each_line in open(GeneContent_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_header_list = line_split
        else:
            branch_id = line_split[0]
            if branch_id in internal_nodes_to_process:
                branch_to_content_dict[branch_id] = []
                branch_child_leaf_set = get_internal_node_leaves(SpeciesTreeRef, branch_id)
                branch_to_leaf_dict[branch_id] = branch_child_leaf_set
                for (id, pa) in zip(col_header_list[1:], line_split[1:]):
                    if pa != '0':
                        branch_to_content_dict[branch_id].append(id)

    for each_branch in branch_to_content_dict:
        branch_faa = '%s/%s.faa' % (op_dir, each_branch)
        branch_content = branch_to_content_dict[each_branch]
        branch_child_set = branch_to_leaf_dict[each_branch]
        branch_child_set_original_name = {gnm_name_dict_new_to_old[i] for i in branch_child_set}
        branch_faa_handle = open(branch_faa, 'w')
        for each_prot_family in branch_content:
            each_prot_family_faa = '%s/%s.faa' % (ale1_op_dir, each_prot_family)
            for each_seq in SeqIO.parse(each_prot_family_faa, 'fasta'):
                seq_gnm = '_'.join(each_seq.id.split('_')[:-1])
                if seq_gnm in branch_child_set_original_name:
                    branch_faa_handle.write('>%s %s\n' % (each_seq.id, each_prot_family))
                    branch_faa_handle.write('%s\n' % each_seq.seq)
        branch_faa_handle.close()


if __name__ == '__main__':

    ALE6_parser = argparse.ArgumentParser()
    ALE6_parser.add_argument('-i1', required=True,                          help='ALE1 output directory')
    ALE6_parser.add_argument('-i3', required=True,                          help='ALE3 output directory')
    ALE6_parser.add_argument('-s',  required=True,                          help='rooted species tree')
    ALE6_parser.add_argument('-n',  required=False, default=None,           help='interested internal node(s)')
    ALE6_parser.add_argument('-o',  required=True,                          help='output directory')
    ALE6_parser.add_argument('-f',  required=False, action="store_true",    help='force overwrite')
    args = vars(ALE6_parser.parse_args())
    ALE6(args)
