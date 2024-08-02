import os
import glob
import argparse
from Bio import SeqIO
from ete3 import Tree


ALE6_usage = '''
====================================== ALE6 example commands ======================================

# This module is developed to faa ancestral genomes based on ALE outputs
TreeSAK ALE6 -1 ALE1_op_dir -3 ALE3_op_dir_30 -s species_tree.rooted.treefile -o ALE6_op_dir_30 -n 380 -cog BioSAK_arCOG_wd -kegg BioSAK_KEGG_wd
TreeSAK ALE6 -1 ALE1_op_dir -3 ALE3_op_dir_30 -s species_tree.rooted.treefile -o ALE6_op_dir_30 -n 294,309,380,404
TreeSAK ALE6 -1 ALE1_op_dir -3 ALE3_op_dir_30 -s species_tree.rooted.treefile -o ALE6_op_dir_30 -n interested_nodes.txt

# Needed input files:
-1: faa files
-3: GeneContent.txt and SpeciesTreeRef.newick
-s: the tree used as input for ALE2

# To be added:
1. A dereplication step to the produced faa file.

===================================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_internal_node_leaves(ale_species_tree_file, internal_node_id):

    ale_species_tree = Tree(ale_species_tree_file, format=1)
    internal_node = ale_species_tree.search_nodes(name=internal_node_id)[0]
    internal_node_leaf_object = internal_node.get_leaves()
    internal_node_leaf_set = set()
    for each_leaf in internal_node_leaf_object:
        internal_node_leaf_set.add(each_leaf.name)

    return internal_node_leaf_set


def ALE6(args):

    ale1_op_dir                 = args['1']
    ale3_op_dir                 = args['3']
    op_dir                      = args['o']
    genome_tree_file_rooted     = args['s']
    force_create_op_dir         = args['f']
    interested_internal_nodes   = args['n']
    cog_annotation_wd           = args['cog']
    kegg_annotation_wd          = args['kegg']

    GeneContent_txt         = '%s/GeneContent.txt'          % ale3_op_dir
    SpeciesTreeRef          = '%s/SpeciesTreeRef.newick'    % ale3_op_dir
    transfer_propensity_txt = '%s/Transfer_propensity.txt'  % ale3_op_dir

    ################################################## define op files #################################################

    faa_dir                     = '%s/faa_files'                                    % op_dir
    cog_dir                     = '%s/annotation_COG'                               % op_dir
    kegg_dir                    = '%s/annotation_KEGG'                              % op_dir
    cog_df_txt                  = '%s/annotation_COG.txt'                           % op_dir
    cog_df_desc_txt             = '%s/annotation_COG_desc.txt'                      % op_dir
    kegg_df_txt                 = '%s/annotation_KEGG.txt'                          % op_dir
    kegg_df_desc_txt            = '%s/annotation_KEGG_desc.txt'                     % op_dir
    fun_transfer_propensity_txt = '%s/function_transfer_propensity_weighted.txt'    % op_dir

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
    os.system('mkdir %s' % faa_dir)

    if cog_annotation_wd is not None:
        os.system('mkdir %s' % cog_dir)

    if kegg_annotation_wd is not None:
        os.system('mkdir %s' % kegg_dir)

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

    branch_to_gene_dict = dict()
    for each_branch in branch_to_content_dict:
        branch_faa = '%s/%s.faa' % (faa_dir, each_branch)
        branch_content = branch_to_content_dict[each_branch]
        branch_child_set = branch_to_leaf_dict[each_branch]
        branch_child_set_original_name = {gnm_name_dict_new_to_old[i] for i in branch_child_set}
        branch_faa_handle = open(branch_faa, 'w')
        branch_to_gene_dict[each_branch] = set()
        for each_prot_family in branch_content:
            each_prot_family_faa = '%s/%s.faa' % (ale1_op_dir, each_prot_family)
            for each_seq in SeqIO.parse(each_prot_family_faa, 'fasta'):
                seq_id = each_seq.id
                seq_gnm = '_'.join(seq_id.split('_')[:-1])
                if seq_gnm in branch_child_set_original_name:
                    branch_faa_handle.write('>%s %s\n' % (each_seq.id, each_prot_family))
                    branch_faa_handle.write('%s\n' % each_seq.seq)
                    branch_to_gene_dict[each_branch].add(each_seq.id)
        branch_faa_handle.close()

    # Read in COG annotation results
    fun_to_gene_dict = dict()
    annotation_dict_cog = dict()
    fun_id_to_desc_dict = dict()
    if cog_annotation_wd is not None:

        print('Reading in COG annotation results')
        file_re = '%s/*COG_wd/*_query_to_cog.txt' % (cog_annotation_wd)
        file_list = glob.glob(file_re)

        if len(file_list) == 0:
            print('COG annotation file not detected, program exited!')
            exit()

        for each_file in file_list:
            gnm_id = each_file.split('/')[-1].split('_query_to_cog')[0]
            if gnm_id not in annotation_dict_cog:
                annotation_dict_cog[gnm_id] = dict()
            line_index = 0
            for each_line in open(each_file):
                if line_index > 0:
                    each_line_split = each_line.strip().split('\t')
                    if len(each_line_split) == 4:
                        gene_id = each_line_split[0]
                        cog_id = each_line_split[1]
                        cog_desc = each_line_split[3]
                        annotation_dict_cog[gnm_id][gene_id] = cog_id
                        fun_id_to_desc_dict[cog_id] = cog_desc
                        if cog_id not in fun_to_gene_dict:
                            fun_to_gene_dict[cog_id] = set()
                        fun_to_gene_dict[cog_id].add(gene_id)
                line_index += 1

    # Read in KEGG annotation results
    annotation_dict_kegg = dict()
    if kegg_annotation_wd is not None:

        print('Reading in KEGG annotation results')
        file_re = '%s/*KEGG_wd/*_ko_assignment_ABCD.txt' % (kegg_annotation_wd)
        file_list = glob.glob(file_re)

        if len(file_list) == 0:
            print('KEGG annotation file not detected, program exited!')
            exit()

        for each_file in file_list:
            gnm_id = each_file.split('/')[-1].split('_ko_assignment_ABCD')[0]
            if gnm_id not in annotation_dict_kegg:
                annotation_dict_kegg[gnm_id] = dict()

            line_index = 0
            for each_line in open(each_file):
                if line_index > 0:
                    each_line_split = each_line.strip().split('\t')
                    if len(each_line_split) == 9:
                        gene_id = each_line_split[0]
                        ko_d_id = each_line_split[4][2:]
                        ko_d_desc = each_line_split[8]
                        annotation_dict_kegg[gnm_id][gene_id] = ko_d_id
                        fun_id_to_desc_dict[ko_d_id] = ko_d_desc
                        if ko_d_id not in fun_to_gene_dict:
                            fun_to_gene_dict[ko_d_id] = set()
                        fun_to_gene_dict[ko_d_id].add(gene_id)
                line_index += 1

    cog_dod = dict()
    kegg_dod = dict()
    all_identified_cog_set = set()
    all_identified_kegg_set = set()
    if (cog_annotation_wd is not None) or (kegg_annotation_wd is not None):
        for each_branch in branch_to_gene_dict:
            branch_gene_content = branch_to_gene_dict[each_branch]
            branch_cog_set = set()
            branch_kegg_set = set()
            for each_gene in branch_gene_content:
                gnm_id = '_'.join(each_gene.split('_')[:-1])
                cog_fun = annotation_dict_cog[gnm_id].get(each_gene, 'na')
                kegg_fun = annotation_dict_kegg[gnm_id].get(each_gene, 'na')
                if cog_fun != 'na':
                    branch_cog_set.add(cog_fun)
                    all_identified_cog_set.add(cog_fun)
                if kegg_fun != 'na':
                    branch_kegg_set.add(kegg_fun)
                    all_identified_kegg_set.add(kegg_fun)

            cog_dod[each_branch] = branch_cog_set
            kegg_dod[each_branch] = branch_kegg_set

            # write out annotation
            if len(branch_cog_set) > 0:
                cog_annotation_txt = '%s/%s_COG.txt' % (cog_dir, each_branch)
                cog_annotation_txt_handle = open(cog_annotation_txt, 'w')
                for each_cog in sorted(list(branch_cog_set)):
                    cog_annotation_txt_handle.write('%s\t%s\n' % (each_cog, fun_id_to_desc_dict[each_cog]))
                cog_annotation_txt_handle.close()

            if len(branch_kegg_set) > 0:
                kegg_annotation_txt = '%s/%s_KEGG.txt' % (kegg_dir, each_branch)
                kegg_annotation_txt_handle = open(kegg_annotation_txt, 'w')
                for each_kegg in sorted(list(branch_kegg_set)):
                    kegg_annotation_txt_handle.write('%s\t%s\n' % (each_kegg, fun_id_to_desc_dict[each_kegg]))
                kegg_annotation_txt_handle.close()

    all_identified_cog_list_sorted       = sorted(list(all_identified_cog_set))
    all_identified_kegg_list_sorted      = sorted(list(all_identified_kegg_set))
    all_identified_cog_list_sorted_desc  = [('%s__%s' % (i, fun_id_to_desc_dict[i])) for i in all_identified_cog_list_sorted]
    all_identified_kegg_list_sorted_desc = [('%s__%s' % (i, fun_id_to_desc_dict[i])) for i in all_identified_kegg_list_sorted]

    # write out COG dataframe
    if len(all_identified_cog_set) > 0:
        cog_df_txt_handle = open(cog_df_txt, 'w')
        cog_df_txt_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list_sorted))
        cog_df_desc_txt_handle = open(cog_df_desc_txt, 'w')
        cog_df_desc_txt_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list_sorted_desc))
        for each_branch in sorted(list(cog_dod.keys())):
            branch_cogs = cog_dod[each_branch]
            cog_pa_list = [each_branch]
            for each_cog in all_identified_cog_list_sorted:
                if each_cog in branch_cogs:
                    cog_pa_list.append('1')
                else:
                    cog_pa_list.append('0')
            cog_df_txt_handle.write('\t'.join(cog_pa_list) + '\n')
            cog_df_desc_txt_handle.write('\t'.join(cog_pa_list) + '\n')
        cog_df_txt_handle.close()
        cog_df_desc_txt_handle.close()
        print('Annotation matrix exported to: %s' % cog_df_txt)

    # write out KEGG dataframe
    if len(all_identified_kegg_set) > 0:
        kegg_df_txt_handle = open(kegg_df_txt, 'w')
        kegg_df_txt_handle.write('\t%s\n' % '\t'.join(all_identified_kegg_list_sorted))
        kegg_df_desc_txt_handle = open(kegg_df_desc_txt, 'w')
        kegg_df_desc_txt_handle.write('\t%s\n' % '\t'.join(all_identified_kegg_list_sorted_desc))
        for each_branch in sorted(list(kegg_dod.keys())):
            branch_keggs = kegg_dod[each_branch]
            kegg_pa_list = [each_branch]
            for each_kegg in all_identified_kegg_list_sorted:
                if each_kegg in branch_keggs:
                    kegg_pa_list.append('1')
                else:
                    kegg_pa_list.append('0')
            kegg_df_txt_handle.write('\t'.join(kegg_pa_list) + '\n')
            kegg_df_desc_txt_handle.write('\t'.join(kegg_pa_list) + '\n')
        kegg_df_txt_handle.close()
        kegg_df_desc_txt_handle.close()
        print('Annotation matrix exported to: %s' % kegg_df_txt)

    ################################# get transfer propensity of individual functions ##################################

    print('Getting transfer propensity of individual function')

    # get gene_to_oma_dict
    faa_file_re = '%s/*.faa' % ale1_op_dir
    faa_file_list = glob.glob(faa_file_re)
    gene_to_oma_dict = dict()
    for faa_file in faa_file_list:
        _, _, faa_base, _ = sep_path_basename_ext(faa_file)
        for each_seq in SeqIO.parse(faa_file, 'fasta'):
            seq_id = each_seq.id
            gene_to_oma_dict[seq_id] = faa_base

    # get oma_to_transfer_propensity_dict
    oma_to_transfer_propensity_dict = dict()
    line_index = 0
    for each_oma in open(transfer_propensity_txt):
        if line_index > 0:
            each_oma_split = each_oma.strip().split('\t')
            oma_id = each_oma_split[0]
            transfer_propensity = float(each_oma_split[1])
            oma_to_transfer_propensity_dict[oma_id] = transfer_propensity
        line_index += 1

    fun_transfer_propensity_txt_handle = open(fun_transfer_propensity_txt, 'w')
    fun_transfer_propensity_txt_handle.write('ID\tWeighted_transfer_propensity\tDescription\n')
    oma_weighted_transfer_propensity_dict = dict()
    for fun_id in sorted(list(fun_to_gene_dict.keys())):
        current_fun_gene_set = fun_to_gene_dict[fun_id]
        current_fun_oma_stats_dict = dict()
        for gene_id in current_fun_gene_set:
            gene_oma = gene_to_oma_dict.get(gene_id, 'na')
            if gene_oma not in current_fun_oma_stats_dict:
                current_fun_oma_stats_dict[gene_oma] = 1
            else:
                current_fun_oma_stats_dict[gene_oma] += 1

        total_transfer_propensity = 0
        total_oma_num = 0
        for oma_id in current_fun_oma_stats_dict:
            oma_num = current_fun_oma_stats_dict[oma_id]
            oma_transfer_propensity = oma_to_transfer_propensity_dict.get(oma_id, 'na')
            if oma_transfer_propensity != 'na':
                total_transfer_propensity += (oma_num*oma_transfer_propensity)
                total_oma_num += oma_num

        oma_transfer_propensity_weighted = 'na'
        if total_oma_num != 0:
            oma_transfer_propensity_weighted = total_transfer_propensity/total_oma_num
            oma_transfer_propensity_weighted = float("{0:.3f}".format(oma_transfer_propensity_weighted))

        if oma_transfer_propensity_weighted != 'na':
            oma_weighted_transfer_propensity_dict[fun_id] = oma_transfer_propensity_weighted
            fun_transfer_propensity_txt_handle.write('%s\t%s\t%s\n' % (fun_id, oma_transfer_propensity_weighted, fun_id_to_desc_dict[fun_id]))

    fun_transfer_propensity_txt_handle.close()

    print('Done!')


if __name__ == '__main__':

    ALE6_parser = argparse.ArgumentParser()
    ALE6_parser.add_argument('-1',      required=True,                          help='ALE1 output directory')
    ALE6_parser.add_argument('-3',      required=True,                          help='ALE3 output directory')
    ALE6_parser.add_argument('-s',      required=True,                          help='rooted species tree')
    ALE6_parser.add_argument('-n',      required=False, default=None,           help='interested internal node(s)')
    ALE6_parser.add_argument('-cog',    required=False, default=None,           help='COG annotation results')
    ALE6_parser.add_argument('-kegg',   required=False, default=None,           help='KEGG annotation results')
    ALE6_parser.add_argument('-o',      required=True,                          help='output directory')
    ALE6_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    args = vars(ALE6_parser.parse_args())
    ALE6(args)
