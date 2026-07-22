import os
import glob
import argparse
from ete3 import Tree


ALE8_usage = '''
=============================== ALE8 example commands ===============================

TreeSAK ALE8 -3 ALE3_op_dir -o -3 ALE8_op_dir
TreeSAK ALE8 -3 ALE3_op_dir -o -3 ALE8_op_dir -l genome_tree_leaf_rename.txt

# needed input files:
-3: GeneContent, SpeciesTreeRef.newick

=====================================================================================
'''

def itol_piechart(gene_gain_loss_stats_txt, leaf_set, itol_txt):

    max_pie_size = 30
    itol_txt_handle = open(itol_txt, 'w')
    itol_txt_handle.write('DATASET_PIECHART\n')
    itol_txt_handle.write('SEPARATOR TAB\n')
    itol_txt_handle.write('DATASET_LABEL\tGene_gain_loss\n')
    itol_txt_handle.write('COLOR\t#ff0000\n')
    itol_txt_handle.write('FIELD_COLORS\t#31C950\t#FB2C36\n')
    itol_txt_handle.write('FIELD_LABELS\tGain\tLoss\n')
    itol_txt_handle.write('MARGIN\t10\n\n')
    itol_txt_handle.write('MAXIMUM_SIZE\t%s\n\n' % max_pie_size)
    itol_txt_handle.write('DATA\n')
    for each in open(gene_gain_loss_stats_txt):
        if not each.startswith('ID\t'):
            each_split = each.strip().split('\t')
            node_id = each_split[0]
            gain_num = int(each_split[1])
            loss_num = int(each_split[2])
            if node_id in leaf_set:
                itol_txt_handle.write('%s\t1\t%s\t%s\t%s\n' % (node_id, (gain_num + loss_num), gain_num, loss_num))
            else:
                itol_txt_handle.write('%s\t0.5\t%s\t%s\t%s\n' % (node_id, (gain_num + loss_num), gain_num, loss_num))
    itol_txt_handle.close()


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def ALE8(args):

    ale3_op_dir         = args['3']
    op_dir              = args['o']
    force_create_op_dir = args['f']
    rename_leaf_txt     = args['l']
    cog_annotation_wd   = args['cog']
    keg_annotation_wd   = args['kegg']

    rename_leaf_dict = dict()
    if rename_leaf_txt is not None:
        if os.path.isfile(rename_leaf_txt) is True:
            for each_line in open(rename_leaf_txt):
                each_line_split = each_line.strip().split('\t')
                rename_leaf_dict[each_line_split[0]] = each_line_split[1]

    # define file name
    gene_gain_loss_dir              = '%s/Gene_gain_loss'                   % op_dir
    gene_gain_loss_stats_txt        = '%s/Gene_gain_loss_stats.txt'         % op_dir
    gene_gain_loss_stats_txt_itol   = '%s/Gene_gain_loss_stats_iTOL.txt'    % op_dir
    species_tree_txt_renamed        = '%s/species_tree_renamed.newick'      % op_dir

    if os.path.isdir(ale3_op_dir) is False:
        print('%s not found, program exited!' % ale3_op_dir)
        exit()

    gnm_content_txt_re   = '%s/GeneContent/*.txt' % (ale3_op_dir)
    gnm_content_txt_list = glob.glob(gnm_content_txt_re)
    if len(gnm_content_txt_list) == 0:
        print('%s/GeneContent is empty, program exited!' % ale3_op_dir)
        exit()

    # read in genome content
    gnm_content_dict = dict()
    for each_content_txt in gnm_content_txt_list:
        _, _, node_id, _ = sep_path_basename_ext(each_content_txt)
        if '(' in node_id and ')' in node_id:
            node_id = node_id.split('(')[0]
        node_id = rename_leaf_dict.get(node_id, node_id)
        gnm_content_dict[node_id] = set()
        for each_line in open(each_content_txt):
            gnm_content_dict[node_id].add(each_line.strip())

    # read in tree (with internal node id)
    ale_species_tree_txt = '%s/SpeciesTreeRef.newick' % ale3_op_dir
    ale_species_tree     = Tree(ale_species_tree_txt, format=1)

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gene_gain_loss_dir)

    # get node_to_parent_dict
    leaf_set = set()
    node_to_parent_dict = dict()
    for each_node in ale_species_tree.traverse():
        node_name = each_node.name
        parent_node = each_node.up
        if parent_node is not None:
            parent_name = parent_node.name
            if each_node.is_leaf() is True:
                node_name_updated = node_name
                node_name_updated = rename_leaf_dict.get(node_name_updated, node_name_updated)
                leaf_set.add(node_name_updated)
            else:
                node_name_updated = 'node%s' % node_name
            parent_name_updated = 'node%s' % parent_name
            node_to_parent_dict[node_name_updated] = parent_name_updated

    for each_node in ale_species_tree.traverse():
        if each_node.is_leaf() is True:
            each_node.name = rename_leaf_dict.get(each_node.name, each_node.name)
        else:
            each_node.name = 'node' + each_node.name

    ale_species_tree.write(outfile=species_tree_txt_renamed, format=1)

    # get gene gain and loss
    gene_gain_num_dict = dict()
    gene_loss_num_dict = dict()
    for each_node in sorted(list(node_to_parent_dict.keys())):
        node_parent     = node_to_parent_dict[each_node]
        if each_node in leaf_set:
            each_node    = rename_leaf_dict.get(each_node, each_node)
            node_content = gnm_content_dict[each_node]
        else:
            node_content = gnm_content_dict[each_node[4:]]
        parent_content  = gnm_content_dict[node_parent[4:]]
        gene_gain_txt   = '%s/%s_gene_gain.txt' % (gene_gain_loss_dir, each_node)
        gene_loss_txt   = '%s/%s_gene_loss.txt' % (gene_gain_loss_dir, each_node)

        shared_protein_set = set(node_content).intersection(parent_content)
        child_uniq = set()
        for e1 in node_content:
            if e1 not in shared_protein_set:
                child_uniq.add(e1)
        parent_uniq = set()
        for e2 in parent_content:
            if e2 not in shared_protein_set:
                parent_uniq.add(e2)

        gene_gain_txt_handle = open(gene_gain_txt, 'w')
        gene_gain_txt_handle.write('\n'.join(sorted(list(child_uniq))))
        gene_gain_txt_handle.close()

        gene_loss_txt_handle = open(gene_loss_txt, 'w')
        gene_loss_txt_handle.write('\n'.join(sorted(list(parent_uniq))))
        gene_loss_txt_handle.close()

        gene_gain_num_dict[each_node] = len(child_uniq)
        gene_loss_num_dict[each_node] = len(parent_uniq)

    gene_gain_loss_stats_txt_handle = open(gene_gain_loss_stats_txt, 'w')
    gene_gain_loss_stats_txt_handle.write('ID\tGain\tLoss\n')
    for each_node in sorted(list(node_to_parent_dict.keys())):
        gene_gain_num = gene_gain_num_dict[each_node]
        gene_loss_num = gene_loss_num_dict[each_node]
        gene_gain_loss_stats_txt_handle.write('%s\t%s\t%s\n' % (each_node, gene_gain_num, gene_loss_num))
    gene_gain_loss_stats_txt_handle.close()

    # get iTOL file
    itol_piechart(gene_gain_loss_stats_txt, leaf_set, gene_gain_loss_stats_txt_itol)

    print('Results exported to %s' % op_dir)
    print('Done!')


if __name__ == '__main__':

    ALE8_parser = argparse.ArgumentParser()
    ALE8_parser.add_argument('-3',      required=True,                          help='ALE3 output directory')
    ALE8_parser.add_argument('-l',      required=False, default=None,           help='rename leaves')
    ALE8_parser.add_argument('-cog',    required=False, default=None,           help='COG annotation dir')
    ALE8_parser.add_argument('-kegg',   required=False, default=None,           help='KEGG annotation dir')
    ALE8_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    ALE8_parser.add_argument('-o',      required=True,                          help='output directory')
    args = vars(ALE8_parser.parse_args())
    ALE8(args)
