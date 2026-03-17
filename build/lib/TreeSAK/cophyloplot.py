import os

from ete3 import Tree

########################################################################################################################

# file in
meta_data_txt           = '/Users/songweizhi/Desktop/Sponge_r226/00_metadata/AOA_metadata_20260130.txt'
tree_file_host          = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/RefSeqs_with_AOA_18S_iden99_g_representatives_with_JL_rooted.treefile'
tree_file_symbiont      = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/AOA_tree_654_with_fixed_topo1_LG_rooted.treefile'
association_table       = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/codiv_sponge_to_AOA_all.txt'
color_code_genome_txt   = '/Users/songweizhi/Desktop/Sponge_r226/00_metadata/color_code_symbiont.txt'
color_code_sponge_txt   = ''
color_link_by_symbiont  = True

# file out
op_association_table    = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/codiv_sponge_to_AOA_op.txt'
op_pdf                  = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/test.pdf'

########################################################################################################################

current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
cophyloplot_Rscript = '%s/cophyloplot.R' % current_file_path

gnm_tax_dict = dict()
col_index = dict()
for each_gnm in open(meta_data_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    if each_gnm.startswith('Genome\t'):
        col_index = {key: i for i, key in enumerate(each_gnm_split)}
    else:
        gnm_id = each_gnm_split[col_index['Genome']]
        gnm_tax = each_gnm_split[col_index['GTDB_Taxon_r226']]
        gnm_g = 'g__'
        for each_r in gnm_tax.split(';'):
            if each_r.startswith('g__'):
                gnm_g = each_r
        gnm_tax_dict[gnm_id] = gnm_g

t_host       = Tree(tree_file_host, quoted_node_names=True, format=1)
t_symbiont   = Tree(tree_file_symbiont, quoted_node_names=True, format=1)
leaf_names_1 = t_host.get_leaf_names()
leaf_names_2 = t_symbiont.get_leaf_names()

host_color_code_dict = dict()
if os.path.isfile(color_code_sponge_txt):
    for each_line in open(color_code_sponge_txt):
        each_line_split = each_line.strip().split()
        host_color_code_dict[each_line_split[0]] = each_line_split[1]

symbiont_color_code_dict = dict()
if os.path.isfile(color_code_genome_txt):
    for each_line in open(color_code_genome_txt):
        each_line_split = each_line.strip().split()
        symbiont_color_code_dict[each_line_split[0]] = each_line_split[1]

op_association_table_handle = open(op_association_table, 'w')
for each_line in open(association_table):
    each_line_split     = each_line.strip().split()
    leaf_l              = each_line_split[0]
    leaf_r              = each_line_split[1]
    color_by_host       = host_color_code_dict.get(leaf_l, 'black')
    symbiont_tax        = gnm_tax_dict.get(leaf_r, 'g__')
    color_by_symbiont   = symbiont_color_code_dict.get(symbiont_tax, 'black')
    color_to_use = color_by_host
    if color_link_by_symbiont is True:
        color_to_use = color_by_symbiont
    if (leaf_l in leaf_names_1) and (leaf_r in leaf_names_2):
        op_association_table_handle.write('%s\t%s\t"%s"\n' % (leaf_l, leaf_r, color_to_use))
op_association_table_handle.close()

# run cophyloplot.R
r_cmd = 'Rscript %s -i %s -s %s -l %s -o %s' % (cophyloplot_Rscript, tree_file_host, tree_file_symbiont, op_association_table, op_pdf)
print(r_cmd)
os.system(r_cmd)




