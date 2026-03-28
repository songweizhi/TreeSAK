import os
from ete3 import Tree

def iTOL_Tanglegram():

    ########################################################################################################################

    # file in
    meta_data_txt           = '/Users/songweizhi/Desktop/Sponge_r226/00_metadata/AOA_metadata_20260130.txt'
    tree_file_host          = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/RefSeqs_with_AOA_18S_iden99_g_representatives_JL_wd/Sponge_phylogeny_Maria_topo_genus_level.tree'
    tree_file_symbiont      = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/18_sym.treefile'
    association_table       = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/codiv_sponge_to_AOA_all.txt'
    color_code_genome_txt   = '/Users/songweizhi/Desktop/Sponge_r226/00_metadata/color_code_symbiont.txt'
    color_code_sponge_txt   = ''
    color_link_by_host      = True
    color_link_by_symbiont  = False

    # tanglegram tree position
    tree_pan_h      = 3000
    tree_pan_v      = 1500

    # tanglegram tree scaling
    tree_scale_h    = 1
    tree_scale_v    = 7

    # tanglegram tree label
    show_leaf_label = 1
    leaf_label_size_factor = 6
    leaf_label_shift = 0
    leaf_label_align = 'left'

    # linkages
    max_link_line_width = 3
    connection_margin_left = 12
    connection_margin_right = 3

    # file out
    itol_tanglegram_txt     = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/iTOL_Tanglegram.txt'

    ########################################################################################################################

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

    tanglegram_tree_str = ''
    with open(tree_file_host, 'r') as f:
        tanglegram_tree_str = f.readline()

    itol_tanglegram_txt_handle = open(itol_tanglegram_txt, 'w')
    itol_tanglegram_txt_handle.write('DATASET_TANGLEGRAM\nSEPARATOR TAB\n\n')
    itol_tanglegram_txt_handle.write('DATASET_LABEL\tTanglegram\nCOLOR\t#00ff00\n\n')
    itol_tanglegram_txt_handle.write('TANGLEGRAM_TREE\n%s\nEND_TANGLEGRAM_TREE\n\n' % tanglegram_tree_str)
    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('# Connection lines\n')
    itol_tanglegram_txt_handle.write('CONNECTION_CURVE\t0\n')
    itol_tanglegram_txt_handle.write('MAXIMUM_LINE_WIDTH\t%s\n' % max_link_line_width)
    itol_tanglegram_txt_handle.write('CONNECTION_MARGIN_LEFT\t%s\nCONNECTION_MARGIN_RIGHT\t%s\n' % (connection_margin_left, connection_margin_right))

    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('\n')

    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('# Move the tanglegram tree horizontally and/or vertically\n')
    itol_tanglegram_txt_handle.write('TREE_PAN_H\t%s\nTREE_PAN_V\t%s\n'     % (tree_pan_h, tree_pan_v))
    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('# Tanglegram tree scaling factors, horizontal and vertical\n')
    itol_tanglegram_txt_handle.write('TREE_SCALE_H\t%s\nTREE_SCALE_V\t%s\n' % (tree_scale_h, tree_scale_v))
    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('# Tanglegram tree label\n')
    itol_tanglegram_txt_handle.write('SHOW_LABELS\t%s\nSIZE_FACTOR\t%s\nLABEL_SHIFT\t%s\nLABEL_ALIGN\t%s\n' % (show_leaf_label, leaf_label_size_factor, leaf_label_shift, leaf_label_align))



    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('\n')
    itol_tanglegram_txt_handle.write('DATA\n')
    itol_tanglegram_txt_handle.write('#connect\tMAIN_TREE_NODE_ID\tTANGLEGRAM_TREE_NODE_ID\tWIDTH\tCOLOR\tSTYLE\tLABEL\n')
    itol_tanglegram_txt_handle.write('#colorstrip\tCOLUMN_INDEX\tTANGLEGRAM_TREE_NODE_ID\tCOLOR LABEL\n')
    itol_tanglegram_txt_handle.write('#style\tTANGLEGRAM_TREE_NODE_ID\tTYPE\tWHAT\tCOLOR\tWIDTH_OR_SIZE_FACTOR\tSTYLE\tBACKGROUND_COLOR\n')
    itol_tanglegram_txt_handle.write('#symbol\tTANGLEGRAM_TREE_NODE_ID\tSYMBOL\tSIZE\tCOLOR\tFILL\tPOSITION\tLABEL\n')
    itol_tanglegram_txt_handle.write('\n')
    for each_line in open(association_table):
        each_line_split     = each_line.strip().split()
        leaf_l              = each_line_split[0]
        leaf_r              = each_line_split[1]
        color_by_host       = host_color_code_dict.get(leaf_l, 'black')
        symbiont_tax        = gnm_tax_dict.get(leaf_r, 'g__')
        color_by_symbiont   = symbiont_color_code_dict.get(symbiont_tax, 'black')

        # link color
        link_color_to_use = 'black'
        if color_link_by_symbiont is True:
            link_color_to_use = color_by_symbiont
        elif color_link_by_host is True:
            link_color_to_use = color_by_host

        itol_tanglegram_txt_handle.write('connect\t%s\t%s\t3\t%s\tnormal\ttest\n' % (leaf_r, leaf_l, link_color_to_use))
    itol_tanglegram_txt_handle.close()

iTOL_Tanglegram()
