import os
import glob
import math
import random
import argparse
import seaborn as sns
from Bio import SeqIO
from ete3 import Tree
from PIL import Image
from itolapi import Itol
import multiprocessing as mp
from PyPDF3.pdf import PageObject
from PyPDF3 import PdfFileWriter, PdfFileReader
from ete3 import TextFace, TreeStyle, NodeStyle


ALE_usage = '''
========================= ALE example commands =========================

TreeSAK ALE1  ->  step 1: get gene tree
TreeSAK ALE2  ->  step 2: run ALE
TreeSAK ALE3  ->  step 3: parse ALE output
TreeSAK ALE4  ->  Infer ancestral genome

cd /Users/songweizhi/Desktop/demo
TreeSAK ALE1 -i OrthologousGroups.txt -s combined_d__Archaea_o_rs.faa -p oma -c genome_taxon.txt -m 50 -n 2 -t 6 -jt 3 -f -o ALE1_op_dir
TreeSAK ALE2 -i ALE1_op_dir -s Marker_set_1_PA_75_C60_PMSF_concatenated_rooted.treefile -c genome_taxon.txt -t 6 -f -o ALE2_op_dir
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ


ALEobserve OMA00003_for_ALE.ufboot


cd /Users/songweizhi/Desktop/demo/subset_gnm_tree_no
docker run -v $PWD:$PWD -w $PWD gregmich/alesuite_new ALEobserve OMA00003_for_ALE.ufboot
docker run -v $PWD:$PWD -w $PWD gregmich/alesuite_new ALEml_undated genome_tree_rooted.treefile OMA00003_for_ALE.ufboot.ale

cd /Users/songweizhi/Desktop/demo/subset_gnm_tree_yes
docker run -v $PWD:$PWD -w $PWD gregmich/alesuite_new ALEobserve OMA00003_for_ALE.ufboot
docker run -v $PWD:$PWD -w $PWD gregmich/alesuite_new ALEml_undated OMA00003_genome_tree_for_ALE.treefile OMA00003_for_ALE.ufboot.ale



========================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def select_seq(arg_list):

    seq_file    = arg_list[0]
    id_file     = arg_list[1]
    output_file = arg_list[2]

    seq_id_set = {i.strip() for i in open(id_file)}
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        if seq_record.id in seq_id_set:
            SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
    output_file_handle.close()


def root_with_out_group(tree_file, out_group_txt, tree_file_rooted):

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=1)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted)


def subset_tree(tree_file_in, leaves_to_keep_list, tree_file_out):

    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(leaves_to_keep_list, preserve_branch_length=True)
    if tree_file_out is None:
        return subset_tree.write()
    else:
        subset_tree.write(outfile=tree_file_out)


def plot_tree(input_tree, tree_title, node_label_dict, node_label_color_dict, align_leaf_label, show_scale, output_plot):

    if os.path.isfile(input_tree) is False:
        print('Tree file not found, program exited!')
        print(input_tree)
        exit()

    t = Tree(input_tree)
    ts = TreeStyle()
    ts.mode = "r"                       # tree model: 'r' for rectangular, 'c' for circular
    ts.show_border = False              # set tree image border
    ts.show_leaf_name = False           # show/hide leaf name, hide here, so you can customise it below with node.add_face()
    ts.title.add_face(TextFace(tree_title, fsize=9, fgcolor='black', ftype='Arial', tight_text=False), column=0)  # add tree title

    # set node style
    for each_node in t.traverse():
        ns = NodeStyle()
        ns["shape"]         = "circle"  # dot shape: circle, square or sphere
        ns["fgcolor"]       = "black"   # color of shape(not label)
        ns['size']          = 0         # node shape size
        ns['hz_line_type']  = 0         # horizontal branch line type: 0 for solid, 1 for dashed, 2 for dotted
        ns['vt_line_type']  = 0         # vertical branch line type:   0 for solid, 1 for dashed, 2 for dotted
        ns['hz_line_width'] = 0.5       # horizontal branch line width
        ns['vt_line_width'] = 0.5       # vertical branch line width

        leaf_label_position = 'branch-right'
        if align_leaf_label is True:
            leaf_label_position = 'aligned'

        if each_node.is_leaf():
            node_id = each_node.name
            node_label_color = node_label_color_dict.get(node_id, 'black')
            node_label_text  = node_label_dict.get(node_id, node_id)
            each_node.add_face(TextFace(node_label_text, fsize=8, fgcolor=node_label_color, tight_text=False, bold=False),
                               column=0, position=leaf_label_position)  # aligned, branch-right
        else:
            pass
        each_node.set_style(ns)

    # set layout
    ts.rotation                 = 0             # from 0 to 360
    ts.margin_top               = 10            # top tree image margin
    ts.margin_bottom            = 10            # bottom tree image margin
    ts.margin_left              = 10            # left tree image margin
    ts.margin_right             = 10            # right tree image margin
    ts.branch_vertical_margin   = 3             # 3 pixels between adjancent branches
    ts.show_scale               = show_scale    # show_scale
    ts.show_border              = False         # set tree image border

    # write out tree
    t.render(output_plot, w=1200, units="px", tree_style=ts)


def merge_image(image_file_list, output_image):

    images = [Image.open(x) for x in image_file_list]
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    new_im = Image.new('RGB', (total_width, max_height), color='white')

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    new_im.save(output_image)


def merge_pdf(pdf_1, pdf_2, margin_size, op_pdf):

    page1 = PdfFileReader(open(pdf_1, "rb"), strict=False).getPage(0)
    page2 = PdfFileReader(open(pdf_2, "rb"), strict=False).getPage(0)

    total_width  = page1.mediaBox.upperRight[0] + page2.mediaBox.upperRight[0] + margin_size*3
    total_height = max([page1.mediaBox.upperRight[1], page2.mediaBox.upperRight[1]]) + margin_size*2

    new_page = PageObject.createBlankPage(None, total_width, total_height)
    new_page.mergeTranslatedPage(page1, margin_size, (total_height-margin_size-page1.mediaBox.upperRight[1]))
    new_page.mergeTranslatedPage(page2, (page1.mediaBox.upperRight[0] + margin_size*2), margin_size)

    output = PdfFileWriter()
    output.addPage(new_page)
    output.write(open(op_pdf, "wb"))


def uts_to_itol_connections(genome_tree_file, ale_formatted_gnm_tree, interal_node_prefix, uts_file, freq_cutoff, ignore_leaf_hgt, ignore_vertical_hgt, donor_node_min_leaf_num, recipient_node_min_leaf_num, itol_connection_txt, dr_separator):

    # get internal_node_to_leaf_dict
    internal_node_to_leaf_dict = get_node_to_leaf_dict(ale_formatted_gnm_tree)

    paired_donor_to_recipient_leaf_dict = dict()
    qualified_hgt_num = 0

    leaf_id_set = []
    if os.path.isfile(genome_tree_file):
        leaf_id_set = [i.name for i in Tree(genome_tree_file, format=3).get_leaves()]
    else:
        print('%s not found!' % genome_tree_file)

    hgt_freq_dict = dict()
    connection_line_to_write_dict = dict()
    with open(itol_connection_txt, 'w') as itol_connection_txt_handle:
        itol_connection_txt_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
        itol_connection_txt_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t1\nARROW_SIZE\t60\nLOOP_SIZE\t100\n')
        itol_connection_txt_handle.write('MAXIMUM_LINE_WIDTH\t10\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t0\nDATA\n')
        for each_line in open(uts_file):
            if not each_line.startswith('#'):
                each_line_split = each_line.strip().split('\t')
                donor = each_line_split[0]
                recipient = each_line_split[1]
                freq = float(each_line_split[2])

                # add prefix to internal donor node
                if donor in leaf_id_set:
                    donor_with_prefix = donor
                else:
                    donor_with_prefix = interal_node_prefix + donor

                # add prefix to internal recipient node
                if recipient in leaf_id_set:
                    recipient_with_prefix = recipient
                else:
                    recipient_with_prefix = interal_node_prefix + recipient

                key_str      = '%s%s%s' % (donor_with_prefix, dr_separator, recipient_with_prefix)

                line_to_write = ''
                if freq >= freq_cutoff:
                    if ignore_leaf_hgt is False:
                        if ignore_vertical_hgt is False:
                            line_to_write = '%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq)
                            qualified_hgt_num += 1
                        else:
                            donor_is_ancestor_of_recipient = check_a_is_ancestor_of_b(ale_formatted_gnm_tree, donor, recipient)
                            donor_is_child_of_recipient    = check_a_is_child_of_b(ale_formatted_gnm_tree, donor, recipient)
                            if (donor_is_ancestor_of_recipient is False) and (donor_is_child_of_recipient is False):
                                line_to_write = '%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq)
                                qualified_hgt_num += 1
                    else:
                        if (each_line_split[0] not in leaf_id_set) and (each_line_split[1] not in leaf_id_set):
                            donor_node_leaf_num = len(internal_node_to_leaf_dict.get(donor, []))
                            recipient_node_leaf_num = len(internal_node_to_leaf_dict.get(recipient, []))
                            if (donor_node_leaf_num >= donor_node_min_leaf_num) and (recipient_node_leaf_num >= recipient_node_min_leaf_num):
                                if ignore_vertical_hgt is False:
                                    line_to_write = '%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq)
                                    qualified_hgt_num += 1
                                else:
                                    donor_is_ancestor_of_recipient = check_a_is_ancestor_of_b(ale_formatted_gnm_tree, donor, recipient)
                                    donor_is_child_of_recipient    = check_a_is_child_of_b(ale_formatted_gnm_tree, donor, recipient)
                                    if (donor_is_ancestor_of_recipient is False) and (donor_is_child_of_recipient is False):
                                        line_to_write = '%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq)
                                        qualified_hgt_num += 1
                                        paired_donor_to_recipient_leaf_dict[key_str] = [internal_node_to_leaf_dict.get(donor, []), internal_node_to_leaf_dict.get(recipient, [])]

                if line_to_write != '':
                    itol_connection_txt_handle.write(line_to_write)
                    connection_line_to_write_dict[key_str] = line_to_write
                    hgt_freq_dict[key_str] = freq

    combined_connection_file_path, combined_connection_file_basename, combined_connection_file_ext = sep_path_basename_ext(itol_connection_txt)

    # write out connections separately
    for each_connection in connection_line_to_write_dict:
        pwd_connection_txt = '%s/%s_%s.txt' % (combined_connection_file_path, combined_connection_file_basename, each_connection)
        pwd_connection_txt_handle = open(pwd_connection_txt, 'w')
        pwd_connection_txt_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
        pwd_connection_txt_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t1\nARROW_SIZE\t60\nLOOP_SIZE\t100\n')
        pwd_connection_txt_handle.write('MAXIMUM_LINE_WIDTH\t10\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t0\nDATA\n')
        pwd_connection_txt_handle.write(connection_line_to_write_dict[each_connection] + '\n')
        pwd_connection_txt_handle.close()

    return internal_node_to_leaf_dict, paired_donor_to_recipient_leaf_dict, hgt_freq_dict


def itol_tree(tree_file, annotation_file_list, project_name, APIkey, display_mode, op_plot):

    # https://github.com/albertyw/itolapi
    # http://itol.embl.de/help.cgi#batch

    op_plot_ext = op_plot.split('.')[-1]

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = APIkey  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = tree_file
    itol_uploader.add_file(tree_file)

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        itol_uploader.add_file(annotation_file)

    status = itol_uploader.upload()
    # import pdb;pdb.set_trace()
    assert status != False

    # the following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    if len(annotation_file_list) == 1:
        datasets_visible_str = '0'
    elif len(annotation_file_list) == 2:
        datasets_visible_str = '0,1'
    elif len(annotation_file_list) == 3:
        datasets_visible_str = '0,1,2'
    else:
        datasets_visible_str = ','.join([str(i) for i in list(range(0, len(annotation_file_list)))])
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', datasets_visible_str)
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '1')
    # itol_exporter.set_export_param_value('current_font_size', '96')
    itol_exporter.set_export_param_value('line_width', '3')
    itol_exporter.set_export_param_value('vertical_shift_factor', '0.9')
    itol_exporter.set_export_param_value('horizontal_scale_factor', '0.9')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


def get_node_to_leaf_dict(tree_file):
    internal_node_to_leaf_dict = dict()
    for node in Tree(tree_file, format=1).traverse():
        if not node.is_leaf():
            node_name = node.name
            node_leaf_list = node.get_leaf_names()
            internal_node_to_leaf_dict[node_name] = node_leaf_list
    return internal_node_to_leaf_dict


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


def root_at_midpoint(tree_in, tree_in_rooted):
    t = Tree(tree_in)
    midpoint = t.get_midpoint_outgroup()
    t.set_outgroup(midpoint)
    t.write(outfile=tree_in_rooted)


def get_color_list(color_num):
    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a',
                               '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']
    else:
        color_num_each = math.ceil(color_num / 8) + 2
        color_list_1 = sns.color_palette('Blues', n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy', n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green', n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds', n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive', n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6,
                           color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


def scale_str_to_size_list(scale_str):
    scale_list = scale_str.split('-')
    scale_list = [float(i) for i in scale_list]

    shape_size_list = []
    if scale_list[0] == 0:
        shape_size_list = [0]
        for each_value in scale_list[1:-1]:
            current_size = each_value / scale_list[-1]
            shape_size_list.append(current_size)
        shape_size_list.append(1)

    if scale_list[0] != 0:
        shape_size_list = [0.1]
        interval_num = len(scale_list) - 1
        interval_value = (1 - 0.1) / interval_num
        n = 1
        for each_value in scale_list[1:-1]:
            shape_size_list.append(interval_value * n + 0.1)
            n += 1
        shape_size_list.append(1)

    return shape_size_list


def get_ortho_to_gene_dict(ortho_groups_txt, og_program):

    ortho_to_gene_dict = dict()
    for each_og in open(ortho_groups_txt):
        if not each_og.startswith('#'):
            og_id = ''
            gene_list = []
            if og_program == 'orthofinder':
                each_og_split = each_og.strip().split(' ')
                og_id = each_og_split[0][:-1]
                gene_list = each_og_split[1:]
            elif og_program == 'oma':
                each_og_split = each_og.strip().split('\t')
                og_id = each_og_split[0]
                group_member_list = each_og_split[1:]
                for each_protein in group_member_list:
                    protein_id = each_protein.split(' ')[0].split(':')[1]
                    gene_list.append(protein_id)
            ortho_to_gene_dict[og_id] = gene_list

    return ortho_to_gene_dict


def prepare_ale_ip_worker(arg_list):

    qualified_og                    = arg_list[0]
    gene_tree_dir                   = arg_list[1]
    ale_wd                          = arg_list[2]
    genome_tree_file_rooted         = arg_list[3]
    gnm_pco_dict                    = arg_list[4]
    gene_tree_ufboot_for_ale        = arg_list[5]
    genome_tree_file_subset_for_ale = arg_list[6]

    genome_tree_file_subset             = '%s_genome_tree.treefile'             % qualified_og
    gene_tree_ufboot                    = '%s.ufboot'                           % qualified_og
    gene_tree_treefile                  = '%s.treefile'                         % qualified_og
    gene_tree_treefile_subset           = '%s_subset.treefile'                  % qualified_og
    pwd_genome_tree_file_subset         = '%s/%s'                               % (gene_tree_dir, genome_tree_file_subset)
    pwd_genome_tree_file_subset_for_ale = '%s/%s'                               % (ale_wd, genome_tree_file_subset_for_ale)
    pwd_gene_tree_ufboot                = '%s/%s'                               % (gene_tree_dir, gene_tree_ufboot)
    pwd_gene_tree_ufboot_for_ale        = '%s/%s'                               % (ale_wd, gene_tree_ufboot_for_ale)
    pwd_gene_tree_treefile              = '%s/%s'                               % (gene_tree_dir, gene_tree_treefile)
    pwd_gene_tree_treefile_subset       = '%s/%s'                               % (gene_tree_dir, gene_tree_treefile_subset)

    # get genomes on gene tree
    gene_gnm_set = set()
    gnm_to_gene_dict = dict()
    for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():

        # get gnm id
        gene_gnm = '_'.join(each_gene.split('_')[:-1])

        gene_gnm_set.add(gene_gnm)
        if gene_gnm not in gnm_to_gene_dict:
            gnm_to_gene_dict[gene_gnm] = {each_gene}
        else:
            gnm_to_gene_dict[gene_gnm].add(each_gene)

    # subset genome tree
    genome_tree_leaf_set = Tree(genome_tree_file_rooted).get_leaf_names()
    gnms_in_both_trees = set(genome_tree_leaf_set).intersection(gene_gnm_set)
    gnm_tree_subset_str = subset_tree(genome_tree_file_rooted, gnms_in_both_trees, None)
    gnm_tree_subset_str_for_ale = gnm_tree_subset_str
    gnm_tree_subset_str_for_ale = gnm_tree_subset_str_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')

    # write out genome tree subset
    with open(pwd_genome_tree_file_subset, 'w') as pwd_genome_tree_file_subset_handle:
        pwd_genome_tree_file_subset_handle.write(gnm_tree_subset_str)

    # write out genome tree subset for running ALE
    with open(pwd_genome_tree_file_subset_for_ale, 'w') as pwd_genome_tree_file_subset_for_ale_handle:
        pwd_genome_tree_file_subset_for_ale_handle.write(gnm_tree_subset_str_for_ale)

    # get genes to keep in gene tree
    gene_set_to_keep = set()
    for each_gnm in gnms_in_both_trees:
        gene_set_to_keep.update(gnm_to_gene_dict.get(each_gnm, set()))

    # subset gene_tree.treefile
    subset_tree(pwd_gene_tree_treefile, gene_set_to_keep, pwd_gene_tree_treefile_subset)

    # subset gene_tree.ufboot and rename leaves for running ALE
    pwd_gene_tree_ufboot_for_ale_handle = open(pwd_gene_tree_ufboot_for_ale, 'w')
    for each_gene_tree in open(pwd_gene_tree_ufboot):
        gene_tree_str = each_gene_tree.strip()
        gene_tree_str_subset_for_ale = subset_tree(gene_tree_str, gene_set_to_keep, None)
        gene_tree_str_subset_for_ale = gene_tree_str_subset_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
        pwd_gene_tree_ufboot_for_ale_handle.write(gene_tree_str_subset_for_ale + '\n')
    pwd_gene_tree_ufboot_for_ale_handle.close()

    # get gene tree leaf name dict (for plot)
    leaf_name_dict = dict()
    for each_gene in Tree(pwd_gene_tree_treefile_subset).get_leaf_names():
        gene_id = each_gene
        gene_genome = '_'.join(gene_id.split('_')[:-1])
        genome_pco = gnm_pco_dict[gene_genome]
        gene_id_with_taxon = '%s_%s' % (genome_pco, gene_id.split('_')[-1])
        leaf_name_dict[gene_id] = gene_id_with_taxon


def ALE1(args):

    orthogroups_op_txt  = args['i']
    combined_faa        = args['s']
    og_program          = args['p']
    genome_taxon_txt    = args['c']
    min_og_genome_num   = args['m']
    min_og_phylum_num   = args['n']
    num_threads         = args['t']
    js_num_threads      = args['jt']
    force_create_op_dir = args['f']
    op_dir              = args['o']
    designate_ogs       = []
    to_ignore_ogs_list  = []

    # define output file name
    get_gene_tree_cmds_txt = '%s_cmds.txt' % op_dir

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    if force_create_op_dir is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
    os.system('mkdir %s' % op_dir)

    # read in genome taxonomy
    gnm_p_dict   = dict()
    gnm_c_dict   = dict()
    gnm_o_dict   = dict()
    gnm_pco_dict = dict()
    for each_gnm in open(genome_taxon_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id     = each_gnm_split[0]
        taxon_str  = each_gnm_split[1]
        gnm_phylum = taxon_str.split(';')[1]
        gnm_class  = taxon_str.split(';')[2]
        gnm_order  = taxon_str.split(';')[3]
        gnm_p_dict[gnm_id] = gnm_phylum
        gnm_c_dict[gnm_id] = gnm_class
        gnm_o_dict[gnm_id] = gnm_order
        gnm_pco_dict[gnm_id] = '%s__%s__%s__%s' % (gnm_phylum[3:], gnm_class[3:], gnm_order[3:], gnm_id)

    # get ortho_to_gene_dict
    ortho_to_gene_dict = get_ortho_to_gene_dict(orthogroups_op_txt, og_program)

    # get qualified orthogroups
    qualified_og_set = set()
    for each_ortho in ortho_to_gene_dict:
        ortho_gene_set = ortho_to_gene_dict[each_ortho]
        ortho_p_set = set()
        ortho_gnm_set = set()
        for each_gene in ortho_gene_set:
            gene_gnm = '_'.join(each_gene.split('_')[:-1])
            gnm_taxon = gnm_p_dict[gene_gnm]
            ortho_gnm_set.add(gene_gnm)
            ortho_p_set.add(gnm_taxon)
        if (len(ortho_gnm_set) >= min_og_genome_num) and (len(ortho_p_set) >= min_og_phylum_num):
            qualified_og_set.add(each_ortho)
    print('The total number of identified orthogroups is %s.' % len(ortho_to_gene_dict))
    print('The number of orthogroups spanning >= %s genomes and >= %s phyla is %s.' % (min_og_genome_num, min_og_phylum_num, len(qualified_og_set)))

    # process qualified OG
    og_to_process = sorted([i for i in qualified_og_set])
    if len(designate_ogs) > 0:
        print('The number of designated OGs to process: %s' % len(designate_ogs))
        og_to_process = designate_ogs

    og_to_process_no_ignored = set()
    for each_og in og_to_process:
        if each_og not in to_ignore_ogs_list:
            og_to_process_no_ignored.add(each_og)

    # extract gene sequences and prepare commands for building gene tree
    print('Preparing commands for building gene trees')
    extract_seq_arg_lol = []
    prepare_ale_ip_worker_arg_lol = []
    get_gene_tree_cmds_txt_handle = open(get_gene_tree_cmds_txt, 'w')
    for qualified_og in sorted(og_to_process_no_ignored):
        qualified_og_gene_set         = ortho_to_gene_dict[qualified_og]
        qualified_og_gene_txt         = '%s/%s.txt'         % (op_dir, qualified_og)
        qualified_og_gene_faa         = '%s/%s.faa'         % (op_dir, qualified_og)
        qualified_og_gene_aln         = '%s/%s.aln'         % (op_dir, qualified_og)
        qualified_og_gene_aln_trimmed = '%s/%s_trimmed.aln' % (op_dir, qualified_og)
        pwd_gene_tree_ufboot          = '%s/%s.ufboot'      % (op_dir, qualified_og)

        # write out the id of genes
        with open(qualified_og_gene_txt, 'w') as qualified_og_gene_txt_handle:
            qualified_og_gene_txt_handle.write('\n'.join(qualified_og_gene_set))

        # add to mp lol
        extract_seq_arg_lol.append([combined_faa, qualified_og_gene_txt, qualified_og_gene_faa])

        # write out js for mafft, trimal and iqtree
        mafft_cmd           = 'mafft-einsi --thread %s --quiet %s.faa > %s.aln'                     % (js_num_threads, qualified_og, qualified_og)
        trimal_cmd          = 'trimal -in %s.aln -out %s -automated1'                               % (qualified_og, qualified_og_gene_aln_trimmed)
        iqtree_cmd          = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s.aln -pre %s'           % (js_num_threads, qualified_og, qualified_og)
        iqtree_cmd_trimmed  = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s_trimmed'       % (js_num_threads, qualified_og_gene_aln_trimmed, qualified_og)
        get_gene_tree_cmds_txt_handle.write('%s; %s\n' % (mafft_cmd, iqtree_cmd))
    get_gene_tree_cmds_txt_handle.close()

    # extract gene sequences with multiprocessing
    print('Extracting gene sequences with %s cores' % num_threads)
    pool = mp.Pool(processes=num_threads)
    pool.map(select_seq, extract_seq_arg_lol)
    pool.close()
    pool.join()


def ALE2(args):

    gene_tree_dir           = args['i']
    genome_tree_file_rooted = args['s']
    genome_taxon_txt        = args['c']
    force_create_ale_wd     = args['f']
    num_threads             = args['t']
    ale_wd                  = args['o']
    designate_ogs           = []
    to_ignore_ogs_list      = []
    run_ale_cmds_txt        = '%s_cmds.txt' % ale_wd

    ufboot_file_re   = '%s/*.ufboot' % gene_tree_dir
    ufboot_file_list = glob.glob(ufboot_file_re)

    og_to_process_list = []
    for each_ufboot in ufboot_file_list:
        _, ufboot_base, _ = sep_path_basename_ext(each_ufboot)
        og_to_process_list.append(ufboot_base)

    # read in genome taxonomy
    gnm_pco_dict = dict()
    for each_gnm in open(genome_taxon_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id     = each_gnm_split[0]
        taxon_str  = each_gnm_split[1]
        gnm_phylum = taxon_str.split(';')[1]
        gnm_class  = taxon_str.split(';')[2]
        gnm_order  = taxon_str.split(';')[3]
        gnm_pco_dict[gnm_id] = '%s__%s__%s__%s' % (gnm_phylum[3:], gnm_class[3:], gnm_order[3:], gnm_id)

    # create ale_wd
    if force_create_ale_wd is True:
        if os.path.isdir(ale_wd) is True:
            os.system('rm -r %s' % ale_wd)
    os.system('mkdir %s' % ale_wd)

    run_ale_cmds_txt_handle = open(run_ale_cmds_txt, 'w')
    prepare_ale_ip_worker_arg_lol = []
    for qualified_og in og_to_process_list:
        pwd_gene_tree_ufboot = '%s/%s.ufboot' % (gene_tree_dir, qualified_og)
        if os.path.isfile(pwd_gene_tree_ufboot) is False:
            print('%s not found, please build gene tree first!' % pwd_gene_tree_ufboot)
        else:
            gene_tree_ufboot_for_ale        = '%s_for_ALE.ufboot'                      % qualified_og
            genome_tree_file_subset_for_ale = '%s_genome_tree_for_ALE.treefile'        % qualified_og
            obtain_ale_file_cmd             = 'ALEobserve %s'                          % gene_tree_ufboot_for_ale
            reconciliation_cmd              = 'ALEml_undated %s %s_for_ALE.ufboot.ale' % (genome_tree_file_subset_for_ale, qualified_og)
            current_arg_list                = [qualified_og, gene_tree_dir, ale_wd, genome_tree_file_rooted, gnm_pco_dict, gene_tree_ufboot_for_ale, genome_tree_file_subset_for_ale]
            run_ale_cmds_txt_handle.write('%s; %s\n' % (obtain_ale_file_cmd, reconciliation_cmd))
            prepare_ale_ip_worker_arg_lol.append(current_arg_list)
    run_ale_cmds_txt_handle.close()

    # prepare input files and job script for running ALE with multiprocessing
    print('Preparing files for running ALE with %s cores for %s OGs' % (num_threads, len(prepare_ale_ip_worker_arg_lol)))
    pool = mp.Pool(processes=num_threads)
    pool.map(prepare_ale_ip_worker, prepare_ale_ip_worker_arg_lol)
    pool.close()
    pool.join()


def iTOL(Leaf_to_Group_dict, Group_to_Color_dict, FileOut):

    Group_set = set()
    for each_leaf in Leaf_to_Group_dict:
        Group_set.add(Leaf_to_Group_dict[each_leaf])

    if len(Group_to_Color_dict) == 0:
        Group_to_Color_dict = dict(zip(Group_set, get_color_list(len(Group_set))))
    else:
        group_without_color_list = []
        for each_group in Group_set:
            if each_group not in Group_to_Color_dict:
                group_without_color_list.append(each_group)
        if len(group_without_color_list) > 0:
            color_list_unprovided = get_color_list(len(group_without_color_list))
            Group_to_Color_dict_unprovided = dict(zip(group_without_color_list, color_list_unprovided))
            for each_group in Group_to_Color_dict_unprovided:
                Group_to_Color_dict[each_group] = Group_to_Color_dict_unprovided[each_group]

    FileOut_handle = open(FileOut, 'w')
    FileOut_handle.write('DATASET_COLORSTRIP\n')
    FileOut_handle.write('SEPARATOR TAB\n')
    FileOut_handle.write('DATASET_LABEL\tTaxonomy\n')
    FileOut_handle.write('\n# customize strip attributes here\n')
    FileOut_handle.write('STRIP_WIDTH\t100\n')
    FileOut_handle.write('MARGIN\t20\n')
    FileOut_handle.write('\n# provide data here\nDATA\n')
    for leaf in Leaf_to_Group_dict:
        leaf_group = Leaf_to_Group_dict[leaf]
        leaf_color = Group_to_Color_dict[leaf_group]
        FileOut_handle.write('%s\t%s\t%s\n' % (leaf, leaf_color, leaf_group))
    FileOut_handle.close()


def parse_ale_op_worker(arg_list):

    qualified_og                = arg_list[0]
    gene_tree_dir               = arg_list[1]
    ale_wd                      = arg_list[2]
    ale_op_dir                  = arg_list[3]
    ale_hgt_plot_dir            = arg_list[4]
    interal_node_prefix         = arg_list[5]
    gnm_pco_dict                = arg_list[6]
    d_color                     = arg_list[7]
    r_color                     = arg_list[8]
    project_name                = arg_list[9]
    API_key                     = arg_list[10]
    display_mode                = arg_list[11]
    hgt_freq_cutoff             = arg_list[12]
    ignore_leaf_hgt             = arg_list[13]
    ignore_vertical_hgt         = arg_list[14]
    donor_node_min_leaf_num     = arg_list[15]
    recipient_node_min_leaf_num = arg_list[16]
    dr_separator                = arg_list[17]
    root_gene_tree_at_midpoint  = arg_list[18]
    p_color_txt                 = arg_list[19]

    ale_uml_rec_file                                    = '%s/%s_for_ALE.ufboot.ale.uml_rec'                    % (ale_wd, qualified_og)
    gene_tree_treefile                                  = '%s.treefile'                                         % qualified_og
    genome_tree_file_subset_for_ale                     = '%s_genome_tree_for_ALE.treefile'                     % qualified_og
    gene_tree_ufboot_for_ale                            = '%s_for_ALE.ufboot'                                   % qualified_og
    uts_file                                            = '%s.ale.uTs'                                          % gene_tree_ufboot_for_ale
    uml_rec_file                                        = '%s.ale.uml_rec'                                      % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree                              = '%s_ALE_formatted_genome_tree.tree'                   % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree_with_len                     = '%s_ALE_formatted_genome_tree_with_len.tree'          % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree_with_len_prefixed            = '%s_ALE_formatted_genome_tree_with_len_prefixed.tree' % gene_tree_ufboot_for_ale
    itol_connection_txt_all                             = '%s_iTOL_connection.txt'                              % qualified_og
    itol_label_txt                                      = '%s_iTOL_genome_pco.txt'                              % qualified_og
    gene_tree_itol_label_txt                            = '%s_iTOL_gene_pco.txt'                                % qualified_og
    gene_tree_treefile_subset                           = '%s_subset.treefile'                                  % qualified_og
    gene_tree_treefile_subset_midpoint_rooted           = '%s_subset_midpoint_rooted.treefile'                  % qualified_og
    gene_tree_itol_colorstrip_txt                       = '%s_iTOL_colorstrip_gene.txt'                         % qualified_og
    genome_tree_itol_colorstrip_txt                     = '%s_iTOL_colorstrip_genome.txt'                       % qualified_og
    pwd_gene_tree_treefile_subset                       = '%s/%s'                                               % (gene_tree_dir, gene_tree_treefile_subset)
    pwd_gene_tree_treefile_subset_midpoint_rooted       = '%s/%s'                                               % (ale_op_dir, gene_tree_treefile_subset_midpoint_rooted)
    pwd_gene_tree_treefile                              = '%s/%s'                                               % (gene_tree_dir, gene_tree_treefile)
    pwd_genome_tree_file_subset_for_ale                 = '%s/%s'                                               % (ale_op_dir, genome_tree_file_subset_for_ale)
    pwd_itol_connection_txt_all                         = '%s/%s'                                               % (ale_hgt_plot_dir, itol_connection_txt_all)
    pwd_itol_label_txt                                  = '%s/%s'                                               % (ale_op_dir, itol_label_txt)
    pwd_gene_tree_itol_label_txt                        = '%s/%s'                                               % (ale_hgt_plot_dir, gene_tree_itol_label_txt)
    pwd_uts_file                                        = '%s/%s'                                               % (ale_op_dir, uts_file)
    pwd_uml_rec_file                                    = '%s/%s'                                               % (ale_op_dir, uml_rec_file)
    pwd_ale_formatted_gnm_tree                          = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree)
    pwd_ale_formatted_gnm_tree_with_len                 = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree_with_len)
    pwd_ale_formatted_gnm_tree_with_len_prefixed        = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree_with_len_prefixed)
    pwd_gene_tree_itol_colorstrip_txt                   = '%s/%s'                                               % (ale_hgt_plot_dir, gene_tree_itol_colorstrip_txt)
    pwd_genome_tree_itol_colorstrip_txt                 = '%s/%s'                                               % (ale_hgt_plot_dir, genome_tree_itol_colorstrip_txt)

    # run ale_splitter
    ale_splitter(ale_uml_rec_file)

    # read in phylum color
    p_color_dict = dict()
    for each_line in open(p_color_txt):
        each_line_split = each_line.strip().split('\t')
        phylum_id = each_line_split[1]
        color_id = each_line_split[0]
        p_color_dict[phylum_id] = color_id

    internal_node_to_leaf_dict = dict()
    paired_donor_to_recipient_leaf_dict = dict()
    hgt_freq_dict = dict()
    if os.path.isfile(pwd_uts_file) is True:

        # write out ALE formatted genome tree
        renamed_genome_tree_str = open(pwd_uml_rec_file).readlines()[2].strip().split('\t')[1]
        with open(pwd_ale_formatted_gnm_tree, 'w') as ale_renamed_species_tree_handle:
            ale_renamed_species_tree_handle.write(renamed_genome_tree_str + '\n')

        internal_node_to_leaf_dict, paired_donor_to_recipient_leaf_dict, hgt_freq_dict = uts_to_itol_connections(pwd_genome_tree_file_subset_for_ale, pwd_ale_formatted_gnm_tree, interal_node_prefix, pwd_uts_file, hgt_freq_cutoff, ignore_leaf_hgt, ignore_vertical_hgt, donor_node_min_leaf_num, recipient_node_min_leaf_num, pwd_itol_connection_txt_all, dr_separator)
    else:
        print('%s: uTs file not found, you need to run ALE first!' % qualified_og)

    # combine_trees
    combine_trees(pwd_genome_tree_file_subset_for_ale, pwd_ale_formatted_gnm_tree, pwd_ale_formatted_gnm_tree_with_len)

    # prefix_internal_nodes of combined tree
    prefix_internal_nodes(pwd_ale_formatted_gnm_tree_with_len, interal_node_prefix, pwd_ale_formatted_gnm_tree_with_len_prefixed)

    # write out iTOL label file for gene and genome tree, also colorstrip for taxonomy
    pwd_itol_label_txt_handle  = open(pwd_itol_label_txt, 'w')
    pwd_itol_label_txt_handle.write('LABELS\nSEPARATOR TAB\n\nDATA\n')
    pwd_gene_tree_itol_label_txt_handle = open(pwd_gene_tree_itol_label_txt, 'w')
    pwd_gene_tree_itol_label_txt_handle.write('LABELS\nSEPARATOR TAB\n\nDATA\n')
    wrote_gnm_set = set()
    gene_to_p_dict = dict()
    genome_to_p_dict = dict()
    for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():
        gene_gnm = '_'.join(each_gene.split('_')[:-1])
        genome_name_for_ale = gene_gnm
        genome_name_for_ale = genome_name_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
        genome_with_taxon = gnm_pco_dict[gene_gnm]
        gene_to_p_dict[each_gene] = genome_with_taxon.split('__')[0]
        if gene_gnm not in wrote_gnm_set:
            genome_to_p_dict[genome_name_for_ale] = genome_with_taxon.split('__')[0]
            pwd_itol_label_txt_handle.write('%s\t%s\n' % (genome_name_for_ale, genome_with_taxon))
            wrote_gnm_set.add(gene_gnm)
        pwd_gene_tree_itol_label_txt_handle.write('%s\t%s_%s\n' % (each_gene, genome_with_taxon, each_gene.split('_')[-1]))
    pwd_itol_label_txt_handle.close()
    pwd_gene_tree_itol_label_txt_handle.close()

    iTOL(gene_to_p_dict, p_color_dict, pwd_gene_tree_itol_colorstrip_txt)
    iTOL(genome_to_p_dict, p_color_dict, pwd_genome_tree_itol_colorstrip_txt)

    # root gene tree at midpoint
    gene_tree_to_plot = pwd_gene_tree_treefile_subset
    if root_gene_tree_at_midpoint is True:
        root_at_midpoint(pwd_gene_tree_treefile_subset, pwd_gene_tree_treefile_subset_midpoint_rooted)
        gene_tree_to_plot = pwd_gene_tree_treefile_subset_midpoint_rooted

    # plot separately
    n = 1
    for each_d2r in paired_donor_to_recipient_leaf_dict:
        each_d2r_freq = hgt_freq_dict[each_d2r]
        each_d2r_d_list = paired_donor_to_recipient_leaf_dict[each_d2r][0]
        each_d2r_r_list = paired_donor_to_recipient_leaf_dict[each_d2r][1]
        pwd_itol_label_txt                                  = '%s/%s_iTOL_genome_pco.txt'               % (ale_op_dir, qualified_og)
        pwd_gene_tree_itol_label_txt                        = '%s/%s_iTOL_gene_pco.txt'                 % (ale_hgt_plot_dir, qualified_og)
        pwd_gnm_tree_label_color_txt                        = '%s/%s_iTOL_label_color_genome_%s.txt'    % (ale_hgt_plot_dir, qualified_og, each_d2r)
        pwd_gene_tree_label_color_txt                       = '%s/%s_iTOL_label_color_gene_%s.txt'      % (ale_hgt_plot_dir, qualified_og, each_d2r)
        pwd_itol_connection_txt                             = '%s/%s_iTOL_connection_%s.txt'            % (ale_hgt_plot_dir, qualified_og, each_d2r)
        pwd_ale_formatted_gnm_tree_with_len_prefixed_pdf    = '%s/%s_genome_tree_with_HGT_%s.pdf'       % (ale_wd, qualified_og, each_d2r)
        pwd_gene_tree_treefile_subset_pdf                   = '%s/%s_subset_%s.pdf'                     % (ale_hgt_plot_dir, qualified_og, each_d2r)
        pwd_gene_tree_treefile_subset_pdf_rooted            = '%s/%s_subset_%s_rooted.pdf'              % (ale_hgt_plot_dir, qualified_og, each_d2r)
        pwd_combined_image_with_ale_hgts                    = '%s/%s_HGT_%s_%s_%s.pdf'                  % (ale_hgt_plot_dir, qualified_og, n, each_d2r, each_d2r_freq)

        # write out gnm_tree_label_color_txt
        pwd_gnm_tree_label_color_txt_handle = open(pwd_gnm_tree_label_color_txt, 'w')
        pwd_gnm_tree_label_color_txt_handle.write('DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\texample_style\nCOLOR\t#ffff00\n\nDATA\n')
        pwd_gnm_tree_label_color_txt_handle.write('%s\tlabel\tclade\t%s\t1\tnormal\n' % (each_d2r.split(dr_separator)[0], d_color))
        pwd_gnm_tree_label_color_txt_handle.write('%s\tlabel\tclade\t%s\t1\tnormal\n' % (each_d2r.split(dr_separator)[1], r_color))
        pwd_gnm_tree_label_color_txt_handle.close()

        # write out iTOL label file for gene and genome tree, also colorstrip for taxonomy
        pwd_gene_tree_label_color_txt_handle = open(pwd_gene_tree_label_color_txt, 'w')
        pwd_gene_tree_label_color_txt_handle.write('DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\texample_style\nCOLOR\t#ffff00\n\nDATA\n')
        for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():

            gene_name_for_ale = '_'.join(each_gene.strip().split('_')[:-1])
            gene_name_for_ale = gene_name_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
            if gene_name_for_ale in each_d2r_d_list:
                pwd_gene_tree_label_color_txt_handle.write('%s\tlabel\tnode\t%s\t1\tnormal\n' % (each_gene, d_color))
            elif gene_name_for_ale in each_d2r_r_list:
                pwd_gene_tree_label_color_txt_handle.write('%s\tlabel\tnode\t%s\t1\tnormal\n' % (each_gene, r_color))
        pwd_gene_tree_label_color_txt_handle.close()

        itol_tree(pwd_ale_formatted_gnm_tree_with_len_prefixed, [pwd_gnm_tree_label_color_txt, pwd_itol_label_txt, pwd_itol_connection_txt, pwd_genome_tree_itol_colorstrip_txt], project_name, API_key, display_mode, pwd_ale_formatted_gnm_tree_with_len_prefixed_pdf)
        itol_tree(gene_tree_to_plot, [pwd_gene_tree_itol_label_txt, pwd_gene_tree_label_color_txt, pwd_gene_tree_itol_colorstrip_txt], project_name, API_key, display_mode, pwd_gene_tree_treefile_subset_pdf)
        merge_pdf(pwd_ale_formatted_gnm_tree_with_len_prefixed_pdf, pwd_gene_tree_treefile_subset_pdf, 66, pwd_combined_image_with_ale_hgts)
        n += 1

        os.system('mv %s %s/annotation_files/' % (pwd_ale_formatted_gnm_tree_with_len_prefixed_pdf, ale_hgt_plot_dir))
        os.system('mv %s %s/annotation_files/' % (pwd_gene_tree_treefile_subset_pdf, ale_hgt_plot_dir))
        os.system('mv %s %s/annotation_files/' % (pwd_gnm_tree_label_color_txt, ale_hgt_plot_dir))
        os.system('mv %s %s/annotation_files/' % (pwd_gene_tree_label_color_txt, ale_hgt_plot_dir))
        os.system('mv %s %s/annotation_files/' % (pwd_itol_connection_txt, ale_hgt_plot_dir))
    os.system('mv %s %s/annotation_files/' % (pwd_itol_label_txt, ale_hgt_plot_dir))
    os.system('mv %s %s/annotation_files/' % (pwd_gene_tree_itol_label_txt, ale_hgt_plot_dir))
    os.system('mv %s %s/annotation_files/' % (pwd_itol_connection_txt_all, ale_hgt_plot_dir))
    os.system('mv %s %s/annotation_files/' % (pwd_gene_tree_itol_colorstrip_txt, ale_hgt_plot_dir))
    os.system('mv %s %s/annotation_files/' % (pwd_genome_tree_itol_colorstrip_txt, ale_hgt_plot_dir))


def ale_splitter(rec_file):

    options = [True, True, True, True]
    with open(rec_file) as f:
        lines               = f.readlines()
        stree               = lines[2].strip()
        ll                  = lines[6].strip().split()[-1]
        rates               = lines[8].strip().split("\t")[1:]
        n_reconciled_trees  = int(lines[10].strip().split()[0])
        reconciled_trees    = lines[12:n_reconciled_trees + 12]
        n_of_events         = lines[12 + n_reconciled_trees + 1].split("\t")[1:]
        table               = lines[12 + n_reconciled_trees + 3:]

    if options[0]:
        with open(rec_file.replace("uml_rec", "stree"), "w") as f:
            f.write(stree.split("\t")[-1])
    if options[1]:
        with open(rec_file.replace("uml_rec", "info"), "w") as f:
            f.write("LL:" + "\t" + ll + "\n")
            f.write("Dp:" + "\t" + rates[0] + "\n")
            f.write("Tp:" + "\t" + rates[1] + "\n")
            f.write("Lp:" + "\t" + rates[2] + "\n")
            f.write("De:" + "\t" + n_of_events[0] + "\n")
            f.write("Te:" + "\t" + n_of_events[1] + "\n")
            f.write("Le:" + "\t" + n_of_events[2] + "\n")
            f.write("Se:" + "\t" + n_of_events[3] + "\n")
    if options[2]:
        with open(rec_file.replace("uml_rec", "recs"), "w") as f:
            for t in reconciled_trees:
                f.write(t)
    if options[3]:
        with open(rec_file.replace("uml_rec", "rec_table"), "w") as f:
            for e in table:
                f.write(e)


def ALE3(args):

    gene_tree_dir               = args['i1']
    ale_wd                      = args['i2']
    genome_taxon_txt            = args['c']
    ar_phylum_color_code_txt    = args['color']
    ale_hgt_plot_dir            = args['o']
    force_create_op_dir         = args['f']
    API_key                     = args['api']
    hgt_freq_cutoff             = args['fc']
    donor_node_min_leaf_num     = args['mld']
    recipient_node_min_leaf_num = args['mlr']
    project_name                = args['itol']

    ignore_vertical_hgt             = True                      # filter ALE predicted HGTs
    ignore_leaf_hgt                 = True                      # filter ALE predicted HGTs
    interal_node_prefix             = 'IN'                      # plot tree with HGT
    display_mode                    = '1'                       # plot tree with HGT, 1=rectangular, 2=circular, 3=unrooted
    align_leaf_name                 = True                      # plot tree with HGT
    show_scale                      = False                     # plot tree with HGT
    d_color                         = '#FF0000'                 # plot tree with HGT
    r_color                         = '#0000FF'                 # plot tree with HGT
    dr_separator                    = '_to_'                    # plot tree with HGT
    root_gene_tree_at_midpoint      = True                      # plot tree with HGT

    ####################################################################################################################

    ufboot_file_re   = '%s/*.ufboot' % gene_tree_dir
    ufboot_file_list = glob.glob(ufboot_file_re)
    og_to_process_list = []
    for each_ufboot in ufboot_file_list:
        _, ufboot_base, _ = sep_path_basename_ext(each_ufboot)
        og_to_process_list.append(ufboot_base)

    # read in genome taxonomy
    gnm_pco_dict = dict()
    for each_gnm in open(genome_taxon_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id     = each_gnm_split[0]
        taxon_str  = each_gnm_split[1]
        gnm_phylum = taxon_str.split(';')[1]
        gnm_class  = taxon_str.split(';')[2]
        gnm_order  = taxon_str.split(';')[3]
        gnm_pco_dict[gnm_id] = '%s__%s__%s__%s' % (gnm_phylum[3:], gnm_class[3:], gnm_order[3:], gnm_id)

    if os.path.isdir(ale_hgt_plot_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % ale_hgt_plot_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % ale_hgt_plot_dir)
    os.system('mkdir %s/annotation_files' % ale_hgt_plot_dir)

    # parse ALE output
    n = 1
    for qualified_og in og_to_process_list:

        print('%s (%s/%s): Parsing ALE outputs' % (qualified_og, n, len(og_to_process_list)))
        current_arg_list = [qualified_og, gene_tree_dir, ale_wd, ale_wd, ale_hgt_plot_dir, interal_node_prefix,
                            gnm_pco_dict, d_color, r_color, project_name, API_key, display_mode, hgt_freq_cutoff,
                            ignore_leaf_hgt, ignore_vertical_hgt, donor_node_min_leaf_num, recipient_node_min_leaf_num,
                            dr_separator, root_gene_tree_at_midpoint, ar_phylum_color_code_txt]
        parse_ale_op_worker(current_arg_list)
        n += 1

    print('Done!')


def ALE4(args):
    pass


if __name__ == '__main__':
    pass

    # ALE1_parser = argparse.ArgumentParser()
    # ALE1_parser.add_argument('-i',   required=True,                         help='orthologous groups, either from orthofinder or oma')
    # ALE1_parser.add_argument('-s',   required=True,                         help='sequence file, e.g., combined.faa')
    # ALE1_parser.add_argument('-p',   required=True,                         help='orthologous identification program, orthofinder or oma')
    # ALE1_parser.add_argument('-m',   required=False, type=int, default=50,  help='min_og_genome_num, default: 50')
    # ALE1_parser.add_argument('-n',   required=False, type=int, default=2,   help='min_og_phylum_num, default: 2')
    # ALE1_parser.add_argument('-o',   required=True,                         help='output dir, i.e., OMA working directory')
    # ALE1_parser.add_argument('-t',   required=False, type=int, default=6,   help='number of threads, default: 6')
    # ALE1_parser.add_argument('-jt',  required=False, type=int, default=3,   help='number of threads for job script, default: 3')
    # ALE1_parser.add_argument('-f',   required=False, action="store_true",   help='force overwrite')
    # ALE1_parser.add_argument('-c',   required=True,                         help='genome_taxon_txt')
    # args = vars(ALE1_parser.parse_args())
    # ALE1(args)

    # ALE2_parser = argparse.ArgumentParser()
    # ALE2_parser.add_argument('-i',   required=True,                         help='ALE1 output directory')
    # ALE2_parser.add_argument('-s',   required=True,                         help='rooted species tree')
    # ALE2_parser.add_argument('-c',   required=True,                         help='genome_taxon_txt')
    # ALE2_parser.add_argument('-o',   required=True,                         help='output dir, i.e., OMA working directory')
    # ALE2_parser.add_argument('-f',   required=False, action="store_true",   help='force overwrite')
    # ALE2_parser.add_argument('-t',   required=False, type=int, default=6,   help='number of threads, default: 6')
    # args = vars(ALE2_parser.parse_args())
    # ALE2(args)

    # ALE3_parser = argparse.ArgumentParser()
    # ALE3_parser.add_argument('-i1',     required=True,                          help='ALE1 output directory')
    # ALE3_parser.add_argument('-i2',     required=True,                          help='ALE2 output directory')
    # ALE3_parser.add_argument('-c',      required=True,                          help='genome_taxon_txt')
    # ALE3_parser.add_argument('-color',  required=True,                          help='phylum_color_code.txt')
    # ALE3_parser.add_argument('-o',      required=True,                          help='output dir, i.e., ALE3_op_dir')
    # ALE3_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    # args = vars(ALE3_parser.parse_args())
    # ALE3(args)
