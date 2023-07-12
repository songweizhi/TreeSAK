import os
from PIL import Image
from ete3 import Tree
from itolapi import Itol
from ete3 import TextFace, TreeStyle, NodeStyle


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


def combine_images(image_file_list, output_image):

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


def uts_to_itol_connections(genome_tree_file, uts_file, freq_cutoff, ignore_leaf_hgt, itol_connection_txt):
    leaf_id_set = [i.name for i in Tree(genome_tree_file).get_leaves()]
    with open(itol_connection_txt, 'w') as itol_connection_txt_handle:
        itol_connection_txt_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
        itol_connection_txt_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t1\nARROW_SIZE\t20\nLOOP_SIZE\t100\n')
        itol_connection_txt_handle.write('MAXIMUM_LINE_WIDTH\t10\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t0\nDATA\n')
        for each_line in open(uts_file):
            if not each_line.startswith('#'):
                each_line_split = each_line.strip().split('\t')
                donor = each_line_split[0]
                recipient = each_line_split[1]

                if donor.isdigit() is True:
                    donor = 'I' + donor

                if recipient.isdigit() is True:
                    recipient = 'I' + recipient

                freq = float(each_line_split[2])
                if freq >= freq_cutoff:

                    if ignore_leaf_hgt is False:
                        itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor, recipient, freq, '#EB984E', 'normal', donor, recipient, freq))
                    else:
                        if (each_line_split[0] not in leaf_id_set) and (each_line_split[1] not in leaf_id_set):
                            itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor, recipient, freq, '#EB984E', 'normal', donor, recipient, freq))


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
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', '0')
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '0')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


########################################################################################################################

'''
This script was write to extract ortho groups with MetaCHIP identified HGTs.
It also prepares job script for running mafft-einsi, trimal and iqtree.
'''

# file in
wd                              = '/Users/songweizhi/Desktop/DateArTree/0_HGT_MetaCHIP'
genome_tree_file                = '%s/concatenated.treefile'                                            % wd
outgroup                        = '%s/out_group.txt'                                                    % wd
taxon_for_MetaCHIP_txt          = '%s/taxon_for_MetaCHIP.txt'                                           % wd
orthogroups_op_txt              = '%s/Orthogroups.txt'                                                  % wd
metachip_BM_op_txt_list         = [('%s/Archaea_133_HGT_pc_p20_HGTs_BM.txt' % wd), ('%s/Archaea_133_HGT_pc_c52_HGTs_BM.txt' % wd)]
combined_faa                    = '%s/Archaea_133_HGT_pc_pc_combined_faa.fasta'                         % wd
ale_splitter_py                 = '/home-user/wzsong/Tests/ALE/ALEtutorial/ale_splitter_modified.py'
min_hgt_gene_per_og             = 6
extract_sequence                = False
align_leaf_name                 = True
show_scale                      = False
js_num_threads                  = 6
hgt_freq_cutoff                 = 0.3
force_create_png_folder         = True
force_create_ale_wd             = True
ignore_leaf_hgt                 = True
designate_ogs                   = ['OG0000032', 'OG0000102', 'OG0000143', 'OG0000168', 'OG0000297', 'OG0000302', 'OG0001758', 'OG0002466']
project_name                    = 'batch_access_tmp'
API_key                         = 'S1kZZuDHc0d5M7J5vLnUNQ'
display_mode                    = '1'  # # 1=rectangular, 2=circular, 3=unrooted

# file out
op_dir                          = '/Users/songweizhi/Desktop/DateArTree/0_HGT_MetaCHIP/op_dir'
png_dir                         = '%s/0_png_dir'                                                        % op_dir
ale_wd                          = '%s/0_ALE_wd'                                                         % op_dir
ale_op_dir                      = '%s/0_ALE_op_dir'                                                     % op_dir
genome_tree_file_rooted         = '%s/concatenated_rooted.treefile'                                     % op_dir
summary_txt                     = '%s/OGs_with_HGT_stats.txt'                                           % op_dir
iq_tree_cmd_txt                 = '%s/iqtree_cmds.txt'                                                  % op_dir

########################################################################################################################

if force_create_png_folder is True:
    if os.path.isdir(png_dir) is True:
        os.system('rm -r %s' % png_dir)
os.system('mkdir %s' % png_dir)

if force_create_ale_wd is True:
    if os.path.isdir(ale_wd) is True:
        os.system('rm -r %s' % ale_wd)
os.system('mkdir %s' % ale_wd)

# read in genome taxonomy
gnm_p_dict = dict()
gnm_c_dict = dict()
gnm_o_dict = dict()
gnm_pco_dict = dict()
for each_gnm in open(taxon_for_MetaCHIP_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id     = each_gnm_split[0]
    taxon_str  = each_gnm_split[1]
    gnm_phylum = taxon_str.split(';')[1]
    gnm_class  = taxon_str.split(';')[2]
    gnm_order  = taxon_str.split(';')[3]
    gnm_id_new = gnm_id
    if '.gtdb' in gnm_id_new:
        gnm_id_new = gnm_id_new.replace('.gtdb', '')
    gnm_p_dict[gnm_id_new] = gnm_phylum
    gnm_c_dict[gnm_id_new] = gnm_class
    gnm_o_dict[gnm_id_new] = gnm_order
    gnm_pco_dict[gnm_id_new] = '%s__%s__%s__%s' % (gnm_phylum[3:], gnm_class[3:], gnm_order[3:], gnm_id_new)

# root genome tree with outgroup
root_with_out_group(genome_tree_file, outgroup, genome_tree_file_rooted)

# read in genome taxonomy
gnm_class_dict = dict()
gnm_phylum_dict = dict()
for each_gnm in open(taxon_for_MetaCHIP_txt):
    each_gnm_split          = each_gnm.strip().split('\t')
    gnm_id                  = each_gnm_split[0]
    taxon_str               = each_gnm_split[1]
    gnm_phylum              = taxon_str.split(';')[1]
    gnm_class               = taxon_str.split(';')[2]
    gnm_class_dict[gnm_id]  = gnm_class
    gnm_phylum_dict[gnm_id] = gnm_phylum

# read in OrthoFinder output
ortho_to_gene_dict = dict()
gene_to_ortho_dict = dict()
for each_og in open(orthogroups_op_txt):
    each_og_split = each_og.strip().split(' ')
    og_id = each_og_split[0][:-1]
    gene_list = each_og_split[1:]
    ortho_to_gene_dict[og_id] = gene_list
    for each_gene in gene_list:
        gene_to_ortho_dict[each_gene] = og_id

# read in MetaCHIP op
ogs_with_hgt = set()
og_to_hgt_dict = dict()
hgt_candidate_set = set()
og_to_paired_hgt_dict = dict()
for each_BM_op_txt in metachip_BM_op_txt_list:
    for each_hgt in open(each_BM_op_txt):
        if not each_hgt.startswith('Gene_1\tGene_2'):
            each_hgt_split = each_hgt.strip().split('\t')
            gene_1         = each_hgt_split[0]
            gene_2         = each_hgt_split[1]
            gnm_1          = '_'.join(gene_1.split('_')[:-1])
            gnm_2          = '_'.join(gene_2.split('_')[:-1])
            gnm_1_phylum   = gnm_phylum_dict[gnm_1]
            gnm_2_phylum   = gnm_phylum_dict[gnm_2]
            gnm_1_class    = gnm_class_dict[gnm_1]
            gnm_2_class    = gnm_class_dict[gnm_2]
            gene_1_og      = gene_to_ortho_dict.get(gene_1, '')
            gene_2_og      = gene_to_ortho_dict.get(gene_2, '')
            if (gene_1_og == gene_2_og) and (gene_1_og != ''):
                ogs_with_hgt.add(gene_1_og)
                hgt_candidate_set.add(gene_1)
                hgt_candidate_set.add(gene_2)
                if gene_1_og not in og_to_hgt_dict:
                    og_to_hgt_dict[gene_1_og] = {gene_1, gene_2}
                    og_to_paired_hgt_dict[gene_1_og] = [[gene_1, gene_2]]
                else:
                    og_to_hgt_dict[gene_1_og].add(gene_1)
                    og_to_hgt_dict[gene_1_og].add(gene_2)
                    og_to_paired_hgt_dict[gene_1_og].append([gene_1, gene_2])

# sort og by the number of identified HGTs
og_to_hgt_num_dict = {i: len(og_to_hgt_dict[i]) for i in og_to_hgt_dict}
og_to_hgt_num_dict_sorted = {k: v for k, v in sorted(og_to_hgt_num_dict.items(), key=lambda item: item[1])[::-1]}

# summarize at gene, genome, class and phylum levels
summary_txt_handle = open(summary_txt, 'w')
summary_txt_handle.write('Index\tOG\tTotal_Gene\tHGT_Gene\tHGT_Genome\tHGT_Class\tHGT_Phylum\n')
qualified_og_set = set()
index = 1
qualified_og_num = 0
for each_hgt_og in og_to_hgt_num_dict_sorted:

    current_og_hgt_txt          = '%s/%s_HGTs.txt'          % (op_dir, each_hgt_og)
    current_og_paired_hgt_txt   = '%s/%s_paired_HGTs.txt'   % (op_dir, each_hgt_og)
    current_og_gene_set         = ortho_to_gene_dict[each_hgt_og]
    current_og_hgt_set          = og_to_hgt_dict[each_hgt_og]
    current_og_paired_hgt_list   = og_to_paired_hgt_dict[each_hgt_og]
    current_og_gnm_set = {'_'.join(i.split('_')[:-1]) for i in current_og_hgt_set}
    current_og_c_set = {gnm_class_dict[i] for i in current_og_gnm_set}
    current_og_p_set = {gnm_phylum_dict[i] for i in current_og_gnm_set}

    with open(current_og_hgt_txt, 'w') as current_og_hgt_txt_handle:
        current_og_hgt_txt_handle.write('%s\n' % '\n'.join(current_og_hgt_set))

    with open(current_og_paired_hgt_txt, 'w') as current_og_paired_hgt_txt_handle:
        for each_pair in current_og_paired_hgt_list:
            current_og_paired_hgt_txt_handle.write('%s\n' % '\t'.join(each_pair))

    if len(current_og_hgt_set) >= min_hgt_gene_per_og:
        summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (index, each_hgt_og, len(current_og_gene_set), len(current_og_hgt_set), len(current_og_gnm_set), len(current_og_c_set), len(current_og_p_set)))
        qualified_og_set.add(each_hgt_og)
        index += 1
        qualified_og_num += 1
summary_txt_handle.close()
print('The number of orthogroups with >=%s genes involved in HGT: %s' % (min_hgt_gene_per_og, qualified_og_num))

# process qualified OG
og_to_process = qualified_og_set
if len(designate_ogs) > 0:
    og_to_process = designate_ogs

for qualified_og in og_to_process:
    qualified_og_gene_set         = ortho_to_gene_dict[qualified_og]
    qualified_og_gene_txt         = '%s/%s.txt'           % (op_dir, qualified_og)
    qualified_og_gene_faa         = '%s/%s.faa'           % (op_dir, qualified_og)
    qualified_og_gene_aln         = '%s/%s.aln'           % (op_dir, qualified_og)
    qualified_og_gene_aln_trimmed = '%s/%s_trimmed.aln'   % (op_dir, qualified_og)
    js_file                       = '%s/zjs_%s.sh'        % (op_dir, qualified_og)

    with open(qualified_og_gene_txt, 'w') as qualified_og_gene_txt_handle:
        qualified_og_gene_txt_handle.write('\n'.join(qualified_og_gene_set))

    # extract sequences
    extract_seq_cmd = 'BioSAK select_seq -seq %s -id %s -out %s -option 1 -oneline' % (combined_faa, qualified_og_gene_txt, qualified_og_gene_faa)
    if extract_sequence is True:
        os.system(extract_seq_cmd)

    # write out js for mafft, trimal and iqtree
    mafft_cmd          = 'mafft-einsi --thread %s --quiet %s > %s'                        % (js_num_threads, qualified_og_gene_faa, qualified_og_gene_aln)
    trimal_cmd         = 'trimal -in %s -out %s -automated1'                              % (qualified_og_gene_aln, qualified_og_gene_aln_trimmed)
    iqtree_cmd         = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'          % (js_num_threads, qualified_og_gene_aln, qualified_og)
    iqtree_cmd_trimmed = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s_trimmed'  % (js_num_threads, qualified_og_gene_aln_trimmed, qualified_og)
    js_file_handle = open(js_file, 'w')
    js_file_handle.write('#!/bin/bash\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task %s\n' % js_num_threads)
    js_file_handle.write(mafft_cmd.replace((op_dir + '/'), '') + '\n')
    js_file_handle.write(trimal_cmd.replace((op_dir + '/'), '') + '\n')
    js_file_handle.write(iqtree_cmd.replace((op_dir + '/'), '') + '\n')
    js_file_handle.write(iqtree_cmd_trimmed.replace((op_dir + '/'), '') + '\n')
    js_file_handle.close()

    # prepare input files for running ALE
    pwd_gene_tree_ufboot                = '%s/%s/%s.ufboot'                         % (op_dir, qualified_og, qualified_og)
    pwd_gene_tree_ufboot_updated        = '%s/%s/%s.ufboot.renamed'                 % (op_dir, qualified_og, qualified_og)
    pwd_gene_tree_treefile              = '%s/%s/%s.treefile'                       % (op_dir, qualified_og, qualified_og)
    pwd_gene_tree_treefile_updated      = '%s/%s/%s.treefile.renamed'               % (op_dir, qualified_og, qualified_og)
    iTOL_binary                         = '%s/%s/%s_iTOL_binary.txt'                % (op_dir, qualified_og, qualified_og)
    genomes_to_keep_txt                 = '%s/%s/%s_genomes_to_keep.txt'            % (op_dir, qualified_og, qualified_og)
    genome_tree_file_subset             = '%s/%s/%s_genome_tree_subset.treefile'    % (op_dir, qualified_og, qualified_og)
    hgt_txt                             = '%s/%s_HGTs.txt'                          % (op_dir, qualified_og)
    gene_tree_updated_plot              = '%s/%s_updated.png'                       % (png_dir, qualified_og)
    genome_tree_subset_plot             = '%s/%s_genome_tree_subset.png'            % (png_dir, qualified_og)
    combined_image                      = '%s/%s_combined_trees.png'                % (png_dir, qualified_og)
    genome_tree_file_subset_for_ALE     = '%s_genome_tree_for_ALE.treefile'         % qualified_og
    genome_tree_file_subset_for_ALE_plot= '%s_genome_tree_for_ALE.pdf'              % qualified_og
    gene_tree_updated_for_ALE           = '%s_gene_tree_for_ALE.ufboot'             % qualified_og
    pwd_genome_tree_file_subset_for_ALE = '%s/%s'                                   % (ale_wd, genome_tree_file_subset_for_ALE)
    pwd_genome_tree_file_subset_for_ALE_plot = '%s/%s'                                   % (ale_op_dir, genome_tree_file_subset_for_ALE_plot)
    pwd_gene_tree_updated_for_ALE       = '%s/%s'                                   % (ale_wd, gene_tree_updated_for_ALE)
    js_ale                              = '%s/js_%s_ALE.sh'                         % (ale_wd, qualified_og)
    itol_connection_txt                 = '%s/%s_iTOL_connection.txt'               % (ale_op_dir, qualified_og)
    uts_file                            = '%s/%s.ale.uTs'                           % (ale_op_dir, gene_tree_updated_for_ALE)
    uml_rec_file                        = '%s/%s.ale.uml_rec'                       % (ale_op_dir, gene_tree_updated_for_ALE)
    ale_renamed_species_tree            = '%s/%s_ALE_renamed_genome_tree.tree'       % (ale_op_dir, gene_tree_updated_for_ALE)

    if os.path.isfile(pwd_gene_tree_ufboot) is False:
        print('%s not found!' % pwd_gene_tree_ufboot)
    else:
        print('Processing %s' % qualified_og)

        # get genome of leaves in gene tree
        gene_gnm_set = set()
        gnm_to_gene_dict = dict()
        gene_tree_str = open(pwd_gene_tree_ufboot).readline().strip()
        for each_gene in Tree(gene_tree_str).get_leaf_names():
            gene_gnm = each_gene.split('.gtdb')[0]
            gene_gnm_set.add(gene_gnm)
            if gene_gnm not in gnm_to_gene_dict:
                gnm_to_gene_dict[gene_gnm] = {each_gene}
            else:
                gnm_to_gene_dict[gene_gnm].add(each_gene)

        # subset genome tree
        genome_tree_leaf_set = Tree(genome_tree_file_rooted).get_leaf_names()
        gene_gnm_set_in_genomes_tree = set(genome_tree_leaf_set).intersection(gene_gnm_set)
        subset_tree(genome_tree_file_rooted, gene_gnm_set_in_genomes_tree, genome_tree_file_subset)

        genome_tree_str_for_ale = subset_tree(genome_tree_file_rooted, gene_gnm_set_in_genomes_tree, None)
        genome_tree_str_for_ale = genome_tree_str_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
        with open(pwd_genome_tree_file_subset_for_ALE, 'w') as genome_tree_file_subset_for_ALE_handle:
            genome_tree_file_subset_for_ALE_handle.write(genome_tree_str_for_ale)

        # get genes to keep in gene tree
        gene_set_to_keep = set()
        for each_gnm in gene_gnm_set_in_genomes_tree:
            gene_set_to_keep.update(gnm_to_gene_dict.get(each_gnm, set()))

        gene_tree_str_updated_for_ALE = subset_tree(pwd_gene_tree_treefile, gene_set_to_keep, None).replace('GCA_', 'GCA').replace('GCF_', 'GCF').replace('.gtdb', '')
        with open(pwd_gene_tree_updated_for_ALE, 'w') as pwd_gene_tree_updated_for_ALE_handle:
            pwd_gene_tree_updated_for_ALE_handle.write(gene_tree_str_updated_for_ALE + '\n')

        # subset gene tree
        pwd_gene_tree_updated_handle = open(pwd_gene_tree_ufboot_updated, 'w')
        for each_gene_tree in open(pwd_gene_tree_ufboot):
            gene_tree_str = each_gene_tree.strip()
            gene_tree_str_updated = subset_tree(gene_tree_str, gene_set_to_keep, None)
            pwd_gene_tree_updated_handle.write(gene_tree_str_updated + '\n')
        pwd_gene_tree_updated_handle.close()

        # get gene tree leaf name dict (for plot)
        leaf_name_dict = dict()
        for each_gene in Tree(pwd_gene_tree_ufboot_updated).get_leaf_names():
            gene_id = each_gene
            gene_id_new = gene_id
            if '.gtdb' in gene_id_new:
                gene_id_new = gene_id_new.replace('.gtdb', '')
            gene_genome = '_'.join(gene_id_new.split('_')[:-1])
            genome_pco = gnm_pco_dict[gene_genome]
            gene_id_with_taxon = '%s_%s' % (genome_pco, gene_id.split('_')[-1])
            leaf_name_dict[gene_id] = gene_id_with_taxon

        # get gene tree leaf color dict
        gene_tree_leaf_color_dict = dict()
        species_tree_leaf_color_dict = dict()
        for each_hgt in open(hgt_txt):
            hgt_id = each_hgt.strip()
            hgt_gnm = '_'.join(hgt_id.split('_')[:-1])
            hgt_gnm_new = hgt_gnm
            if '.gtdb' in hgt_gnm_new:
                hgt_gnm_new = hgt_gnm_new.replace('.gtdb', '')
            gene_tree_leaf_color_dict[hgt_id] = 'red'
            species_tree_leaf_color_dict[hgt_gnm_new] = 'red'

        # plot tree
        plot_tree(genome_tree_file_subset, 'Species tree (rooted)', gnm_pco_dict, species_tree_leaf_color_dict, align_leaf_name, show_scale, genome_tree_subset_plot)
        plot_tree(pwd_gene_tree_treefile, 'Gene tree (unrooted)', leaf_name_dict, gene_tree_leaf_color_dict, align_leaf_name, show_scale, gene_tree_updated_plot)
        combine_images([genome_tree_subset_plot, gene_tree_updated_plot], combined_image)
        os.system('rm %s' % genome_tree_subset_plot)
        os.system('rm %s' % gene_tree_updated_plot)

        # run ALE
        with open(js_ale, 'w') as js_ale_handle:
            obtain_ale_file_cmd = 'ALEobserve %s'                       % gene_tree_updated_for_ALE
            reconciliation_cmd  = 'ALEml_undated %s %s.ale'             % (genome_tree_file_subset_for_ALE, gene_tree_updated_for_ALE)
            ale_splitter_cmd    = 'python3 %s -i %s.ale.uml_rec -sftr'  % (ale_splitter_py, gene_tree_updated_for_ALE)
            js_ale_handle.write('#!/bin/bash\n')
            js_ale_handle.write(obtain_ale_file_cmd + '\n')
            js_ale_handle.write(reconciliation_cmd + '\n')
            js_ale_handle.write(ale_splitter_cmd + '\n')

        # uts_to_itol_connections
        if os.path.isfile(uts_file) is True:
            uts_to_itol_connections(pwd_genome_tree_file_subset_for_ALE, uts_file, hgt_freq_cutoff, ignore_leaf_hgt, itol_connection_txt)

            # get ALE renamed genome tree
            uml_rec_file = '%s/%s.ale.uml_rec' % (ale_op_dir, gene_tree_updated_for_ALE)
            ale_renamed_species_tree = '%s/%s_ALE_renamed_genome_tree.nwk' % (ale_op_dir, gene_tree_updated_for_ALE)
            renamed_genome_tree_str = open(uml_rec_file).readlines()[2].strip().split('\t')[1]
            with open(ale_renamed_species_tree, 'w') as ale_renamed_species_tree_handle:
                ale_renamed_species_tree_handle.write(renamed_genome_tree_str + '\n')

            itol_tree(ale_renamed_species_tree, [itol_connection_txt], project_name, API_key, display_mode, pwd_genome_tree_file_subset_for_ALE_plot)

        else:
            print('uTs file not found, you need to run ALE first!')

