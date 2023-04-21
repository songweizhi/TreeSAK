import os
from PIL import Image
from ete3 import Tree
from ete3 import TextFace, TreeStyle, NodeStyle


def subset_tree(tree_file_in, leaves_to_keep_list, tree_file_out):

    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(leaves_to_keep_list, preserve_branch_length=True)
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


########################################################################################################################

# file in
species_tree_file      = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/concatenated_rooted_subset.treefile'
gene_tree_file         = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/OG0001758.treefile'
taxon_for_MetaCHIP_txt = '/Users/songweizhi/Desktop/DateArTree/0_HGT_MetaCHIP/taxon_for_MetaCHIP.txt'
hgt_txt                = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/OG0001758_HGTs.txt'
align_leaf_name        = True
show_scale             = False

# file out
species_tree_plot      = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/concatenated_rooted_subset.png'
gene_tree_plot         = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/OG0001758.png'
combined_image         = '/Users/songweizhi/Desktop/reconcile_gene_genome_tree/OG0001758_species_gene_tree.png'

########################################################################################################################

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

# get gene tree leaf name dict
leaf_name_dict = dict()
for each_gene in Tree(gene_tree_file).get_leaf_names():
    gene_id = each_gene
    gene_id_new = gene_id
    if '.gtdb' in gene_id_new:
        gene_id_new = gene_id_new.replace('.gtdb', '')
    gene_genome             = '_'.join(gene_id_new.split('_')[:-1])
    genome_pco              = gnm_pco_dict[gene_genome]
    gene_id_with_taxon      = '%s_%s' % (genome_pco, gene_id.split('_')[-1])
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

plot_tree(species_tree_file, 'Species tree (rooted)', gnm_pco_dict,   species_tree_leaf_color_dict, align_leaf_name, show_scale, species_tree_plot)
plot_tree(gene_tree_file,    'Gene tree (unrooted)',  leaf_name_dict, gene_tree_leaf_color_dict,    align_leaf_name, show_scale, gene_tree_plot)
merge_image([species_tree_plot, gene_tree_plot], combined_image)
os.system('rm %s' % species_tree_plot)
os.system('rm %s' % gene_tree_plot)

