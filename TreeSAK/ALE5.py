import os
import glob
import operator
from ete3 import Tree
from itertools import chain
from itertools import combinations


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def powerset(iterable):

    " powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3) "

    s = list(iterable)  # allows duplicate elements
    chain_obj = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    combo_lol = []
    for _, combo in enumerate(chain_obj, 1):
        if len(list(combo)) > 0:
            combo_lol.append(list(combo))

    return combo_lol


def lca_to_two_leaves(species_tree_from_ale, internal_node_id):

    # read in ale species tree
    stree_ale = Tree(species_tree_from_ale, format=1)

    # get all leaves of the internal node
    internal_node = stree_ale.search_nodes(name=internal_node_id)[0]
    internal_node_leaf_object = internal_node.get_leaves()
    internal_node_leaf_set = set()
    for each_leaf in internal_node_leaf_object:
        internal_node_leaf_set.add(each_leaf.name)

    # get the two leaves needed
    targeted_two_leaves = []
    leaves_found = False
    for leaf_1 in internal_node_leaf_set:
        for leaf_2 in internal_node_leaf_set:
            if leaf_1 != leaf_2:
                if leaves_found is False:
                    current_lca_id = stree_ale.get_common_ancestor(leaf_1, leaf_2).name
                    if current_lca_id == internal_node_id:
                        targeted_two_leaves.append(leaf_1)
                        targeted_two_leaves.append(leaf_2)
                        leaves_found = True

    return targeted_two_leaves[0], targeted_two_leaves[1]


def keep_highest_rrtc(rrtc_in, rrtc_out):

    rrtc_highest_prob_dict = dict()
    for each_rrtc in open(rrtc_in):
        rrtc_r = each_rrtc.strip().split(':')[0].split('\t')[0]
        rrtc_d = each_rrtc.strip().split(':')[0].split('\t')[1]
        rrtc_v = float(each_rrtc.strip().split(':')[1])
        rrtc_key = '%s___%s' % (rrtc_r, rrtc_d)
        if rrtc_key not in rrtc_highest_prob_dict:
            rrtc_highest_prob_dict[rrtc_key] = rrtc_v
        else:
            if rrtc_v > rrtc_highest_prob_dict[rrtc_key]:
                rrtc_highest_prob_dict[rrtc_key] = rrtc_v

    with open(rrtc_out, 'w') as rrtc_out_handle:
        for each_rrtc in sorted(rrtc_highest_prob_dict.items(), key=operator.itemgetter(1))[::-1]:
            rrtc_r = each_rrtc[0].split('___')[0]
            rrtc_d = each_rrtc[0].split('___')[1]
            rrtc_v = each_rrtc[1]
            rrtc_out_handle.write('%s\t%s:%s\n' % (rrtc_r, rrtc_d, rrtc_v))


########################################################################################################################

# file in
ip_dir                  = '/Users/songweizhi/Desktop/DateArTree/01_HGT_ALE_with_OMA'
species_tree_from_ale   = '/Users/songweizhi/Desktop/DateArTree/05_pRTC_wd/genome_tree.newick.ufboot.ale.stree'
round_list              = [1, 2, 3, 4, 5]
color_list              = ['dodgerblue', 'goldenrod1', 'darkorange1', 'seagreen3', 'orchid3']
min_detected_times      = 2

# file out
op_dir                  = '/Users/songweizhi/Desktop/DateArTree/05_pRTC_wd/HGTs_5_rds_ALE'

########################################################################################################################

rscript                 = '%s/rscript.R'                        % op_dir
plot_file               = '%s/Venn.pdf'                         % op_dir
rtc_txt                 = '%s/rrtc.txt'                         % op_dir
rtc_txt_highest         = '%s/rrtc_uniq_by_highest_prob.txt'    % op_dir

if os.path.isdir(op_dir):
    os.system('rm -r %s' % op_dir)
os.system('mkdir %s' % op_dir)

########################################################################################################################

hgt_dict = dict()
rd_to_hgt_dict= dict()
for each_rd in round_list:

    rd_id               = each_rd
    current_rd_op_dir   = '%s/ALE4_op_dir_%s_0.3'   % (ip_dir, each_rd)
    pdf_file_re         = '%s/*.%s'                 % (current_rd_op_dir, 'pdf')
    pdf_file_list       = glob.glob(pdf_file_re)

    rd_to_hgt_dict[rd_id] = set()

    for each_pdf in pdf_file_list:
        f_path, f_base, f_ext = sep_path_basename_ext(each_pdf)
        f_base_split = f_base.split('_')
        id_by_d_to_r = '%s_to_%s' % (f_base_split[3], f_base_split[5])
        rd_og        = '%s_%s'    % (each_rd, f_base_split[0])
        rd_og_value  = '%s_%s_%s' % (each_rd, f_base_split[0], f_base_split[6])

        rd_to_hgt_dict[rd_id].add(id_by_d_to_r)

        if id_by_d_to_r not in hgt_dict:
            hgt_dict[id_by_d_to_r] = []
        hgt_dict[id_by_d_to_r].append(rd_og_value)

################################################### get Venn diagram ###################################################

combination_list = powerset(round_list)

value_str_list = []
for each_cmbo in combination_list:
    current_str = ''
    if len(each_cmbo) == 1:
        current_value = rd_to_hgt_dict[each_cmbo[0]]
        current_str = 'area%s=%s' % (each_cmbo[0], len(current_value))
        value_str_list.append(current_str)
    else:
        value_lol = []
        for each_element in each_cmbo:
            ele_value = rd_to_hgt_dict[each_element]
            value_lol.append(ele_value)
        shared = set(value_lol[0]).intersection(*value_lol)
        current_str = 'n%s=%s' % (''.join([str(i) for i in each_cmbo]), len(shared))
        value_str_list.append(current_str)

value_str     = ', '.join(value_str_list)
label_str     = '"' + '", "'.join([str(i) for i in round_list]) + '"'
color_str     = '"' + '", "'.join([str(i) for i in color_list]) + '"'
font_size_str = ', '.join(['1.2']*len(combination_list))

rscript_handle = open(rscript, 'w')
rscript_handle.write('library(futile.logger)\n')
rscript_handle.write('library(gridBase)\n')
rscript_handle.write('library(VennDiagram)\n')
rscript_handle.write('pdf(file="%s")\n' % plot_file)
rscript_handle.write('venn.plot <- draw.quintuple.venn(%s, category=c(%s), fill=c(%s), cat.col=c(%s), cat.cex=1.2, cat.dist=0.3, margin=0.05, cex=c(%s), ind=TRUE)\n' % (value_str, label_str, color_str, color_str, font_size_str))
rscript_handle.write('dev.off()\n')
rscript_handle.close()

os.system('Rscript %s' % rscript)

########################################################################################################################

rtc_txt_handle = open(rtc_txt, 'w')
qualified_hgt_num = 0
for each_hgt in hgt_dict:

    occurence_list = hgt_dict[each_hgt]
    pdf_dir = '%s/%s_%s' % (op_dir, each_hgt, len(occurence_list))
    if len(occurence_list) >= min_detected_times:

        #################### prepare rtc file ####################

        donor_id           = each_hgt.split('_to_')[0][2:]
        recipient_id       = each_hgt.split('_to_')[1][2:]
        d_leaf_1, d_leaf_2 = lca_to_two_leaves(species_tree_from_ale, donor_id)
        r_leaf_1, r_leaf_2 = lca_to_two_leaves(species_tree_from_ale, recipient_id)

        for each_occurence in occurence_list:
            value = each_occurence.split('_')[-1]
            rtc_str = '%s,%s\t%s,%s:%s' % (r_leaf_1, r_leaf_2, d_leaf_1, d_leaf_2, value)
            rtc_txt_handle.write(rtc_str + '\n')

        ##########################################################

        qualified_hgt_num += 1
        os.system('mkdir %s' % pdf_dir)
        for each_h in occurence_list:
            rd_id             = each_h.split('_')[0]
            og_id             = each_h.split('_')[1]
            value             = each_h.split('_')[2]
            pwd_input_pdf_in  = '%s/ALE4_op_dir_%s_0.3/%s_HGT_*_%s_%s.pdf'  % (ip_dir, rd_id, og_id, each_hgt,value)
            pwd_input_pdf_out = '%s/%s_%s_%s_%s.pdf'                        % (pdf_dir, rd_id, og_id, each_hgt,value)
            os.system('cp %s %s' % (pwd_input_pdf_in, pwd_input_pdf_out))
rtc_txt_handle.close()

# remove redundant HGTs, keep the one with the highest probability
keep_highest_rrtc(rtc_txt, rtc_txt_highest)

print('The number of HGTs detected in >= %s runs is %s.' % (min_detected_times, qualified_hgt_num))
print('Done!')
