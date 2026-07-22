import os
import glob
import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


ALE9_usage = '''
=============================== ALE9 example commands ===============================

TreeSAK ALE9 -1 ALE1_op_dir -8 ALE8_op_dir -o ALE9_op_dir -cog AOA_genomes_faa_arCOG_wd -kegg AOA_genomes_faa_KEGG_wd -cog_db /Users/songweizhi/DB/arCOG18 -f 

=====================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def get_bar_plot_vertical(num_list_1, num_list_1_label, num_list_2, num_list_2_label, label_list, title_str, fig_width, fig_height, output_plot):

    x = np.arange(len(label_list))  # the label locations
    width = 0.3                     # the width of the bars

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.bar(x - width/2, num_list_1, color='#90EE90', width=width, label=num_list_1_label)
    ax.bar(x + width/2, num_list_2, color='#F87C63', width=width, label=num_list_2_label)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Gene number')
    ax.set_title(title_str)
    ax.set_xticks(x)
    ax.set_xticklabels(label_list, rotation=270, ha='left', rotation_mode='anchor')
    ax.legend()
    fig.tight_layout()
    plt.savefig(output_plot)
    plt.close()


def get_bar_plot(num_list_1, num_list_1_label, num_list_2, num_list_2_label, label_list, title_str, output_plot):

    num_list_1 = num_list_1[::-1]
    num_list_2 = num_list_2[::-1]
    label_list = label_list[::-1]

    y = np.arange(len(label_list))  # The label locations
    height = 0.3  # The height of the bars
    fig, ax = plt.subplots(figsize=(8, (len(label_list)*0.3)))

    # Plot horizontal bars
    ax.barh(y - height / 2, num_list_1, color='#90EE90', height=height, label=num_list_1_label)
    ax.barh(y + height / 2, num_list_2, color='#F87C63', height=height, label=num_list_2_label)

    # Add text for labels and title
    ax.set_xlabel('Gene number')
    ax.set_title(title_str)

    # Configure the Y-axis ticks on the right side
    ax.set_yticks(y)
    ax.set_yticklabels(label_list)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    # CRITICAL: Invert the X-axis to align bars to the right and put 0 on the right
    ax.invert_xaxis()

    ax.legend()
    fig.tight_layout()
    plt.savefig(output_plot)
    plt.close()

def get_og_fun(faa_dir, faa_ext, cog_annotation_dir, kegg_annotation_dir, op_dir):

    faa_file_re   = '%s/*.%s' % (faa_dir, faa_ext)
    faa_file_list = glob.glob(faa_file_re)
    og_to_gene_dict = dict()
    for each_faa_file in faa_file_list:
        _, _, f_base, _ = sep_path_basename_ext(each_faa_file)
        og_to_gene_dict[f_base] = set()
        for each_seq in SeqIO.parse(each_faa_file, 'fasta'):
            seq_id = each_seq.id
            og_to_gene_dict[f_base].add(seq_id)

    query_to_cog_txt_re = '%s/*COG_wd/*_query_to_cog.txt' % cog_annotation_dir
    query_to_cog_txt_list = glob.glob(query_to_cog_txt_re)
    gene_to_cog_dict = dict()
    cog_cate_dict = dict()
    cog_desc_dict = dict()
    for each_cog_file in query_to_cog_txt_list:
        for each_line in open(each_cog_file):
            if not each_line.startswith('Query\t'):
                each_line_split = each_line.strip().split('\t')
                if len(each_line_split) == 4:
                    gene_id  = each_line_split[0]
                    cog_id   = each_line_split[1]
                    cog_cate = each_line_split[2]
                    cog_desc = each_line_split[3]
                    gene_to_cog_dict[gene_id] = cog_id
                    cog_cate_dict[cog_id] = cog_cate
                    cog_desc_dict[cog_id] = cog_desc

    ko_assignment_abcd_txt_re   = '%s/*KEGG_wd/*_ko_assignment_ABCD.txt' % kegg_annotation_dir
    ko_assignment_abcd_txt_list = glob.glob(ko_assignment_abcd_txt_re)
    gene_to_kegg_dict = dict()
    ko_desc_dict      = dict()
    ko_d_to_b_dict    = dict()
    ko_d_to_c_dict    = dict()
    for each_kegg_file in ko_assignment_abcd_txt_list:
        for each_line in open(each_kegg_file):
            if not each_line.startswith('Gene_id\t'):
                each_line_split = each_line.strip().split('\t')
                if len(each_line_split) == 9:
                    gene_id   = each_line_split[0]
                    ko_id_b   = each_line_split[2][2:]
                    ko_id_c   = each_line_split[3][2:]
                    ko_id_d   = each_line_split[4][2:]
                    ko_desc_b = each_line_split[6]
                    ko_desc_c = each_line_split[7]
                    ko_desc_d = each_line_split[8]
                    gene_to_kegg_dict[gene_id] = ko_id_d
                    ko_desc_dict[ko_id_b]   = ko_desc_b
                    ko_desc_dict[ko_id_c]   = ko_desc_c
                    ko_desc_dict[ko_id_d]   = ko_desc_d
                    ko_d_to_b_dict[ko_id_d] = ko_id_b
                    ko_d_to_c_dict[ko_id_d] = ko_id_c

    og_to_cog_dict = dict()
    og_to_kegg_dict = dict()
    for each_og in og_to_gene_dict:
        gene_member = og_to_gene_dict[each_og]
        matched_cog_set = set()
        matched_kegg_set_d = set()
        for each_gene in gene_member:
            assigned_cog  = gene_to_cog_dict.get(each_gene, None)
            assigned_kegg = gene_to_kegg_dict.get(each_gene, None)
            if assigned_cog is not None:
                matched_cog_set.add(assigned_cog)
            if assigned_kegg is not None:
                matched_kegg_set_d.add(assigned_kegg)
        if len(matched_cog_set) > 0:
            og_to_cog_dict[each_og] = matched_cog_set
        if len(matched_kegg_set_d) > 0:
            og_to_kegg_dict[each_og] = matched_kegg_set_d

    # filter og_to_cog_dict
    og_to_cog_dict_filtered = dict()
    for each_og in og_to_cog_dict:
        matched_cog_set = og_to_cog_dict[each_og]
        if len(matched_cog_set) == 1:
            matched_cog = list(matched_cog_set)[0]
            matched_cog_desc = cog_desc_dict[matched_cog]
            if matched_cog_desc not in ['Uncharacterized protein']:
                og_to_cog_dict_filtered[each_og] = matched_cog_set
        else:
            current_desc_to_cog_dict = dict()
            for each_matched_cog in matched_cog_set:
                each_matched_cog_desc = cog_desc_dict[each_matched_cog]
                if each_matched_cog_desc not in ['Uncharacterized protein']:
                    if each_matched_cog_desc not in current_desc_to_cog_dict:
                        current_desc_to_cog_dict[each_matched_cog_desc] = set()
                    current_desc_to_cog_dict[each_matched_cog_desc].add(each_matched_cog)

            matched_cog_set_filtered = set()
            for each_desc in current_desc_to_cog_dict:
                same_desc_cog_list = sorted(list(current_desc_to_cog_dict[each_desc]))
                matched_cog_set_filtered.add(same_desc_cog_list[0])

            if len(matched_cog_set_filtered) > 0:
                og_to_cog_dict_filtered[each_og] = matched_cog_set_filtered

    # filter og_to_kegg_dict
    og_to_kegg_dict_filtered_b = dict()
    og_to_kegg_dict_filtered_c = dict()
    og_to_kegg_dict_filtered_d = dict()
    for each_og in og_to_kegg_dict:
        matched_kegg_set_d = og_to_kegg_dict[each_og]
        matched_kegg_set_b_filtered = set()
        matched_kegg_set_c_filtered = set()
        matched_kegg_set_d_filtered = set()
        for each_matched_ko_d in matched_kegg_set_d:
            each_matched_ko_b = ko_d_to_b_dict[each_matched_ko_d]
            each_matched_ko_c = ko_d_to_c_dict[each_matched_ko_d]
            each_matched_ko_d_desc = ko_desc_dict[each_matched_ko_d]
            if 'uncharacterized protein' not in each_matched_ko_d_desc:
                matched_kegg_set_b_filtered.add(each_matched_ko_b)
                matched_kegg_set_c_filtered.add(each_matched_ko_c)
                matched_kegg_set_d_filtered.add(each_matched_ko_d)
        if len(matched_kegg_set_d_filtered) > 0:
            og_to_kegg_dict_filtered_b[each_og] = matched_kegg_set_b_filtered
            og_to_kegg_dict_filtered_c[each_og] = matched_kegg_set_c_filtered
            og_to_kegg_dict_filtered_d[each_og] = matched_kegg_set_d_filtered

    op_txt_cog  = '%s/OG_COG.txt'  % op_dir
    op_txt_ko_b = '%s/OG_ko_B.txt' % op_dir
    op_txt_ko_c = '%s/OG_ko_C.txt' % op_dir
    op_txt_ko_d = '%s/OG_ko_D.txt' % op_dir

    op_txt_cog_handle = open(op_txt_cog, 'w')
    for i in sorted(list(og_to_cog_dict_filtered.keys())):
        matched_cog_set = og_to_cog_dict_filtered[i]
        for each_matched_cog in matched_cog_set:
            cog_desc = cog_desc_dict[each_matched_cog]
            op_txt_cog_handle.write('%s\t%s\t%s\n' % (i, each_matched_cog, cog_desc))
    op_txt_cog_handle.close()

    op_txt_ko_d_handle = open(op_txt_ko_d, 'w')
    for i in sorted(list(og_to_kegg_dict_filtered_d.keys())):
        matched_ko_d_set = og_to_kegg_dict_filtered_d[i]
        for each_ko_d in matched_ko_d_set:
            ko_d_desc = ko_desc_dict[each_ko_d]
            op_txt_ko_d_handle.write('%s\t%s\t%s\n' % (i, each_ko_d, ko_d_desc))
    op_txt_ko_d_handle.close()

    op_txt_ko_c_handle = open(op_txt_ko_c, 'w')
    for i in sorted(list(og_to_kegg_dict_filtered_c.keys())):
        matched_ko_c_set = og_to_kegg_dict_filtered_c[i]
        for each_ko_c in matched_ko_c_set:
            ko_c_desc = ko_desc_dict[each_ko_c]
            op_txt_ko_c_handle.write('%s\t%s\t%s\n' % (i, each_ko_c, ko_c_desc))
    op_txt_ko_c_handle.close()

    op_txt_ko_b_handle = open(op_txt_ko_b, 'w')
    for i in sorted(list(og_to_kegg_dict_filtered_b.keys())):
        matched_ko_b_set = og_to_kegg_dict_filtered_b[i]
        for each_ko_b in matched_ko_b_set:
            ko_b_desc = ko_desc_dict[each_ko_b]
            op_txt_ko_b_handle.write('%s\t%s\t%s\n' % (i, each_ko_b, ko_b_desc))
    op_txt_ko_b_handle.close()

    return og_to_cog_dict_filtered, cog_cate_dict, cog_desc_dict, og_to_kegg_dict_filtered_b, og_to_kegg_dict_filtered_c, og_to_kegg_dict_filtered_d, ko_desc_dict


def plot_gain_loss_plot_cog(og_set_gain, og_set_loss, og_cog_dict, cog_cate_dict, cog_cate_desc_dict, title_text, op_plot_cog):

    cog_cate_list = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'X']

    og_gain_cog_cate_stats_dict = dict()
    for each_gain in og_set_gain:
        gain_cog_set = og_cog_dict.get(each_gain, set())
        for each_gain_cog in gain_cog_set:
            each_gain_cate_str = cog_cate_dict[each_gain_cog]
            for each_gain_cate in list(each_gain_cate_str):
                if each_gain_cate not in og_gain_cog_cate_stats_dict:
                    og_gain_cog_cate_stats_dict[each_gain_cate] = 0
                og_gain_cog_cate_stats_dict[each_gain_cate] += 1

    og_loss_cog_cate_stats_dict = dict()
    for each_loss in og_set_loss:
        loss_cog_set = og_cog_dict.get(each_loss, set())
        for each_loss_cog in loss_cog_set:
            each_loss_cate_str = cog_cate_dict[each_loss_cog]
            for each_loss_cate in list(each_loss_cate_str):
                if each_loss_cate not in og_loss_cog_cate_stats_dict:
                    og_loss_cog_cate_stats_dict[each_loss_cate] = 0
                og_loss_cog_cate_stats_dict[each_loss_cate] += 1

    label_list = []
    num_list_gain = []
    num_list_loss = []
    for l in cog_cate_list:
        num_list_gain.append(og_gain_cog_cate_stats_dict.get(l, 0))
        num_list_loss.append(og_loss_cog_cate_stats_dict.get(l, 0))
        long_label = '%s - %s' % (l, cog_cate_desc_dict[l])
        long_label = long_label[:50]
        label_list.append(long_label)

    get_bar_plot(num_list_gain, 'Gain', num_list_loss, 'Loss', label_list, title_text, op_plot_cog)


def plot_gain_loss_plot_kegg(og_set_gain, og_set_loss, og_ko_b_dict, og_ko_c_dict, og_ko_d_dict, ko_desc_dict, gnm_id, op_plot_ko_b, op_plot_ko_c):

    ko_b_set_all = set()
    ko_c_set_all = set()
    og_gain_stats_dict_b = dict()
    og_gain_stats_dict_c = dict()
    for each_gain in og_set_gain:
        gain_ko_b_set = og_ko_b_dict.get(each_gain, set())
        gain_ko_c_set = og_ko_c_dict.get(each_gain, set())
        for each_ko_b in gain_ko_b_set:
            ko_b_set_all.add(each_ko_b)
            if each_ko_b not in og_gain_stats_dict_b:
                og_gain_stats_dict_b[each_ko_b] = 0
            og_gain_stats_dict_b[each_ko_b] += 1
        for each_ko_c in gain_ko_c_set:
            ko_c_set_all.add(each_ko_c)
            if each_ko_c not in og_gain_stats_dict_c:
                og_gain_stats_dict_c[each_ko_c] = 0
            og_gain_stats_dict_c[each_ko_c] += 1

    og_loss_stats_dict_b = dict()
    og_loss_stats_dict_c = dict()
    for each_loss in og_set_loss:
        loss_ko_b_set = og_ko_b_dict.get(each_loss, set())
        loss_ko_c_set = og_ko_c_dict.get(each_loss, set())
        for each_ko_b in loss_ko_b_set:
            ko_b_set_all.add(each_ko_b)
            if each_ko_b not in og_loss_stats_dict_b:
                og_loss_stats_dict_b[each_ko_b] = 0
            og_loss_stats_dict_b[each_ko_b] += 1
        for each_ko_c in loss_ko_c_set:
            ko_c_set_all.add(each_ko_c)
            if each_ko_c not in og_loss_stats_dict_c:
                og_loss_stats_dict_c[each_ko_c] = 0
            og_loss_stats_dict_c[each_ko_c] += 1

    ko_b_list_all_sorted = sorted(list(ko_b_set_all))
    ko_c_list_all_sorted = sorted(list(ko_c_set_all))

    ko_b_label_list = []
    ko_b_num_list_gain = []
    ko_b_num_list_loss = []
    for b in ko_b_list_all_sorted:
        ko_b_num_list_gain.append(og_gain_stats_dict_b.get(b, 0))
        ko_b_num_list_loss.append(og_loss_stats_dict_b.get(b, 0))
        long_label = '%s - %s' % (b, ko_desc_dict[b])
        long_label = long_label[:50]
        ko_b_label_list.append(long_label)
    get_bar_plot(ko_b_num_list_gain, 'Gain', ko_b_num_list_loss, 'Loss', ko_b_label_list, gnm_id, op_plot_ko_b)

    ko_c_label_list = []
    ko_c_num_list_gain = []
    ko_c_num_list_loss = []
    for c in ko_c_list_all_sorted:
        ko_c_num_list_gain.append(og_gain_stats_dict_c.get(c, 0))
        ko_c_num_list_loss.append(og_loss_stats_dict_c.get(c, 0))
        long_label = '%s - %s' % (c, ko_desc_dict[c])
        long_label = long_label[:50]
        ko_c_label_list.append(long_label)
    get_bar_plot(ko_c_num_list_gain, 'Gain', ko_c_num_list_loss, 'Loss', ko_c_label_list, gnm_id, op_plot_ko_c)


def ALE9(args):

    ale1_op_dir             = args['1']
    ale8_op_dir             = args['8']
    force_create_op_dir     = args['f']
    op_dir                  = args['o']
    cog_annotation_dir      = args['cog']
    kegg_annotation_dir     = args['kegg']
    cog_db_dir              = args['cog_db']

    gene_gain_loss_txt_dir  = '%s/Gene_gain_loss'   % ale8_op_dir
    fun_20_tab              = '%s/fun-20.tab'       % cog_db_dir

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    cog_cate_desc_dict = dict()
    for each_line in open(fun_20_tab):
        each_line_split = each_line.strip().split('\t')
        cog_cate_desc_dict[each_line_split[0]] = each_line_split[2]

    og_cog_dict, cog_cate_dict, cog_desc_dict, og_ko_b_dict, og_ko_c_dict, og_ko_d_dict, ko_desc_dict = get_og_fun(ale1_op_dir, 'faa', cog_annotation_dir, kegg_annotation_dir, op_dir)

    gene_gain_txt_re   = '%s/*_gene_gain.txt' % gene_gain_loss_txt_dir
    gene_gain_txt_list = glob.glob(gene_gain_txt_re)

    for gene_gain_txt in gene_gain_txt_list:
        _, _, f_base, _ = sep_path_basename_ext(gene_gain_txt)
        gnm_id          = f_base.split('_gene_gain')[0]
        gene_loss_txt   = gene_gain_txt.replace('_gene_gain.txt', '_gene_loss.txt')
        op_plot_cog     = '%s/%s_COG.pdf'    % (op_dir, gnm_id)
        op_plot_ko_b    = '%s/%s_KEGG_B.pdf' % (op_dir, gnm_id)
        op_plot_ko_c    = '%s/%s_KEGG_C.pdf' % (op_dir, gnm_id)

        og_set_gain = set()
        for each_line in open(gene_gain_txt):
            og_set_gain.add(each_line.strip())

        og_set_loss = set()
        for each_line in open(gene_loss_txt):
            og_set_loss.add(each_line.strip())

        plot_gain_loss_plot_cog(og_set_gain, og_set_loss, og_cog_dict, cog_cate_dict, cog_cate_desc_dict, gnm_id, op_plot_cog)
        plot_gain_loss_plot_kegg(og_set_gain, og_set_loss, og_ko_b_dict, og_ko_c_dict, og_ko_d_dict, ko_desc_dict, gnm_id, op_plot_ko_b, op_plot_ko_c)

        _, _, gain_txt_base, _ = sep_path_basename_ext(gene_gain_txt)
        _, _, loss_txt_base, _ = sep_path_basename_ext(gene_loss_txt)

        op_gene_gain_txt_cog   = '%s/%s_COG.txt'  % (op_dir, gain_txt_base)
        op_gene_loss_txt_cog   = '%s/%s_COG.txt'  % (op_dir, loss_txt_base)
        op_gene_gain_txt_kegg  = '%s/%s_KEGG.txt' % (op_dir, gain_txt_base)
        op_gene_loss_txt_kegg  = '%s/%s_KEGG.txt' % (op_dir, loss_txt_base)

        op_gene_gain_txt_cog_handle = open(op_gene_gain_txt_cog, 'w')
        op_gene_gain_txt_kegg_handle = open(op_gene_gain_txt_kegg, 'w')
        for each_gain in sorted(list(og_set_gain)):
            # write out cog
            gain_cog_set = og_cog_dict.get(each_gain, set())
            for gain_cog in gain_cog_set:
                gain_cog_desc = cog_desc_dict.get(gain_cog, '')
                op_gene_gain_txt_cog_handle.write('%s\t%s\t%s\n' % (each_gain, gain_cog, gain_cog_desc))
            # write out kegg
            gain_ko_d_set = og_ko_d_dict.get(each_gain, set())
            for gain_ko_d in gain_ko_d_set:
                gain_ko_d_desc = ko_desc_dict.get(gain_ko_d, '')
                op_gene_gain_txt_kegg_handle.write('%s\t%s\t%s\n' % (each_gain, gain_ko_d, gain_ko_d_desc))
        op_gene_gain_txt_cog_handle.close()
        op_gene_gain_txt_kegg_handle.close()

        op_gene_loss_txt_cog_handle = open(op_gene_loss_txt_cog, 'w')
        op_gene_loss_txt_kegg_handle = open(op_gene_loss_txt_kegg, 'w')
        for each_loss in sorted(list(og_set_loss)):
            # write out cog
            loss_cog_set = og_cog_dict.get(each_loss, set())
            for loss_cog in loss_cog_set:
                loss_cog_desc = cog_desc_dict.get(loss_cog, '')
                op_gene_loss_txt_cog_handle.write('%s\t%s\t%s\n' % (each_loss, loss_cog, loss_cog_desc))
            # write out kegg
            loss_ko_d_set = og_ko_d_dict.get(each_loss, set())
            for loss_ko_d in loss_ko_d_set:
                loss_ko_d_desc = ko_desc_dict.get(loss_ko_d, '')
                op_gene_loss_txt_kegg_handle.write('%s\t%s\t%s\n' % (each_loss, loss_ko_d, loss_ko_d_desc))
        op_gene_loss_txt_cog_handle.close()
        op_gene_loss_txt_kegg_handle.close()


if __name__ == '__main__':

    ALE9_parser = argparse.ArgumentParser()
    ALE9_parser.add_argument('-1',      required=True,                          help='ALE1 output directory')
    ALE9_parser.add_argument('-8',      required=True,                          help='ALE8 output directory')
    ALE9_parser.add_argument('-cog',    required=False, default=None,           help='COG annotation dir')
    ALE9_parser.add_argument('-kegg',   required=False, default=None,           help='KEGG annotation dir')
    ALE9_parser.add_argument('-cog_db', required=False, default=None,           help='COG database dir')
    ALE9_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    ALE9_parser.add_argument('-o',      required=True,                          help='output directory')
    args = vars(ALE9_parser.parse_args())
    ALE9(args)
