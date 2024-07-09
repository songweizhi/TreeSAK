import os


def group_marker(taxa_counts_tats_op_txt, marker_seq_dir, marker_rank_cutoff_str, op_dir):

    marker_rank_cutoff_list = marker_rank_cutoff_str.split(',')

    marker_score_dict = dict()
    header_index_dict = dict()
    for each_marker in open(taxa_counts_tats_op_txt):
        each_marker_split = each_marker.replace('\n', '').split('\t')
        if each_marker.startswith('MarkerID\t'):
            header_index_dict = {k: v for v, k in enumerate(each_marker_split)}
        else:
            marker_id    = each_marker_split[header_index_dict['MarkerID']]
            marker_score = int(each_marker_split[header_index_dict['RankA_B']])
            marker_score_dict[marker_id] = marker_score

    marker_list_sorted_best_to_wrost = [i[0] for i in sorted(marker_score_dict.items(), key=lambda x: x[1])]
    marker_list_sorted_wrost_to_best = [i[0] for i in sorted(marker_score_dict.items(), key=lambda x: x[1])][::-1]

    for each_cutoff in marker_rank_cutoff_list:

        marker_num_rounded = round(len(marker_list_sorted_wrost_to_best)*float(each_cutoff)/100)
        marker_list_best   = marker_list_sorted_best_to_wrost[:marker_num_rounded]
        marker_list_worst  = marker_list_sorted_wrost_to_best[:marker_num_rounded]
        seq_dir_best       = '%s/best%s'  % (op_dir, each_cutoff)
        seq_dir_worst      = '%s/worst%s' % (op_dir, each_cutoff)

        os.system('mkdir %s' % seq_dir_best)
        os.system('mkdir %s' % seq_dir_worst)

        # get the best markers
        for bm in marker_list_best:
            os.system('cp %s/%s.fa %s/' % (marker_seq_dir, bm, seq_dir_best))

        # get the worst markers
        for wm in marker_list_worst:
            os.system('cp %s/%s.fa %s/' % (marker_seq_dir, wm, seq_dir_worst))


# TaxaCountStats_output_txt = '/Users/songweizhi/Desktop/TaxaCountStats_output.txt'
# qualified_og_seq_dir      = '/Users/songweizhi/Desktop/qualified_OGs'
# step_2_op_dir             = '/Users/songweizhi/Desktop/SplitScore2_op_dir'
# marker_rank_cutoff_str    = '25,50,75'

TaxaCountStats_output_txt = '/scratch/PI/ocessongwz/Sponge_r220/4_OMA_wd_115_genomes/OMA_wd/Output/OMA_75_5_dRep85_115_cov90/SplitScore2_op_dir/get_taxa_count_stats_wd/TaxaCountStats_output.txt'
qualified_og_seq_dir      = '/scratch/PI/ocessongwz/Sponge_r220/4_OMA_wd_115_genomes/OMA_wd/Output/OMA_75_5_dRep85_115_cov90/SplitScore1_op_dir/qualified_OGs'
step_2_op_dir             = '/scratch/PI/ocessongwz/Sponge_r220/4_OMA_wd_115_genomes/OMA_wd/Output/OMA_75_5_dRep85_115_cov90/SplitScore2_op_dir'
marker_rank_cutoff_str    = '10'

group_marker(TaxaCountStats_output_txt, qualified_og_seq_dir, marker_rank_cutoff_str, step_2_op_dir)

