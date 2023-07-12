import os


def read_in_assessment_pa_2_txt(assessment_pa_2_txt):

    pa_pct_list = []
    pa_pct_to_marker_dict = dict()
    for each_marker in open(assessment_pa_2_txt):
        each_marker_split = each_marker.strip().split('\t')
        if each_marker.startswith('Marker\t'):
            pa_pct_list = [int(i) for i in each_marker_split[1:]]

            # initialize pa_pct_to_marker_dict
            for pa_pct in pa_pct_list:
                pa_pct_to_marker_dict[pa_pct] = set()

        elif not each_marker.startswith('Total\t'):
            marker_id = each_marker_split[0]
            for (pct, pa) in zip(pa_pct_list, each_marker_split[1:]):
                if pa == '1':
                    pa_pct_to_marker_dict[pct].add(marker_id)

    return pa_pct_to_marker_dict


def get_marker_set_dict(assessment_pa_2_txt, assessment_deltall_txt, deltall_keep_pct_list, min_marker_num):
    # read in assessment_pa_2_txt
    pa_pct_to_marker_dict = read_in_assessment_pa_2_txt(assessment_pa_2_txt)
    # for each in pa_pct_to_marker_dict:
    #     print('%s(%s)\t%s' % (each, len(pa_pct_to_marker_dict[each]), pa_pct_to_marker_dict[each]))

    # store marker in list according to their DeltaLl scores
    marker_list_by_score = []
    for each_marker in open(assessment_deltall_txt):
        if not each_marker.startswith('Marker\t'):
            each_marker_split = each_marker.strip().split('\t')
            marker_id = each_marker_split[0]
            marker_list_by_score.append(marker_id)

    # get intersections
    marker_set_dict = dict()
    for each_keep_pct in deltall_keep_pct_list:
        keep_num = round(len(marker_list_by_score) * each_keep_pct / 100)
        if keep_num >= min_marker_num:
            marker_to_keep = marker_list_by_score[:keep_num]
            # print('deltall_keep_pct(top %s%s)(%s)\t%s' % (each_keep_pct, '%', keep_num, marker_to_keep))
            for each_pa_pct in pa_pct_to_marker_dict:
                current_pa_pct_marker_set = pa_pct_to_marker_dict[each_pa_pct]
                marker_set_key = 'deltall%s_pa%s' % (each_keep_pct, each_pa_pct)
                marker_shared = set(marker_to_keep).intersection(current_pa_pct_marker_set)
                marker_set_dict[marker_set_key] = marker_shared

    return marker_set_dict


# input
# wd                     = '/Users/songweizhi/Desktop/DateArTree/01_marker_gene_selection'
# trimmed_aln_dir        = '%s/best_hit_by_marker_5_aln_trimmed'                             % wd
# assessment_pa_2_txt    = '%s/op_marker_assessment_2.txt'                                   % wd
# assessment_deltall_txt = '%s/deltall_op_reformatted.txt'                                   % wd
# deltall_keep_pct_list  = [25, 50, 75, 100]
# min_marker_num         = 20
# force_create_dir       = True
# catfasta2phyml_pl      = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/catfasta2phyml.pl'

wd                     = '/Users/songweizhi/Desktop/DateArTree/01_marker_gene_selection'
trimmed_aln_dir        = '%s/best_hit_by_marker_5_aln_trimmed'                             % wd
assess_summary_2_txt_by_marker    = '%s/op_marker_assessment_2.txt'                                   % wd
assessment_deltall_txt = '%s/deltall_op_reformatted.txt'                                   % wd
deltall_keep_pct_list  = [25, 50, 75, 100]
min_marker_num         = 20
force_create_dir       = True
catfasta2phyml_pl      = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/catfasta2phyml.pl'


marker_set_dict = get_marker_set_dict(assess_summary_2_txt_by_marker, assessment_deltall_txt, deltall_keep_pct_list, min_marker_num)
for each_marker_set in marker_set_dict:
    current_marker_set = marker_set_dict[each_marker_set]
    if len(current_marker_set) >= min_marker_num:

        #print('%s\t%s\t%s' % (each_marker_set, len(current_marker_set), current_marker_set))
        pwd_marker_id_txt = '%s/%s.txt'         % (wd, each_marker_set)
        pwd_aln_dir       = '%s/%s_aln_trimmed' % (wd, each_marker_set)

        # write out marker id
        with open(pwd_marker_id_txt, 'w') as pwd_marker_id_txt_handle:
            pwd_marker_id_txt_handle.write('\n'.join(current_marker_set))

        # copy marker alignments into corresponding dir
        if force_create_dir is True:
            if os.path.isdir(pwd_aln_dir) is True:
                os.system('rm -r %s' % pwd_aln_dir)
        os.system('mkdir %s' % pwd_aln_dir)
        for each_aln in current_marker_set:
            pwd_aln = '%s/%s.aln' % (trimmed_aln_dir, each_aln)
            cp_cmd = 'cp %s %s/' % (pwd_aln, pwd_aln_dir)
            os.system(cp_cmd)

        # concatenate alignment
        pwd_concatenate_aln             = '%s/%s_aln_trimmed_concatenated.phy'                      % (wd, each_marker_set)
        pwd_concatenate_aln_partitions  = '%s/%s_aln_trimmed_concatenated_partitions.txt'           % (wd, each_marker_set)
        catfasta2phyml_cmd              = 'perl %s --sequential --concatenate %s/*.aln > %s 2> %s'  % (catfasta2phyml_pl, pwd_aln_dir, pwd_concatenate_aln, pwd_concatenate_aln_partitions)
        print(catfasta2phyml_cmd)
        os.system(catfasta2phyml_cmd)










