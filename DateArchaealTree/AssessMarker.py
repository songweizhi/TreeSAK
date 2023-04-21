import os
import glob
from Bio import SeqIO


def parse_deltall_stdout(deltall_stdout_txt, summary_txt):

    deltall_op_dict = dict()
    for each_line in open(deltall_stdout_txt):
        if not ((each_line.startswith('WARNING:')) or (each_line.startswith('awk:'))):
            each_line_split = each_line.strip().split('\t')
            marker_id = each_line_split[0]
            value = float(each_line_split[1])
            if marker_id not in deltall_op_dict:
                deltall_op_dict[marker_id] = [value]
            else:
                deltall_op_dict[marker_id].append(value)

    metric_1_dict = dict()
    metric_2_dict = dict()
    for each_marker in deltall_op_dict:
        metric_1_value = float("{0:.2f}".format(deltall_op_dict[each_marker][0]))
        metric_2_value = float("{0:.2f}".format(deltall_op_dict[each_marker][1]))
        metric_1_dict[each_marker] = metric_1_value
        metric_2_dict[each_marker] = metric_2_value

    metric_1_dict_sorted = {k: v for k, v in sorted(metric_1_dict.items(), key=lambda item: item[1])[::-1]}
    metric_2_dict_sorted = {k: v for k, v in sorted(metric_2_dict.items(), key=lambda item: item[1])}

    metric_1_score_dict = dict()
    metric_1_score = 1
    for each_marker_1 in metric_1_dict_sorted:
        metric_1_score_dict[each_marker_1] = metric_1_score
        metric_1_score += 1

    metric_2_score_dict = dict()
    metric_2_score = 1
    for each_marker_2 in metric_2_dict_sorted:
        metric_2_score_dict[each_marker_2] = metric_2_score
        metric_2_score += 1

    overall_score_dict = dict()
    for each_marker in deltall_op_dict:
        metric_score_1 = metric_1_score_dict[each_marker]
        metric_score_2 = metric_2_score_dict[each_marker]
        metric_score_overall = metric_score_1 + metric_score_2
        overall_score_dict[each_marker] = metric_score_overall

    overall_score_dict_sorted = {k: v for k, v in sorted(overall_score_dict.items(), key=lambda item: item[1])}

    summary_txt_handle = open(summary_txt, 'w')
    summary_txt_handle.write('Marker\tmetric1\tmetric1_score\tmetric2\tmetric2_score\toverall_score\n')
    for each_marker in overall_score_dict_sorted:
        metric_value_1 = metric_1_dict[each_marker]
        metric_value_2 = metric_2_dict[each_marker]
        metric_score_1 = metric_1_score_dict[each_marker]
        metric_score_2 = metric_2_score_dict[each_marker]
        metric_score_overall = overall_score_dict_sorted[each_marker]
        summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_marker, metric_value_1, metric_score_1, metric_value_2, metric_score_2, metric_score_overall))
    summary_txt_handle.close()


def assess_markers_pa(trimmed_aln_dir, gnm_meta_txt, present_pct_cutoff_list, assess_summary_1_txt, assess_summary_2_txt):

    trimmed_aln_file_re = '%s/*.aln' % trimmed_aln_dir
    trimmed_aln_file_list = [os.path.basename(file_name) for file_name in glob.glob(trimmed_aln_file_re)]

    # read in genome metadata
    domain_to_gnm_dict = dict()
    gnm_to_domain_dict = dict()
    for each_gnm in open(gnm_meta_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        domain_name = each_gnm_split[1]
        gnm_to_domain_dict[gnm_id] = domain_name

        if domain_name not in domain_to_gnm_dict:
            domain_to_gnm_dict[domain_name] = {gnm_id}
        else:
            domain_to_gnm_dict[domain_name].add(gnm_id)

    assess_summary_1_txt_handle = open(assess_summary_1_txt, 'w')
    assess_summary_2_txt_handle = open(assess_summary_2_txt, 'w')
    gnm_to_marker_dict = dict()
    marker_to_gnm_dict = dict()
    cutoff_to_qualified_marker_dict = dict()
    assess_summary_1_txt_handle.write('Marker\tArchaea\tEukaryota\n')
    assess_summary_2_txt_handle.write('Marker\t%s\n' % '\t'.join([str(i) for i in present_pct_cutoff_list]))
    for each_aln in trimmed_aln_file_list:
        marker_id = each_aln.split('.aln')[0]
        pwd_aln = '%s/%s' % (trimmed_aln_dir, each_aln)
        gnm_set = set()
        for each_seq in SeqIO.parse(pwd_aln, 'fasta'):
            gnm_id = each_seq.id
            gnm_set.add(gnm_id)
            if gnm_id not in gnm_to_marker_dict:
                gnm_to_marker_dict[gnm_id] = {marker_id}
            else:
                gnm_to_marker_dict[gnm_id].add(marker_id)
        marker_to_gnm_dict[marker_id] = gnm_set

        gnm_num_ar = 0
        gnm_num_eu = 0
        for each_g in gnm_set:
            g_domain = gnm_to_domain_dict[each_g]
            if g_domain == 'Archaea':
                gnm_num_ar += 1
            if g_domain == 'Eukaryota':
                gnm_num_eu += 1
        gnm_pct_ar = float("{0:.2f}".format(gnm_num_ar / 133 * 100))
        gnm_pct_eu = float("{0:.2f}".format(gnm_num_eu / 27 * 100))

        # assessment
        assessment_result_list = []
        for present_pct_cutoff in present_pct_cutoff_list:
            if (gnm_pct_ar >= present_pct_cutoff) and (gnm_pct_eu >= present_pct_cutoff):
                assessment_result_list.append('1')
                if str(present_pct_cutoff) not in cutoff_to_qualified_marker_dict:
                    cutoff_to_qualified_marker_dict[str(present_pct_cutoff)] = [marker_id]
                else:
                    cutoff_to_qualified_marker_dict[str(present_pct_cutoff)].append(marker_id)
            else:
                assessment_result_list.append('0')
        assess_summary_1_txt_handle.write('%s\t%s\t%s\n' % (marker_id, gnm_pct_ar, gnm_pct_eu))
        assess_summary_2_txt_handle.write('%s\t%s\n' % (marker_id, '\t'.join(assessment_result_list)))

    summary_list = [len(cutoff_to_qualified_marker_dict.get(str(i), [])) for i in present_pct_cutoff_list]
    summary_list_str = [str(j) for j in summary_list]
    assess_summary_2_txt_handle.write('Total\t%s\n' % ('\t'.join(summary_list_str)))
    assess_summary_1_txt_handle.close()
    assess_summary_2_txt_handle.close()

    return cutoff_to_qualified_marker_dict, gnm_to_marker_dict, marker_to_gnm_dict


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


def AssessMarker(assess_marker_wd, trimmed_aln_dir, gnm_group_txt, deltall_stdout_txt, present_pct_cutoff_list, deltall_keep_pct_list, min_marker_pct_per_gnm, min_marker_num, force_create_dir, catfasta2phyml_pl):

    # define output file name
    assess_summary_deltall_txt     = '%s/assessment_deltall.txt'      % assess_marker_wd
    assess_summary_1_txt_by_marker = '%s/assessment_pa_1.txt'         % assess_marker_wd
    assess_summary_2_txt_by_marker = '%s/assessment_pa_2.txt'         % assess_marker_wd
    assess_summary_txt_by_genome   = '%s/assessment_pa_by_genome.txt' % assess_marker_wd

    # parse deltall stdout
    parse_deltall_stdout(deltall_stdout_txt, assess_summary_deltall_txt)

    # assess markers
    cutoff_to_qualified_marker_dict, gnm_to_marker_dict, marker_to_gnm_dict = assess_markers_pa(trimmed_aln_dir, gnm_group_txt, present_pct_cutoff_list, assess_summary_1_txt_by_marker, assess_summary_2_txt_by_marker)

    # write out qualified markers
    for each_cutoff in cutoff_to_qualified_marker_dict:
        qualified_m_list = sorted(cutoff_to_qualified_marker_dict[each_cutoff])
        pwd_op_txt = '%s/assessment_pa_qualified_marker_%s.txt' % (assess_marker_wd, each_cutoff)
        with open(pwd_op_txt, 'w') as pwd_op_txt_handle:
            pwd_op_txt_handle.write('\n'.join(qualified_m_list))

    # write out summary by genomes
    assess_summary_txt_by_genome_handle = open(assess_summary_txt_by_genome, 'w')
    assess_summary_txt_by_genome_handle.write('Cutoff\tGenome\tMarker_all\tMarker_qualified\tMarker_qualified_pct(cutoff:%s)\n' % min_marker_pct_per_gnm)
    for each_cutoff in present_pct_cutoff_list:
        current_cutoff_qualified_marker_list = cutoff_to_qualified_marker_dict.get(str(each_cutoff), [])
        if len(current_cutoff_qualified_marker_list) > 0:
            qualified_gnm_set = set()
            for each_gnm in gnm_to_marker_dict:
                gnm_identified_marker_set = gnm_to_marker_dict[each_gnm]
                gnm_identified_marker_set_qualified = set()
                for identified_marker in gnm_identified_marker_set:
                    if identified_marker in current_cutoff_qualified_marker_list:
                        gnm_identified_marker_set_qualified.add(identified_marker)
                gnm_identified_marker_set_qualified_pct = len(gnm_identified_marker_set_qualified)*100/len(current_cutoff_qualified_marker_list)
                gnm_identified_marker_set_qualified_pct = float("{0:.2f}".format(gnm_identified_marker_set_qualified_pct))
                if gnm_identified_marker_set_qualified_pct >= min_marker_pct_per_gnm:
                    qualified_gnm_set.add(each_gnm)
                else:
                    assess_summary_txt_by_genome_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_cutoff, each_gnm, len(gnm_identified_marker_set), len(gnm_identified_marker_set_qualified), gnm_identified_marker_set_qualified_pct))
        assess_summary_txt_by_genome_handle.write('\n')
    assess_summary_txt_by_genome_handle.close()

    # select marker gene and concatenate the alignments
    marker_set_dict = get_marker_set_dict(assess_summary_2_txt_by_marker, assess_summary_deltall_txt, deltall_keep_pct_list, min_marker_num)
    for each_marker_set in marker_set_dict:
        current_marker_set = marker_set_dict[each_marker_set]
        if len(current_marker_set) >= min_marker_num:

            #print('%s\t%s\t%s' % (each_marker_set, len(current_marker_set), current_marker_set))
            pwd_iqtree_dir    = '%s/%s_iqtree_wd'       % (assess_marker_wd, each_marker_set)
            pwd_marker_id_txt = '%s/%s_marker_id.txt'   % (pwd_iqtree_dir, each_marker_set)
            pwd_aln_dir       = '%s/%s_aln_trimmed'     % (pwd_iqtree_dir, each_marker_set)

            # copy marker alignments into corresponding dir
            if force_create_dir is True:
                if os.path.isdir(pwd_iqtree_dir) is True:
                    os.system('rm -r %s' % pwd_iqtree_dir)
            os.system('mkdir %s' % pwd_iqtree_dir)
            os.system('mkdir %s' % pwd_aln_dir)

            # write out marker id
            with open(pwd_marker_id_txt, 'w') as pwd_marker_id_txt_handle:
                pwd_marker_id_txt_handle.write('\n'.join(current_marker_set))

            for each_aln in current_marker_set:
                pwd_aln = '%s/%s.aln' % (trimmed_aln_dir, each_aln)
                cp_cmd = 'cp %s %s/' % (pwd_aln, pwd_aln_dir)
                os.system(cp_cmd)

            # concatenate alignment
            pwd_concatenate_aln            = '%s/%s_concatenated.phy'                                 % (pwd_iqtree_dir, each_marker_set)
            pwd_concatenate_aln_partitions = '%s/%s_partitions.txt'                                   % (pwd_iqtree_dir, each_marker_set)
            catfasta2phyml_cmd             = 'perl %s --sequential --concatenate %s/*.aln > %s 2> %s' % (catfasta2phyml_pl, pwd_aln_dir, pwd_concatenate_aln, pwd_concatenate_aln_partitions)
            #print(catfasta2phyml_cmd)
            os.system(catfasta2phyml_cmd)

            # get guide tree
            get_guide_cmd = 'iqtree -s %s_concatenated.phy --prefix %s_guide_tree --seqtype AA -m LG -T 12 -B 1000 --alrt 1000' % (each_marker_set, each_marker_set)
            # print(get_guide_cmd)

            # run C60 + PMSF
            c60_pmsf_cmd = 'iqtree -s %s_concatenated.phy --prefix %s --seqtype AA -m LG+G+F+C60 -T 12 -B 1000 --alrt 1000 -ft %s_guide_tree.treefile' % (each_marker_set, each_marker_set, each_marker_set)
            # print(c60_pmsf_cmd)

            # generate job script
            pwd_js = '%s/js_%s_iqtree.sh' % (assess_marker_wd, each_marker_set)
            with open(pwd_js, 'w') as pwd_js_handle:
                #pwd_js_handle.write('#!/bin/bash\n#SBATCH --nodelist cl007\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task 12\n\n')
                pwd_js_handle.write('#!/bin/bash\n#SBATCH\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task 12\n\n')
                pwd_js_handle.write('cd %s_iqtree_wd\n' % each_marker_set)
                pwd_js_handle.write(get_guide_cmd + '\n')
                pwd_js_handle.write(c60_pmsf_cmd + '\n')


# inputs
assess_marker_wd        = '/home-user/wzsong/DateArTree/04_dating_Williams_2017_45_arCOG_assess_marker'
trimmed_aln_dir         = '/home-user/wzsong/DateArTree/02_identify_marker_gene_Williams_2017_45_arCOG/best_hit_by_marker_5_aln_trimmed'
gnm_group_txt           = '/home-user/wzsong/DateArTree/01_genome_selection/gnm_metadata.txt'
deltall_stdout_txt      = '/home-user/wzsong/DateArTree/02_identify_marker_gene_Williams_2017_45_arCOG_DeltaLL/nohup.out'
present_pct_cutoff_list = [25, 50, 75, 85, 100]
deltall_keep_pct_list   = [25, 50, 75, 100]
min_marker_pct_per_gnm  = 75
min_marker_num          = 20
force_create_dir        = True
#catfasta2phyml_pl       = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/catfasta2phyml.pl'
catfasta2phyml_pl       = '/home-user/wzsong/Scripts/catfasta2phyml.pl'

AssessMarker(assess_marker_wd, trimmed_aln_dir, gnm_group_txt, deltall_stdout_txt, present_pct_cutoff_list, deltall_keep_pct_list, min_marker_pct_per_gnm, min_marker_num, force_create_dir, catfasta2phyml_pl)


'''
Note
1. Extra genomes in gnm_metadata.txt won't affect assessment results.
2. Genomes can not be found in gnm_metadata.txt will trigger an error.
3. Alignments in {trimmed_aln_dir} need to be trimmed before assessment
4. Sequences in MSAs need to be named by genome id.
'''
