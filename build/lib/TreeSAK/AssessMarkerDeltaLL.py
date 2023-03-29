import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO


AssessMarkerDeltaLL_usage = '''
============================= AssessMarkerDeltaLL example commands =============================

Dependencies: iqtree

# example commands
TreeSAK AssessMarkerDeltaLL -deltall DeltaLL_stdout.txt -o op_dir -c 25-50-75-100 -mmn 20 -aln trimmed_aln_dir -jst 6 -qsub -f 

================================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, partition_file):

    concatenated_msa_fasta = '%s.fasta' % concatenated_msa_phy
    msa_file_re            = '%s/*.%s'  % (msa_dir, msa_ext)
    msa_file_list          = [os.path.basename(file_name) for file_name in glob.glob(msa_file_re)]
    msa_file_list_sorted   = sorted(msa_file_list)

    complete_gnm_set = set()
    for each_msa_file in msa_file_list:
        pwd_msa = '%s/%s' % (msa_dir, each_msa_file)
        for each_seq in SeqIO.parse(pwd_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)

    complete_gnm_list_sorted = sorted([i for i in complete_gnm_set])

    # initialize concatenated msa dict
    gnm_to_seq_dict = {i: '' for i in complete_gnm_list_sorted}
    msa_len_dict = dict()
    for each_msa_file in msa_file_list_sorted:
        gene_id = each_msa_file.split('.' + msa_ext)[0]

        # read in msa
        current_msa_len = 0
        current_msa_len_set = set()
        pwd_current_msa = '%s/%s' % (msa_dir, each_msa_file)
        current_msa_seq_dict = dict()
        for each_seq in SeqIO.parse(pwd_current_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)
            current_msa_seq_dict[each_seq.id] = str(each_seq.seq)
            current_msa_len_set.add(len(each_seq.seq))
            current_msa_len = len(each_seq.seq)

        if len(current_msa_len_set) != 1:
            print('Sequences with different length were found in %s, program exited!' % each_msa_file)
            exit()

        msa_len_dict[gene_id] = current_msa_len

        # add sequence to concatenated msa dict
        for each_gnm in complete_gnm_list_sorted:
            msa_seq = current_msa_seq_dict.get(each_gnm, current_msa_len*'-')
            gnm_to_seq_dict[each_gnm] += msa_seq

    # write out concatenated msa
    concatenated_msa_handle = open(concatenated_msa_fasta, 'w')
    for each_gnm in complete_gnm_list_sorted:
        concatenated_msa_handle.write('>%s\n' % each_gnm)
        concatenated_msa_handle.write('%s\n' % gnm_to_seq_dict[each_gnm])
    concatenated_msa_handle.close()

    # write out partition file
    end_pos = 0
    partition_file_handle = open(partition_file, 'w')
    for each_m in msa_file_list_sorted:
        gene_id = each_m.split('.' + msa_ext)[0]
        current_m_len = msa_len_dict[gene_id]
        partition_file_handle.write('%s = %s-%s\n' % (each_m, (end_pos + 1), (end_pos + current_m_len)))
        end_pos += current_m_len
    partition_file_handle.close()

    # convert msa in fasta to phy
    AlignIO.convert(concatenated_msa_fasta, 'fasta', concatenated_msa_phy, 'phylip-relaxed')


def submit_js(js):
    current_wd = os.getcwd()
    js_path, js_basename, js_ext = sep_path_basename_ext(js)
    os.chdir(js_path)
    os.system('qsub %s%s' % (js_basename, js_ext))
    os.chdir(current_wd)


def AssessMarkerDeltaLL(args):

    deltall_stdout_txt      = args['deltall']
    op_dir                  = args['o']
    deltall_keep_pct_str    = args['c']
    min_marker_num          = args['mmn']
    trimmed_aln_dir         = args['aln']
    js_cpu_num              = args['jst']
    submit_job              = args['qsub']
    force_overwrite         = args['f']

    deltall_keep_pct_list = [int(i) for i in deltall_keep_pct_str.split('-')]

    deltall_stdout_path, deltall_stdout_basename, deltall_stdout_ext = sep_path_basename_ext(deltall_stdout_txt)

    # create dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # define file name
    deltall_stdout_summary_txt = '%s/%s_summary.txt' % (op_dir, deltall_stdout_basename)

    # read in deltall_stdout_txt
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

    # assigned score to marker
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

    overall_score_dict_sorted     = {k: v for k, v in sorted(overall_score_dict.items(), key=lambda item: item[1])}
    marker_list_sorted_by_deltall = [k for k, v in sorted(overall_score_dict.items(), key=lambda item: item[1])]

    # write out summary_txt
    summary_txt_handle = open(deltall_stdout_summary_txt, 'w')
    summary_txt_handle.write('Marker\tmetric1\tmetric1_score\tmetric2\tmetric2_score\toverall_score\n')
    for each_marker in overall_score_dict_sorted:
        metric_value_1 = metric_1_dict[each_marker]
        metric_value_2 = metric_2_dict[each_marker]
        metric_score_1 = metric_1_score_dict[each_marker]
        metric_score_2 = metric_2_score_dict[each_marker]
        metric_score_overall = overall_score_dict_sorted[each_marker]
        summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_marker, metric_value_1, metric_score_1, metric_value_2, metric_score_2, metric_score_overall))
    summary_txt_handle.close()

    # get qualified marker list
    for each_keep_pct in deltall_keep_pct_list:
        marker_num_to_keep = round(len(marker_list_sorted_by_deltall)*each_keep_pct/100)
        markers_to_keep_id_list = marker_list_sorted_by_deltall[:marker_num_to_keep]

        if marker_num_to_keep < min_marker_num:
            print('Ignored DeltaLL cutoff at %s , the number of qualified markers (%s) less than %s' % (each_keep_pct, marker_num_to_keep, min_marker_num))
        else:
            pwd_aln_dir                     = '%s/%s_DeltaLL_%s_trimmed_aln'                % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_aln_concatenated            = '%s/%s_DeltaLL_%s_concatenated.phy'           % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_aln_concatenated_partitions = '%s/%s_DeltaLL_%s_concatenated_partition.txt' % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_iqtree_guide_tree_wd        = '%s/%s_DeltaLL_%s_iqtree_guide_tree'          % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_iqtree_c60_pmsf_wd          = '%s/%s_DeltaLL_%s_iqtree_C60_PMSF'            % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_iqtree_js                   = '%s/js_%s_DeltaLL_%s_iqtree.sh'               % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)

            # create dir
            os.mkdir(pwd_aln_dir)

            # copy msa of qualified markers
            for each_marker in markers_to_keep_id_list:
                pwd_marker_aln = '%s/%s.aln' % (trimmed_aln_dir, each_marker)
                cp_cmd = 'cp %s %s/' % (pwd_marker_aln, pwd_aln_dir)
                os.system(cp_cmd)

            # concatenate msa
            catfasta2phy(pwd_aln_dir, 'aln', pwd_aln_concatenated, pwd_aln_concatenated_partitions)

            # create dir
            os.mkdir(pwd_iqtree_guide_tree_wd)
            os.mkdir(pwd_iqtree_c60_pmsf_wd)

            # run iqtree
            get_guide_tree_cmd = 'iqtree -s ../%s_DeltaLL_%s_concatenated.phy --prefix guide_tree --seqtype AA -m LG -T %s -B 1000 --alrt 1000 --quiet'                                                                      % (deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct, js_cpu_num)
            get_c60_tree_cmd   = 'iqtree -s ../%s_DeltaLL_%s_concatenated.phy --prefix concatenated --seqtype AA -m LG+G+F+C60 -T %s -B 1000 --alrt 1000 --quiet -ft ../%s_DeltaLL_%s_iqtree_guide_tree/guide_tree.treefile' % (deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct, js_cpu_num, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)

            # write job script
            with open(pwd_iqtree_js, 'w') as pwd_iqtree_js_handle:
                pwd_iqtree_js_handle.write('#!/bin/bash\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task %s\n\n' % js_cpu_num)
                pwd_iqtree_js_handle.write('cd %s/%s\n' % (os.getcwd(), pwd_iqtree_guide_tree_wd))
                pwd_iqtree_js_handle.write(get_guide_tree_cmd + '\n\n')
                pwd_iqtree_js_handle.write('cd %s/%s\n' % (os.getcwd(), pwd_iqtree_c60_pmsf_wd))
                pwd_iqtree_js_handle.write(get_c60_tree_cmd + '\n')

            if submit_job is True:
                print(pwd_iqtree_js)
                submit_js(pwd_iqtree_js)
            else:
                print('Job script for running iqtree exported to %s' % pwd_iqtree_js)

    # prepare files for performing dating
    print('Preparing files for performing dating')



if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-deltall', required=True,                          help='DeltaLL stdout')
    parser.add_argument('-o',       required=True,                          help='output dir')
    parser.add_argument('-c',       required=False, default='25-50-75-100', help='cutoffs, default: 25-50-75-100')
    parser.add_argument('-mmn',     required=False, default=20, type=int,   help='minimal marker number, default: 20')
    parser.add_argument('-aln',     required=True,                          help='faa file dir')
    parser.add_argument('-jst',     required=False, default='6',            help='threads to request in job script, for running iqtree')
    parser.add_argument('-qsub',    required=False, action="store_true",    help='submit job scripts')
    parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    AssessMarkerDeltaLL(args)
