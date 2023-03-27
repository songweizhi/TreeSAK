import os
import argparse


Dating_usage = '''
============================= Dating example commands =============================

Dependencies:

# example commands
TreeSAK Dating -deltall DeltaLL_stdout.txt -aod s11_marker_sets_by_DeltaLL -o s12_dating_wd -c 25-50-75-100 -mmn 20 -f

===================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def submit_js(js):
    current_wd = os.getcwd()
    js_path, js_basename, js_ext = sep_path_basename_ext(js)
    os.chdir(js_path)
    os.system('qsub %s%s' % (js_basename, js_ext))
    os.chdir(current_wd)


def Dating(args):

    deltall_stdout_txt      = args['deltall']
    aod                     = args['aod']
    op_dir                  = args['o']
    deltall_keep_pct_str    = args['c']
    min_marker_num          = args['mmn']
    js_cpu_num              = args['jst']
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
    marker_list_sorted_by_deltall = [k for k, v in sorted(overall_score_dict.items(), key=lambda item: item[1])]

    # get qualified marker list
    for each_keep_pct in deltall_keep_pct_list:
        marker_num_to_keep = round(len(marker_list_sorted_by_deltall)*each_keep_pct/100)

        if marker_num_to_keep < min_marker_num:
            print('Ignored DeltaLL cutoff at %s , the number of qualified markers (%s) less than %s' % (each_keep_pct, marker_num_to_keep, min_marker_num))
        else:
            pwd_aln_concatenated       = '%s/%s_DeltaLL_%s_concatenated.phy'                         % (aod, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file          = '%s/%s_DeltaLL_%s_iqtree_C60_PMSF/concatenated.treefile'    % (aod, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            current_dating_wd          = '%s/%s_DeltaLL_%s'                                          % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file_renamed  = '%s/concatenated_raw.treefile'                              % current_dating_wd
            pwd_to_do_note_txt         = '%s/to_do.txt'                                              % current_dating_wd
            pwd_dating_js              = '%s/js_%s_DeltaLL_%s_dating.sh'                             % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)

            print('Alignment:\t%s' % pwd_aln_concatenated)
            print('Tree file:\t%s' % pwd_c60_tree_file)
            print('Dating wd:\t%s' % current_dating_wd)
            print()

            # create dating wd and copy file into it
            os.mkdir(current_dating_wd)
            os.system('cp %s %s' % (pwd_c60_tree_file, pwd_c60_tree_file_renamed))
            os.system('cp %s %s/' % (pwd_aln_concatenated, current_dating_wd))

            # write out to do
            with open(pwd_to_do_note_txt, 'w') as pwd_to_do_note_txt_handle:
                pwd_to_do_note_txt_handle.write('1. Root the tree in iTOL\n')
                pwd_to_do_note_txt_handle.write('2. Replace the eukaryotic clade with the one has EU time constraints using TreeGraph2\n')
                pwd_to_do_note_txt_handle.write('3. Add time constraints to the tree\n')
                pwd_to_do_note_txt_handle.write('4. Renamed the modified tree to concatenated_modified.treefile\n')

            # write job script
            with open(pwd_dating_js, 'w') as pwd_dating_js_handle:
                dating_cmd = ''
                pwd_dating_js_handle.write('cd %s/%s\n' % (os.getcwd(), current_dating_wd))
                pwd_dating_js_handle.write(dating_cmd + '\n')
            print('Job script for dating exported to %s' % pwd_dating_js)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-deltall', required=True,                          help='DeltaLL stdout')
    parser.add_argument('-aod',     required=True,                          help='AssessMarkerDeltaLL output dir')
    parser.add_argument('-o',       required=True,                          help='dating wd')
    parser.add_argument('-c',       required=False, default='25-50-75-100', help='cutoffs, default: 25-50-75-100')
    parser.add_argument('-mmn',     required=False, default=20, type=int,   help='minimal marker number, default: 20')
    parser.add_argument('-jst',     required=False, default='6',            help='threads to request in job script, for performing dating')
    parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    Dating(args)


'''

cd /Users/songweizhi/Desktop/dating_test
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/Dating.py -deltall Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s10_assess_marker_deltaLL/PA_75_DeltaLL_stdout.txt -aod Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s11_marker_sets_by_DeltaLL -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s12_dating_wd -c 25-50-75-100 -mmn 20 -f

'''
