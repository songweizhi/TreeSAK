import os
import argparse
import itertools
from ete3 import Tree


Dating_usage = '''
============================= Dating example commands =============================

Dependencies:

# example commands
TreeSAK Dating -deltall DeltaLL_stdout.txt -aod s11_marker_sets_by_DeltaLL -o s12_dating_wd -c 25-50-75-100 -mmn 20 -f

===================================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def submit_js(js):
    current_wd = os.getcwd()
    js_path, js_basename, js_ext = sep_path_basename_ext(js)
    os.chdir(js_path)
    os.system('qsub %s%s' % (js_basename, js_ext))
    os.chdir(current_wd)


def root_with_out_group(tree_file, out_group_txt, tree_file_rooted):

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=1)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted)


def replace_clades(main_tree, sub_tree, tree_out):

    tre_sub = Tree(sub_tree, format=1)
    subtree_leaf_name_list = tre_sub.get_leaf_names()
    tre_main = Tree(main_tree)
    lca = tre_main.get_common_ancestor(subtree_leaf_name_list)

    if len(lca.get_leaf_names()) != len(subtree_leaf_name_list):
        print('LCA of subtree leaves in main tree contain extra leaves, program exited!')
        exit()

    lca_p = lca.up
    lca_p.remove_child(lca)
    lca_p.add_child(tre_sub)
    tre_main.write(outfile=tree_out, format=8)


def prep_mcmctree_ctl(ctl_para_dict, mcmctree_ctl_file):
    with open(mcmctree_ctl_file, 'w') as ctl_file_handle:
        ctl_file_handle.write('      finetune = %s\n' % ctl_para_dict.get('seed',           '-1'))
        ctl_file_handle.write('       seqfile = %s\n' % ctl_para_dict['seqfile'])
        ctl_file_handle.write('      treefile = %s\n' % ctl_para_dict['treefile'])
        ctl_file_handle.write('      mcmcfile = %s\n' % ctl_para_dict['mcmcfile'])
        ctl_file_handle.write('       outfile = %s\n' % ctl_para_dict['outfile'])
        ctl_file_handle.write('         ndata = %s\n' % ctl_para_dict.get('ndata',          1))
        ctl_file_handle.write('       seqtype = %s\n' % ctl_para_dict['seqtype'])
        ctl_file_handle.write('       usedata = %s\n' % ctl_para_dict['usedata'])
        ctl_file_handle.write('         clock = %s\n' % ctl_para_dict['clock'])
        ctl_file_handle.write('       RootAge = %s\n' % ctl_para_dict.get('RootAge',        '<1.0'))
        ctl_file_handle.write('         model = %s\n' % ctl_para_dict.get('model',          0))
        ctl_file_handle.write('         alpha = %s\n' % ctl_para_dict.get('alpha',          0.5))
        ctl_file_handle.write('         ncatG = %s\n' % ctl_para_dict.get('ncatG',          4))
        ctl_file_handle.write('     cleandata = %s\n' % ctl_para_dict.get('cleandata',      0))
        ctl_file_handle.write('       BDparas = %s\n' % ctl_para_dict.get('BDparas',        '1 1 0.1'))
        ctl_file_handle.write('   kappa_gamma = %s\n' % ctl_para_dict.get('kappa_gamma',    '6 2'))
        ctl_file_handle.write('   alpha_gamma = %s\n' % ctl_para_dict.get('alpha_gamma',    '1 1'))
        ctl_file_handle.write('   rgene_gamma = %s\n' % ctl_para_dict.get('rgene_gamma',    '1 50 1'))
        ctl_file_handle.write('  sigma2_gamma = %s\n' % ctl_para_dict.get('sigma2_gamma',   '1 10 1'))
        ctl_file_handle.write('      finetune = %s\n' % ctl_para_dict.get('finetune',       '1: .1 .1 .1 .1 .1 .1'))
        ctl_file_handle.write('         print = %s\n' % ctl_para_dict.get('print',          1))
        ctl_file_handle.write('        burnin = %s\n' % ctl_para_dict.get('burnin',         50000))
        ctl_file_handle.write('      sampfreq = %s\n' % ctl_para_dict.get('sampfreq',       5))
        ctl_file_handle.write('       nsample = %s\n' % ctl_para_dict.get('nsample',        50000))


def get_parameter_combinations(para_to_test_dict):

    para_lol_name = []
    para_lol_value = []
    para_lol_name_with_value = []
    for each_para in sorted(list(para_to_test_dict.keys())):
        para_setting_list_name = []
        para_setting_list_value = []
        para_setting_list_name_with_value = []
        for each_setting in sorted(para_to_test_dict[each_para]):
            name_str = ('%s%s' % (each_para, each_setting)).replace(' ', '_')
            para_setting_list_name.append(each_para)
            para_setting_list_value.append(each_setting)
            para_setting_list_name_with_value.append(name_str)
        para_lol_name.append(para_setting_list_name)
        para_lol_value.append(para_setting_list_value)
        para_lol_name_with_value.append(para_setting_list_name_with_value)

    all_combination_list_name = [p for p in itertools.product(*para_lol_name)]
    all_combination_list_value = [p for p in itertools.product(*para_lol_value)]
    all_combination_list_name_with_value = [p for p in itertools.product(*para_lol_name_with_value)]
    all_combination_list_name_with_value_str = ['_'.join(i) for i in all_combination_list_name_with_value]

    para_dod = dict()
    element_index = 0
    for each_combination in all_combination_list_name_with_value_str:
        current_name_list   = all_combination_list_name[element_index]
        current_value_list  = all_combination_list_value[element_index]
        current_para_dict = dict()
        for key, value in zip(current_name_list, current_value_list):
            current_para_dict[key] = value
        para_dod[each_combination] = current_para_dict
        element_index += 1

    return para_dod


def Dating(args):

    deltall_stdout_txt      = args['deltall']
    aod                     = args['aod']
    out_group_txt           = args['og']
    eu_tree                 = args['eu']
    op_dir                  = args['o']
    deltall_keep_pct_str    = args['c']
    min_marker_num          = args['mmn']
    js_cpu_num              = args['jst']
    force_overwrite         = args['f']
    root_age                = args['ra']
    submit_job              = args['qsub']

    if os.path.isfile(eu_tree) is False:
        print('%s not found, program exited!' % eu_tree)
        exit()

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
            prefix_base                              = '%s_DeltaLL_%s'                                          % (deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            aln_concatenated                         = '%s_DeltaLL_%s_concatenated.phy'                         % (deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            c60_tree_file_rooted_with_time_final     = '%s_DeltaLL_%s_rooted_with_time_final.treefile'          % (deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file                        = '%s/%s_DeltaLL_%s_iqtree_C60_PMSF/concatenated.treefile' % (aod, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file_renamed                = '%s/%s_DeltaLL_%s_raw.treefile'                          % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file_rooted                 = '%s/%s_DeltaLL_%s_rooted.treefile'                       % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_c60_tree_file_rooted_with_time       = '%s/%s_DeltaLL_%s_rooted_with_time.treefile'             % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct)
            pwd_aln_concatenated                     = '%s/%s'                                                  % (aod, aln_concatenated)
            pwd_c60_tree_file_rooted_with_time_final = '%s/%s'                                                  % (op_dir, c60_tree_file_rooted_with_time_final)
            get_BV_wd                                = '%s/%s_get_BV_wd'                                        % (op_dir, prefix_base)

            os.system('cp %s %s/' % (pwd_aln_concatenated, op_dir))
            os.system('cp %s %s'  % (pwd_c60_tree_file, pwd_c60_tree_file_renamed))

            # root genome tree with outgroup
            root_with_out_group(pwd_c60_tree_file_renamed, out_group_txt, pwd_c60_tree_file_rooted)

            # add time constraints
            replace_clades(pwd_c60_tree_file_rooted, eu_tree, pwd_c60_tree_file_rooted_with_time)

            # remove "NoName" from the rooted tree with time constraints
            tree_str = open(pwd_c60_tree_file_rooted_with_time).readline().strip().replace('NoName', '')

            # add root age
            tree_str = tree_str.replace(';', '<%s;' % root_age)
            tre_object = Tree(tree_str, format=1)
            with open(pwd_c60_tree_file_rooted_with_time_final, 'w') as pwd_c60_tree_file_rooted_with_time_final_hanlde:
                pwd_c60_tree_file_rooted_with_time_final_hanlde.write('%s\t1\n' % len(tre_object.get_leaf_names()))
                pwd_c60_tree_file_rooted_with_time_final_hanlde.write(tree_str + '\n')

            # rm tmp tree files
            os.system('rm %s' % pwd_c60_tree_file_renamed)
            os.system('rm %s' % pwd_c60_tree_file_rooted)
            os.system('rm %s' % pwd_c60_tree_file_rooted_with_time)

            # get BV file
            os.mkdir(get_BV_wd)
            os.system('cp %s %s/' % (pwd_aln_concatenated, get_BV_wd))
            os.system('cp %s %s/' % (pwd_c60_tree_file_rooted_with_time_final, get_BV_wd))

            get_BV_js                       = '%s/%s_get_BV.sh'         % (op_dir, prefix_base)
            get_BV_mcmctree_ctl             = '%s_get_BV_mcmctree.ctl'  % (prefix_base)
            pwd_get_BV_mcmctree_ctl         = '%s/%s'                   % (get_BV_wd, get_BV_mcmctree_ctl)

            get_BV_para_dict = dict()
            get_BV_para_dict['seqfile']     = aln_concatenated
            get_BV_para_dict['treefile']    = c60_tree_file_rooted_with_time_final
            get_BV_para_dict['mcmcfile']    = '%s_mcmc.txt'             % prefix_base
            get_BV_para_dict['outfile']     = '%s_out.txt'              % prefix_base
            get_BV_para_dict['seqtype']     = '2'
            get_BV_para_dict['usedata']     = '3'
            get_BV_para_dict['clock']       = '3'
            prep_mcmctree_ctl(get_BV_para_dict, pwd_get_BV_mcmctree_ctl)

            with open(get_BV_js, 'w') as get_BV_js_handle:
                get_BV_js_handle.write('#!/bin/bash\n\n')
                get_BV_js_handle.write('cd %s/%s\n' % (os.getcwd(), get_BV_wd))
                get_BV_js_handle.write('mcmctree %s\n' % get_BV_mcmctree_ctl)

            # prepare files for dating
            para_to_test_dict = {'clock': [2, 3], 'model': [0, 4]}
            para_to_test_dict = {'clock': [2, 3]}
            para_dod = get_parameter_combinations(para_to_test_dict)
            for para_combination in para_dod:
                mcmctree_ctl        = '%s_%s_mcmctree.ctl'              % (prefix_base, para_combination)
                current_dating_wd   = '%s/%s_DeltaLL_%s_%s_dating_wd'   % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct, para_combination)
                pwd_mcmctree_ctl    = '%s/%s_%s_mcmctree.ctl'           % (current_dating_wd, prefix_base, para_combination)
                js_mcmctree         = '%s/%s_DeltaLL_%s_%s.sh'          % (op_dir, deltall_stdout_basename.split('_DeltaLL_stdout')[0], each_keep_pct, para_combination)

                # create dating wd and copy tree and alignment files into it
                os.mkdir(current_dating_wd)
                os.system('cp %s %s/'       % (pwd_aln_concatenated, current_dating_wd))
                os.system('cp %s %s/'       % (pwd_c60_tree_file_rooted_with_time_final, current_dating_wd))

                current_para_dict = para_dod[para_combination]
                current_para_dict['seqfile']  = aln_concatenated
                current_para_dict['treefile'] = c60_tree_file_rooted_with_time_final
                current_para_dict['mcmcfile'] = '%s_%s_mcmc.txt' % (prefix_base, para_combination)
                current_para_dict['outfile']  = '%s_%s_out.txt'  % (prefix_base, para_combination)
                current_para_dict['seqtype']  = '2'
                current_para_dict['usedata']  = '2'

                prep_mcmctree_ctl(current_para_dict, pwd_mcmctree_ctl)

                with open(js_mcmctree, 'w') as js_mcmctree_handle:
                    js_mcmctree_handle.write('#!/bin/bash\n\n')
                    js_mcmctree_handle.write('cd %s/%s\n' % (os.getcwd(), current_dating_wd))
                    js_mcmctree_handle.write('cp ../%s_get_BV_wd/out.BV in.BV\n' % prefix_base)
                    js_mcmctree_handle.write('mcmctree %s\n' % mcmctree_ctl)
                print('Job script for performing dating exported to %s' % js_mcmctree)

            if submit_job is True:
                submit_js(get_BV_js)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-deltall', required=True,                          help='DeltaLL stdout')
    parser.add_argument('-aod',     required=True,                          help='AssessMarkerDeltaLL output dir')
    parser.add_argument('-og',      required=True,                          help='outgroup leaves, one leaf id per line')
    parser.add_argument('-eu',      required=True,                          help='EU tree with time constraints')
    parser.add_argument('-o',       required=True,                          help='dating wd')
    parser.add_argument('-c',       required=False, default='25-50-75-100', help='cutoffs, default: 25-50-75-100')
    parser.add_argument('-mmn',     required=False, default=20, type=int,   help='minimal marker number, default: 20')
    parser.add_argument('-ra',      required=False, default=45, type=int,   help='root age, default: 45')
    parser.add_argument('-jst',     required=False, default='6',            help='threads to request in job script, for performing dating')
    parser.add_argument('-qsub',    required=False, action="store_true",    help='submit job scripts for getting in.BV')
    parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    Dating(args)


'''

cd /Users/songweizhi/Desktop/dating_test
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/Dating.py -deltall Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s10_assess_marker_deltaLL/PA_75_DeltaLL_stdout.txt -aod Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s11_marker_sets_by_DeltaLL -og out_group.txt -eu 27.nwk -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s12_dating_wd -c 25-50-75-100 -mmn 20 -f 

cd /home-user/wzsong/DateArTree
python3 Dating.py -deltall Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s10_assess_marker_deltaLL/PA_75_DeltaLL_stdout.txt -aod Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s11_marker_sets_by_DeltaLL -og out_group.txt -eu 27.nwk -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30/s12_dating_wd -c 25-50-75-100 -mmn 20 -f 


'''
