import os
import argparse
import itertools


dating_usage = '''
======================== Dating example commands ========================

# requireMENT: mcmctree, baseml and codeml (all from PAML)

# example commands
TreeSAK Dating -tree tree.treefile -msa markers.phylip -o dating_wd -f

=========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def prep_mcmctree_ctl(ctl_para_dict, mcmctree_ctl_file):

    ctl_file_handle = open(mcmctree_ctl_file, 'w')
    ctl_file_handle.write('          seed = %s\n' % ctl_para_dict.get('seed', '-1'))
    ctl_file_handle.write('       seqfile = %s\n' % ctl_para_dict['seqfile'])
    ctl_file_handle.write('      treefile = %s\n' % ctl_para_dict['treefile'])
    ctl_file_handle.write('      mcmcfile = %s\n' % ctl_para_dict['mcmcfile'])
    ctl_file_handle.write('       outfile = %s\n' % ctl_para_dict['outfile'])
    ctl_file_handle.write('         ndata = %s\n'                                                                           % ctl_para_dict.get('ndata',        1))
    ctl_file_handle.write('       seqtype = %s    	* 0: nucleotides; 1:codons; 2:AAs\n'                                    % ctl_para_dict['seqtype'])
    ctl_file_handle.write('       usedata = %s    	* 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)\n'   % ctl_para_dict['usedata'])
    ctl_file_handle.write('         clock = %s    	* 1: global clock; 2: independent rates; 3: correlated rates\n'         % ctl_para_dict.get('clock',        2))
    ctl_file_handle.write('       RootAge = %s      * safe constraint on root age, used if no fossil for root.\n'           % ctl_para_dict.get('RootAge',      '<1.0'))
    ctl_file_handle.write('         model = %s    	* 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n'                               % ctl_para_dict.get('model',        0))
    ctl_file_handle.write('         alpha = %s  	* alpha for gamma rates at sites\n'                                     % ctl_para_dict.get('alpha',        0.5))
    ctl_file_handle.write('         ncatG = %s    	* No. categories in discrete gamma\n'                                   % ctl_para_dict.get('ncatG',        4))
    ctl_file_handle.write('     cleandata = %s    	* remove sites with ambiguity data (1:yes, 0:no)?\n'                    % ctl_para_dict.get('cleandata',    0))
    ctl_file_handle.write('       BDparas = %s      * birth, death, sampling\n'                                             % ctl_para_dict.get('BDparas',      '1 1 0.1'))
    ctl_file_handle.write('   kappa_gamma = %s      * gamma prior for kappa\n'                                              % ctl_para_dict.get('kappa_gamma',  '6 2'))
    ctl_file_handle.write('   alpha_gamma = %s      * gamma prior for alpha\n'                                              % ctl_para_dict.get('alpha_gamma',  '1 1'))
    ctl_file_handle.write('   rgene_gamma = %s      * gammaDir prior for rate for genes\n'                                  % ctl_para_dict.get('rgene_gamma',  '1 50 1'))
    ctl_file_handle.write('  sigma2_gamma = %s      * gammaDir prior for sigma^2     (for clock=2 or 3)\n'                  % ctl_para_dict.get('sigma2_gamma', '1 10 1'))
    ctl_file_handle.write('      finetune = %s      * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr\n'    % ctl_para_dict.get('finetune',     '1: .1 .1 .1 .1 .1 .1'))
    ctl_file_handle.write('         print = %s      * 0: no mcmc sample; 1: everything except branch rates 2: everything\n' % ctl_para_dict.get('print',        1))
    ctl_file_handle.write('        burnin = %s\n'                                                                           % ctl_para_dict.get('burnin',       50000))
    ctl_file_handle.write('      sampfreq = %s\n'                                                                           % ctl_para_dict.get('sampfreq',     50))
    ctl_file_handle.write('       nsample = %s\n'                                                                           % ctl_para_dict.get('nsample',      10000))
    ctl_file_handle.close()


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


def dating(args):

    op_dir              = args['o']
    tree_file           = args['tree']
    msa_file            = args['msa']
    force_overwrite     = args['f']
    seq_type            = args['st']
    settings_to_compare = args['s']

    para_to_test_dict = dict()
    for each_para in open(settings_to_compare):
        each_para_split = each_para.strip().split()
        para_list = each_para_split[1].split(',')
        para_to_test_dict[each_para_split[0]] = para_list

    ####################################################################################################################

    tree_f_name, tree_f_path, tree_f_base, tree_f_ext = sep_path_basename_ext(tree_file)
    msa_f_name,  msa_f_path,  msa_f_base,  msa_f_ext  = sep_path_basename_ext(msa_file)

    get_bv_wd         = '%s/get_bv_wd'         % op_dir
    mcmctree_ctl_bv   = '%s/mcmctree.ctl'      % get_bv_wd
    mcmctree_cmds_txt = '%s/mcmctree_cmds.txt' % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder exist, program exited!')
            exit()

    os.system('mkdir %s' % op_dir)

    ####################################################################################################################

    # prepare files for getting bv file
    os.system('mkdir %s'  % get_bv_wd)
    os.system('cp %s %s/' % (tree_file, get_bv_wd))
    os.system('cp %s %s/' % (msa_file, get_bv_wd))

    get_bv_para_dict = dict()
    get_bv_para_dict['seqfile']  = msa_f_name
    get_bv_para_dict['treefile'] = tree_f_name
    get_bv_para_dict['mcmcfile'] = 'mcmc.txt'
    get_bv_para_dict['outfile']  = 'out.txt'
    get_bv_para_dict['seqtype']  = seq_type
    get_bv_para_dict['usedata']  = '3'

    prep_mcmctree_ctl(get_bv_para_dict, mcmctree_ctl_bv)

    # write out command
    mcmctree_cmds_txt_handle = open(mcmctree_cmds_txt, 'w')
    mcmctree_cmds_txt_handle.write('cd %s; mcmctree\n' % get_bv_wd.split('/')[-1])
    mcmctree_cmds_txt_handle.write('\n# Please wait for the above command to be finished before you move to the following ones!\n\n')
    mcmctree_cmds_txt_handle.close()

    ####################################################################################################################

    para_comb_dict = get_parameter_combinations(para_to_test_dict)

    mcmctree_cmds_txt_handle = open(mcmctree_cmds_txt, 'a')
    for para_comb in sorted(list(para_comb_dict.keys())):

        # create dir
        current_dating_wd_1 = '%s/%s_run1'  % (op_dir, para_comb)
        current_dating_wd_2 = '%s/%s_run2'  % (op_dir, para_comb)
        os.system('mkdir %s' % current_dating_wd_1)
        os.system('mkdir %s' % current_dating_wd_2)

        # copy tree and msa file
        os.system('cp %s %s/' % (tree_file, current_dating_wd_1))
        os.system('cp %s %s/' % (tree_file, current_dating_wd_2))
        os.system('cp %s %s/' % (msa_file,  current_dating_wd_1))
        os.system('cp %s %s/' % (msa_file,  current_dating_wd_2))

        # prepare mcmctree.ctl file
        mcmctree_ctl_1  = '%s/mcmctree.ctl' % current_dating_wd_1
        mcmctree_ctl_2  = '%s/mcmctree.ctl' % current_dating_wd_2

        current_para_dict = para_comb_dict[para_comb]
        current_para_dict['seqfile']  = msa_f_name
        current_para_dict['treefile'] = tree_f_name
        current_para_dict['mcmcfile'] = '%s_mcmc.txt' % para_comb
        current_para_dict['outfile']  = '%s_out.txt'  % para_comb
        current_para_dict['seqtype']  = seq_type
        current_para_dict['usedata']  = '2'

        prep_mcmctree_ctl(current_para_dict, mcmctree_ctl_1)
        prep_mcmctree_ctl(current_para_dict, mcmctree_ctl_2)

        # write out commands
        mcmctree_cmds_txt_handle.write('cd %s; cp ../get_bv_wd/out.BV in.BV; mcmctree\n' % (current_dating_wd_1.split('/')[-1]))
        mcmctree_cmds_txt_handle.write('cd %s; cp ../get_bv_wd/out.BV in.BV; mcmctree\n' % (current_dating_wd_2.split('/')[-1]))
    mcmctree_cmds_txt_handle.close()

    print('Job script for performing dating exported to %s' % mcmctree_cmds_txt)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-tree',    required=True,                          help='outgroup leaves, one id per line')
    parser.add_argument('-msa',     required=True,                          help='EU tree with time constraints')
    parser.add_argument('-o',       required=True,                          help='dating wd')
    parser.add_argument('-s',       required=True,                          help='settings to compare')
    parser.add_argument('-st',      required=False, default='2',            help='sequence type, 0 for nucleotides, 1 for codons, 2 for AAs, default: 2')
    parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    dating(args)

