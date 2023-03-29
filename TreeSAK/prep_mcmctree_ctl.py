import itertools


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


mcmctree_ctl_dict = {'seqfile' : 'concatenated.phy',
                     'treefile': 'deltall75_pa75_rooted_with_calibrations.nwk',
                     'mcmcfile': 'mcmc.txt',
                     'outfile' : 'DateArTree_out.txt',
                     'seqtype' : 2,
                     'usedata' : 3,
                     'clock'   : 3}


prep_mcmctree_ctl(mcmctree_ctl_dict, '/Users/songweizhi/Desktop/aaa.txt')


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


para_to_test_dict = {'clock': [2, 3], 'nsample': [20000, 50000], 'model': [0, 4], 'kappa_gamma': ['6 2', '5 1']}
para_dod = get_parameter_combinations(para_to_test_dict)
print(para_dod)

# all_combination_list_in_str = ['_'.join(i) for i in all_combination_list]
# print(all_combination_list_in_str)
# print(len(all_combination_list_in_str))




