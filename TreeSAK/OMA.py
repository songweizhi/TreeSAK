import os
import glob
import argparse


OMA_usage = '''
======================= OMA example commands =======================

TreeSAK OMA -i faa_files -x faa -og og_gnm.txt -o OMA_wd -f -t 32

====================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def get_default_para_dict():

    default_para_str = '''
    OutputFolder := 'Output';
    ReuseCachedResults := true;
    AlignBatchSize := 1e6;
    MinScore := 181;
    LengthTol := 0.61;
    StablePairTol := 1.81;
    InparalogTol := 3.00;
    ParalogTol := -2.5*StablePairTol;
    VerifiedPairTol := 1.53;
    MinSeqLen := 50;
    UseOnlyOneSplicingVariant := true;
    UseExperimentalHomologousClusters := false;
    QuasiCliquesCutoff := 1.0:
    StableIdsForGroups := false;
    GuessIdType := false;
    DoHierarchicalGroups := 'bottom-up';
    SpeciesTree := 'estimate';
    MinEdgeCompletenessFraction := 0.65;
    ReachabilityCutoff := 0.65;
    MaxTimePerLevel := 1200;  # 20min
    DoGroupFunctionPrediction := true;
    GroupFunctionCutoff := 0.5;
    CladeDefinition := 'default';
    UseEsprit := false;
    DistConfLevel := 2;
    MinProbContig := 0.4;
    MaxContigOverlap := 5;
    MinSeqLenContig := 20;
    MinBestScore := 250;
    '''

    default_para_dict = dict()
    for each_line in default_para_str.split('    '):
        para_line = each_line.replace(' ', '').replace('\n', '').split(';')[0]
        if para_line != '':
            para_line_split = para_line.split(':=')
            default_para_dict[para_line_split[0]] = para_line_split[1]

    return default_para_dict


def OMA(args):

    gnm_dir         = args['i']
    file_ext        = args['x']
    seq_type        = args['st']
    og_gnm_txt      = args['og']
    op_dir          = args['o']
    force_overwrite = args['f']
    num_threads     = args['t']

    # define file name
    pwd_gnm_rename_txt = '%s/rename.txt'     % op_dir
    pwd_parameter_file = '%s/parameters.drw' % op_dir
    oma_input_dir      = '%s/DB'             % op_dir

    # create dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % oma_input_dir)

    # check genome files
    gnm_file_re   = '%s/*.%s' % (gnm_dir, file_ext)
    gnm_file_list = glob.glob(gnm_file_re)
    if len(gnm_file_list) == 0:
        print('No genome detected, program exited!')
        exit()

    # check og_gnm_txt
    if os.path.isfile(og_gnm_txt) is False:
        print('Out group genome id file not detected, program exited!')
        exit()

    # copy genome files into DB folder
    gnm_id_rename_dict = dict()
    rename_list = []
    for each_gnm in gnm_file_list:
        gnm_path, gnm_base, gnm_ext = sep_path_basename_ext(each_gnm)
        gnm_base_renamed = gnm_base.replace('.', '_')
        pwd_gnm_db = '%s/%s.fa' % (oma_input_dir, gnm_base_renamed)
        if gnm_base != gnm_base_renamed:
            rename_list.append('%s\t%s' % (gnm_base, gnm_base_renamed))
        gnm_id_rename_dict[gnm_base] = gnm_base_renamed
        os.system('cp %s %s' % (each_gnm, pwd_gnm_db))

    # write out rename file
    if len(rename_list) > 0:
        pwd_gnm_rename_txt_handle = open(pwd_gnm_rename_txt, 'w')
        for each_e in sorted(rename_list):
            pwd_gnm_rename_txt_handle.write(each_e + '\n')
        pwd_gnm_rename_txt_handle.close()
    else:
        print('Format of file names passed checking')

    # get default_para_dict
    default_para_dict = get_default_para_dict()

    # read in og_gnm_txt
    renamed_og_gnm_list = []
    for each_og_gnm in open(og_gnm_txt):
        og_gnm_renamed = gnm_id_rename_dict[each_og_gnm.strip()]
        renamed_og_gnm_list.append(og_gnm_renamed)

    # write out parameter file
    with open(pwd_parameter_file, 'w') as pwd_parameter_file_handle:

        # write InputDataType line
        if seq_type in ['AA', 'aa', 'Aa']:
            pwd_parameter_file_handle.write("InputDataType := 'AA';\n")
        if seq_type in ['DNA', 'dna', 'Dna']:
            pwd_parameter_file_handle.write("InputDataType := 'DNA';\n")

        # write OutgroupSpecies line
        OutgroupSpecies_value_str = "['%s']" % "', '".join(renamed_og_gnm_list)
        pwd_parameter_file_handle.write("OutgroupSpecies := %s;\n" % OutgroupSpecies_value_str)

        # write out the rest lines
        for each_para in default_para_dict:
            para_value = default_para_dict[each_para]
            pwd_parameter_file_handle.write("%s := %s;\n" % (each_para, para_value))

    # final report
    print('You can run OMA with:')
    print('cd %s' % op_dir)
    print('oma -n %s' % num_threads)
    print('# You may want to customize parameters specified in %s ' % pwd_parameter_file)


if __name__ == '__main__':

    OMA_parser = argparse.ArgumentParser()
    OMA_parser.add_argument('-i',   required=True,                       help='genome folder')
    OMA_parser.add_argument('-x',   required=True,                       help='genome file extension')
    OMA_parser.add_argument('-st',  required=False, default='AA',        help='sequence type, AA or DNA, default: AA')
    OMA_parser.add_argument('-og',  required=True,                       help='outgroup genomes, without file extension')
    OMA_parser.add_argument('-o',   required=True,  default=None,        help='output dir, i.e., OMA working directory')
    OMA_parser.add_argument('-f',   required=False, action="store_true", help='force overwrite')
    OMA_parser.add_argument('-t',   required=False, type=int, default=6, help='number of threads for running OMA, default: 6')
    args = vars(OMA_parser.parse_args())
    OMA(args)
