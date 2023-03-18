import os
import glob
import argparse
from Bio import SeqIO


AssessMarkerPA_usage = '''
=========================== AssessMarkerPA example commands ===========================

BioSAK AssessMarkerPA -ta trimmed_aln -tax aln -aa faa_files -aax faa -g gnm_group.txt -c 25-50-75-100 -o s10_assess_marker_PA -pl catfasta2phyml.pl

Note
1. Extra genomes in gnm_metadata.txt won't affect assessment results.
2. Genomes can not be found in gnm_metadata.txt will trigger an error.
3. Alignments in {trimmed_aln_dir} need to be trimmed before assessment
4. Sequences in MSAs need to be named by genome id.

=======================================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def AssessMarkerPA(args):

    trimmed_aln_dir   = args['ta']
    trimmed_aln_ext   = args['tax']
    faa_file_dir      = args['aa']
    faa_file_ext      = args['aax']
    gnm_group_txt     = args['g']
    cutoff_str        = args['c']
    op_dir            = args['o']
    force_overwriting = args['f']
    catfasta2phyml_pl = args['pl']

    # get gnm id list
    faa_file_re   = '%s/*.%s' % (faa_file_dir, faa_file_ext)
    faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
    gnm_set = set()
    for each_faa_file in faa_file_list:
        faa_path, faa_basename, faa_ext = sep_path_basename_ext(each_faa_file)
        gnm_set.add(faa_basename)

    # read in genome metadata
    group_to_gnm_dict       = dict()
    group_to_gnm_num_dict   = dict()
    gnm_to_group_dict       = dict()
    for each_gnm in open(gnm_group_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        domain_name = each_gnm_split[1]

        if gnm_id in gnm_set:
            gnm_to_group_dict[gnm_id] = domain_name

            if domain_name not in group_to_gnm_num_dict:
                group_to_gnm_num_dict[domain_name] = 1
            else:
                group_to_gnm_num_dict[domain_name] += 1

            if domain_name not in group_to_gnm_dict:
                group_to_gnm_dict[domain_name] = {gnm_id}
            else:
                group_to_gnm_dict[domain_name].add(gnm_id)

    group_id_list_sorted = sorted(list(group_to_gnm_dict.keys()))

    # exit program if group information is missing
    gnms_without_group_info = set()
    for gnm in gnm_set:
        if gnm not in gnm_to_group_dict:
            gnms_without_group_info.add(gnm)

    if len(gnms_without_group_info) > 0:
        print('Group information for the following genomes are missing from %s, program exited!' % gnm_group_txt)
        print(','.join(gnms_without_group_info))
        print('Group information for the above genomes are missing from %s, program exited!' % gnm_group_txt)
        exit()

    # create folder
    if force_overwriting is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
    else:
        if os.path.isdir(op_dir) is True:
            print('Output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # read in provided cutoffs
    present_pct_cutoff_list = [int(i) for i in cutoff_str.split('-')]
    assess_summary_1_txt    = '%s/assessment_PA.txt'                % op_dir
    assess_summary_2_txt    = '%s/assessment_PA_summary.txt'        % op_dir
    itol_binary_txt         = '%s/assessment_PA_iTOL_binary.txt'    % op_dir

    trimmed_aln_file_re = '%s/*.%s' % (trimmed_aln_dir, trimmed_aln_ext)
    trimmed_aln_file_list = [os.path.basename(file_name) for file_name in glob.glob(trimmed_aln_file_re)]

    assess_summary_1_txt_handle = open(assess_summary_1_txt, 'w')
    assess_summary_1_txt_handle.write('Marker\t%s\n' % '\t'.join([str(i) for i in group_id_list_sorted]))
    assess_summary_2_txt_handle = open(assess_summary_2_txt, 'w')
    assess_summary_2_txt_handle.write('Marker\t%s\n' % '\t'.join([str(i) for i in present_pct_cutoff_list]))
    cutoff_to_qualified_marker_dict = dict()
    gnm_to_identified_marker_dict = dict()
    marker_id_list = []
    for each_aln in trimmed_aln_file_list:

        marker_id = each_aln.split(('.%s' % trimmed_aln_ext))[0]
        marker_id_list.append(marker_id)
        pwd_aln = '%s/%s' % (trimmed_aln_dir, each_aln)

        current_marker_num_by_group_dict = dict()
        for each_seq in SeqIO.parse(pwd_aln, 'fasta'):
            gnm_id = each_seq.id

            # get genome to marker dist
            if gnm_id not in gnm_to_identified_marker_dict:
                gnm_to_identified_marker_dict[gnm_id] = {marker_id}
            else:
                gnm_to_identified_marker_dict[gnm_id].add(marker_id)

            if gnm_id in gnm_to_group_dict:
                gnm_group = gnm_to_group_dict[gnm_id]
            else:
                print('Not all genomes used to generate the MSA being found in -aa, program exited!')
                exit()

            if gnm_group not in current_marker_num_by_group_dict:
                current_marker_num_by_group_dict[gnm_group] = 1
            else:
                current_marker_num_by_group_dict[gnm_group] += 1

        # write out assess_summary_1_txt
        pct_list = []
        for each_grp in group_id_list_sorted:
            grp_pct = current_marker_num_by_group_dict.get(each_grp, 0)*100/group_to_gnm_num_dict[each_grp]
            grp_pct = float("{0:.2f}".format(grp_pct))
            pct_list.append(grp_pct)
        assess_summary_1_txt_handle.write('%s\t%s\n' % (marker_id, '\t'.join([str(i) for i in pct_list])))

        # write out assess_summary_2_txt
        assess_list = []
        for each_cutoff in present_pct_cutoff_list:

            good_marker = True
            for each_pct in pct_list:
                if each_pct < each_cutoff:
                    good_marker = False

            if each_cutoff not in cutoff_to_qualified_marker_dict:
                cutoff_to_qualified_marker_dict[each_cutoff] = {marker_id}

            if good_marker is True:
                assess_list.append('1')
                cutoff_to_qualified_marker_dict[each_cutoff].add(marker_id)
            else:
                assess_list.append('0')
        assess_summary_2_txt_handle.write('%s\t%s\n' % (marker_id, '\t'.join(assess_list)))

    # write out total in assess_summary_2_txt
    total_stats_list = [str(len(cutoff_to_qualified_marker_dict[each_c])) for each_c in present_pct_cutoff_list]
    assess_summary_2_txt_handle.write('Total\t%s\n' % ('\t'.join(total_stats_list)))
    assess_summary_1_txt_handle.close()
    assess_summary_2_txt_handle.close()

    # copy alignments of qualified marker to corresponding folders
    for each_cutoff in cutoff_to_qualified_marker_dict:
        qualified_marker_set = cutoff_to_qualified_marker_dict[each_cutoff]
        pwd_qualified_marker_dir       = '%s/qualified_marker_PA_%s'                             % (op_dir, each_cutoff)
        pwd_qualified_marker_id_txt    = '%s/qualified_marker_PA_%s_id.txt'                      % (op_dir, each_cutoff)
        pwd_qualified_marker_phy       = '%s/qualified_marker_PA_%s_concatenated.phy'            % (op_dir, each_cutoff)
        pwd_qualified_marker_partition = '%s/qualified_marker_PA_%s_concatenated_partition.txt'  % (op_dir, each_cutoff)

        os.system('mkdir %s' % pwd_qualified_marker_dir)
        for each_marker in qualified_marker_set:
            pwd_marker_aln = '%s/%s.%s' % (trimmed_aln_dir, each_marker, trimmed_aln_ext)
            cp_cmd = 'cp %s %s/' % (pwd_marker_aln, pwd_qualified_marker_dir)
            os.system(cp_cmd)

        # write out id
        with open(pwd_qualified_marker_id_txt, 'w') as pwd_qualified_marker_id_txt_handle:
            pwd_qualified_marker_id_txt_handle.write('%s\n' % '\n'.join(qualified_marker_set))

        # concatenate qualified alignments
        catfasta2phyml_cmd = 'perl %s --sequential --concatenate %s/*.aln > %s 2> %s' % (catfasta2phyml_pl, pwd_qualified_marker_dir, pwd_qualified_marker_phy, pwd_qualified_marker_partition)
        print('running: ' + catfasta2phyml_cmd)
        os.system(catfasta2phyml_cmd)

    # write out iTOL file
    itol_binary_txt_handle = open(itol_binary_txt, 'w')
    itol_binary_txt_handle.write('DATASET_BINARY\n\nSEPARATOR TAB\nDATASET_LABEL\tlabel1\nCOLOR\t#85C1E9\n')
    itol_binary_txt_handle.write('SHOW_LABELS\t1\nLABEL_ROTATION\t45\nLABEL_SHIFT\t5\n')
    itol_binary_txt_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(sorted(marker_id_list)))
    itol_binary_txt_handle.write('FIELD_SHAPES\t%s\n' % '\t'.join(['1']*len(marker_id_list)))
    itol_binary_txt_handle.write('\nDATA\n')
    for each_g in gnm_to_identified_marker_dict:
        g_identified_marker_set = gnm_to_identified_marker_dict[each_g]

        pa_list = []
        for each_m in sorted(marker_id_list):
            if each_m in g_identified_marker_set:
                pa_list.append('1')
            else:
                pa_list.append('-1')
        itol_binary_txt_handle.write('%s\t%s\n' % (each_g, '\t'.join(pa_list)))
    itol_binary_txt_handle.close()

    print('Assessment results exported to:\n%s\n%s' % (assess_summary_1_txt, assess_summary_2_txt))
    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-ta',         required=True,                               help='trimmed alignments')
    parser.add_argument('-tax',        required=True,                               help='extension of trimmed alignments')
    parser.add_argument('-aa',         required=True,                               help='faa file dir')
    parser.add_argument('-aax',        required=True,                               help='faa file ext')
    parser.add_argument('-g',          required=True,                               help='genome group')
    parser.add_argument('-c',          required=False, default='25-50-75-85-100',   help='cutoffs, default: 25-50-75-85-100')
    parser.add_argument('-o',          required=True,                               help='output dir')
    parser.add_argument('-f',          required=False, action="store_true",         help='force overwrite existing output folder')
    parser.add_argument('-pl',         required=True,                               help='path to catfasta2phyml.pl')
    args = vars(parser.parse_args())
    AssessMarkerPA(args)
