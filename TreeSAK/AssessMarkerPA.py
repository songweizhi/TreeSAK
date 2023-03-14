import os
import glob
import argparse
from Bio import SeqIO


AssessMarkerPA_usage = '''
=========================== AssessMarkerPA example commands ===========================

BioSAK AssessMarkerPA -a trimmed_aln -x aln -g gnm_group.txt -c 25-50-75-100 -o s10_assess_marker_PA

Note
1. Extra genomes in gnm_metadata.txt won't affect assessment results.
2. Genomes can not be found in gnm_metadata.txt will trigger an error.
3. Alignments in {trimmed_aln_dir} need to be trimmed before assessment
4. Sequences in MSAs need to be named by genome id.

=======================================================================================
'''


def AssessMarkerPA(args):

    trimmed_aln_dir     = args['a']
    trimmed_aln_ext     = args['x']
    gnm_meta_txt        = args['g']
    cutoff_str          = args['c']
    op_dir              = args['o']
    force_overwriting   = args['f']

    # create folder
    if force_overwriting is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
    else:
        if os.path.isdir(op_dir) is True:
            print('Output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)


    present_pct_cutoff_list = [int(i) for i in cutoff_str.split('-')]
    assess_summary_1_txt    = '%s/assessment_PA.txt'         % op_dir
    assess_summary_2_txt    = '%s/assessment_PA_summary.txt' % op_dir

    trimmed_aln_file_re = '%s/*.%s' % (trimmed_aln_dir, trimmed_aln_ext)
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a',          required=True,                               help='trimmed alignments')
    parser.add_argument('-x',          required=True,                               help='extension of trimmed alignments')
    parser.add_argument('-g',          required=True,                               help='genome group')
    parser.add_argument('-c',          required=False, default='25-50-75-85-100',   help='cutoffs, default: 25-50-75-85-100')
    parser.add_argument('-o',          required=True,                               help='output dir')
    parser.add_argument('-f',          required=False, action="store_true",         help='force overwrite existing output folder')
    args = vars(parser.parse_args())
    AssessMarkerPA(args)


# assess_marker_wd        = '/home-user/wzsong/DateArTree/04_dating_Williams_2017_45_arCOG_assess_marker'
# trimmed_aln_dir         = '/home-user/wzsong/DateArTree/02_identify_marker_gene_Williams_2017_45_arCOG/best_hit_by_marker_5_aln_trimmed'
# gnm_group_txt           = '/home-user/wzsong/DateArTree/01_genome_selection/gnm_metadata.txt'
# deltall_stdout_txt      = '/home-user/wzsong/DateArTree/02_identify_marker_gene_Williams_2017_45_arCOG_DeltaLL/nohup.out'
# present_pct_cutoff_list = [25, 50, 75, 85, 100]
# deltall_keep_pct_list   = [25, 50, 75, 100]
# min_marker_pct_per_gnm  = 75
# min_marker_num          = 20
# force_create_dir        = True
# #catfasta2phyml_pl       = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/catfasta2phyml.pl'
# catfasta2phyml_pl       = '/home-user/wzsong/Scripts/catfasta2phyml.pl'

'''
cd 
BioSAK AssessMarkerPA -a op2/s06_identified_marker_aln_trimmed -x aln -g gnm_group.txt -c 25-50-75-100 -o op2/s10_assess_marker_PA

'''