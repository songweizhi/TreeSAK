import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
from distutils.spawn import find_executable


ExtractMarkerSeq_usage = '''
============================ ExtractMarkerSeq example commands ============================

Dependencies: blastp

BioSAK ExtractMarkerSeq -m marker_ref_seq -mx fa -aa faa_files -aax faa -o op_dir -e "1e-30" -t 6

===========================================================================================
'''


def check_dependencies(program_list):
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not found, program exited!' % ','.join(not_detected_programs))
        exit()


def exe_cmds(cmd_list, num_threads):
    print('Running %s commands with %s cores' % (len(cmd_list), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, cmd_list)
    pool.close()
    pool.join()


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def select_seq(seq_file, id_file,select_option, output_file, one_line, in_fastq):

    # get provided id list
    seq_id_list = set()
    for seq_id in open(id_file):
        seq_id_list.add(seq_id.strip())

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    # extract sequences
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        if select_option == 1:
            if seq_id in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')

        if select_option == 0:
            if seq_id not in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')
    output_file_handle.close()


def ExtractMarkerSeq(args):

    marker_seq_dir          = args['m']
    marker_seq_ext          = args['mx']
    faa_file_dir            = args['aa']
    faa_file_ext            = args['aax']
    op_dir                  = args['o']
    e_value                 = args['e']
    num_of_threads          = args['t']
    force_overwrite         = args['f']

    # check dependencies
    check_dependencies(['blastp'])

    # get marker id set
    marker_seq_re   = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
    marker_seq_list = [os.path.basename(file_name) for file_name in glob.glob(marker_seq_re)]
    marker_id_set = set()
    for each_marker_seq_file in marker_seq_list:
        marker_seq_path, marker_seq_basename, marker_seq_ext = sep_path_basename_ext(each_marker_seq_file)
        marker_id_set.add(marker_seq_basename)

    # get gnm id list
    faa_file_re   = '%s/*.%s' % (faa_file_dir, faa_file_ext)
    faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
    gnm_set = set()
    for each_faa_file in faa_file_list:
        faa_path, faa_basename, faa_ext = sep_path_basename_ext(each_faa_file)
        gnm_set.add(faa_basename)
    gnm_id_list_sorted = sorted([i for i in gnm_set])

    # define output dir
    blastp_cmd_txt                       = '%s/blastp_cmds_%s.txt'                  % (op_dir, (len(gnm_id_list_sorted)*len(marker_id_set)))
    pwd_combined_protein                 = '%s/combined.faa'                        % op_dir
    blast_op_dir                         = '%s/s01_blast_op'                        % op_dir
    best_hit_id_by_marker_dir            = '%s/s02_identified_marker_id'            % op_dir
    best_hit_seq_by_marker_dir           = '%s/s03_identified_marker_seq'           % op_dir
    best_hit_seq_by_marker_dir_renamed   = '%s/s04_identified_marker_seq_renamed'   % op_dir

    # create folder
    if force_overwrite is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
        os.system('mkdir %s' % op_dir)
        os.system('mkdir %s' % blast_op_dir)
    else:
        if os.path.isdir(op_dir) is False:
            os.system('mkdir %s' % op_dir)
        if os.path.isdir(blast_op_dir) is False:
            os.system('mkdir %s' % blast_op_dir)

    os.system('cat %s/*.%s > %s' % (faa_file_dir, faa_file_ext, pwd_combined_protein))

    # get blastp command
    blast_cmd_list = []
    blast_op_to_cmd_dict = dict()
    blastp_cmd_txt_handle = open(blastp_cmd_txt, 'w')
    for gnm_id in gnm_id_list_sorted:
        for each_cog in marker_id_set:
            pwd_blast_op = '%s/%s_vs_%s_blastp.txt'                                                 % (blast_op_dir, gnm_id, each_cog)
            blastp_cmd   = 'blastp -subject %s/%s.fa -evalue %s -outfmt 6 -query %s/%s.faa -out %s' % (marker_seq_dir, each_cog, e_value, faa_file_dir, gnm_id, pwd_blast_op)
            blast_op_to_cmd_dict[pwd_blast_op] = blastp_cmd
            blastp_cmd_txt_handle.write(blastp_cmd + '\n')
            blast_cmd_list.append(blastp_cmd)
    blastp_cmd_txt_handle.close()

    # run blastp
    if force_overwrite is True:
        exe_cmds(blast_cmd_list, num_of_threads)
    else:
        cmds_to_rerun = []
        num_of_good_ones = 0
        for each_blast_op in blast_op_to_cmd_dict:

            look_good = False
            if os.path.isfile(each_blast_op) is True:
                look_good = True
                num_of_good_ones += 1

            if look_good is False:
                cmds_to_rerun.append(blast_op_to_cmd_dict[each_blast_op])

        print('Detected blastp outputs: %s' % num_of_good_ones)
        exe_cmds(cmds_to_rerun, num_of_threads)

    # get best_hit_dict_by_marker
    best_hit_to_gnm_dict = dict()
    best_hit_dict_by_marker = dict()
    for gnm_id in gnm_id_list_sorted:
        for each_cog in marker_id_set:
            current_blastp_op = '%s/%s_vs_%s_blastp.txt' % (blast_op_dir, gnm_id, each_cog)
            # get best hit
            if os.path.isfile(current_blastp_op) is True:
                best_hit_gene = ''
                best_hit_score = 0
                for each_line in open(current_blastp_op):
                    each_line_split = each_line.strip().split('\t')
                    query_id = each_line_split[0]
                    bit_score = float(each_line_split[11])
                    if bit_score > best_hit_score:
                        best_hit_score = bit_score
                        best_hit_gene = query_id

                if best_hit_gene != '':
                    best_hit_to_gnm_dict[best_hit_gene] = gnm_id

                    if each_cog not in best_hit_dict_by_marker:
                        best_hit_dict_by_marker[each_cog] = [best_hit_gene]
                    else:
                        best_hit_dict_by_marker[each_cog].append(best_hit_gene)

    # create output dir
    if os.path.isdir(best_hit_id_by_marker_dir) is False:
        os.system('mkdir %s' % best_hit_id_by_marker_dir)
    if os.path.isdir(best_hit_seq_by_marker_dir) is False:
        os.system('mkdir %s' % best_hit_seq_by_marker_dir)
    if os.path.isdir(best_hit_seq_by_marker_dir_renamed) is False:
        os.system('mkdir %s' % best_hit_seq_by_marker_dir_renamed)

    # write out best hits and extract sequences
    processing_index = 1
    for each_marker in best_hit_dict_by_marker:
        print('Extracting marker sequence %s/%s: %s' % (processing_index, len(best_hit_dict_by_marker), each_marker))
        processing_index += 1

        current_m_hit_list      = best_hit_dict_by_marker[each_marker]
        marker_hits_txt         = ('%s/%s.txt' % (best_hit_id_by_marker_dir,          each_marker)).replace(':', '')
        marker_hits_seq         = ('%s/%s.fa'  % (best_hit_seq_by_marker_dir,         each_marker)).replace(':', '')
        marker_hits_seq_renamed = ('%s/%s.fa'  % (best_hit_seq_by_marker_dir_renamed, each_marker)).replace(':', '')

        with open(marker_hits_txt, 'w') as marker_hits_txt_handle:
            marker_hits_txt_handle.write('\n'.join(current_m_hit_list))

        # extract sequences
        select_seq(pwd_combined_protein, marker_hits_txt, 1, marker_hits_seq, True, False)

        # rename sequences
        marker_hits_seq_renamed_handle = open(marker_hits_seq_renamed, 'w')
        for each_seq in SeqIO.parse(marker_hits_seq, 'fasta'):
            seq_id = each_seq.id
            seq_gnm = best_hit_to_gnm_dict[seq_id]
            marker_hits_seq_renamed_handle.write('>%s\n' % seq_gnm)
            marker_hits_seq_renamed_handle.write('%s\n' % str(each_seq.seq))
        marker_hits_seq_renamed_handle.close()

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',               required=True,                          help='marker seq dir')
    parser.add_argument('-mx',              required=True,                          help='marker seq ext')
    parser.add_argument('-aa',              required=True,                          help='faa file dir')
    parser.add_argument('-aax',             required=True,                          help='faa file ext')
    parser.add_argument('-o',               required=True,                          help='output dir')
    parser.add_argument('-e',               required=True,  default=1e-30,          help='e-value cutoff, default: 1e-30')
    parser.add_argument('-t',               required=True,  type=int,               help='num of threads')
    parser.add_argument('-f',               required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    ExtractMarkerSeq(args)


'''

conda activate mypy3env
cd /home-user/wzsong/DateArTree
python3 MarkerRef2Tree.py -m Marker_set_2_Betts_2018_29_arCOG -mx fa -aa /home-user/wzsong/DateArTree/01_genome_selection_Prokka/d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files -aax faa -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30_demo -e 30 -t 24 -pl /home-user/wzsong/Scripts/catfasta2phyml.pl
submitHPC.sh --cmd "python3 MarkerRef2Tree.py -m Marker_set_2_Betts_2018_29_arCOG -mx fa -aa /home-user/wzsong/DateArTree/01_genome_selection_Prokka/d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files -aax faa -g /home-user/wzsong/DateArTree/gnm_group.txt -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30_demo -e 30 -t 24 -pl /home-user/wzsong/Scripts/catfasta2phyml.pl" -n 24 -c Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30_demo

cd /home-user/wzsong/DateArTree
python3 MarkerRef2Tree.py -m Marker_set_2_Betts_2018_29_arCOG -mx fa -aa /home-user/wzsong/DateArTree/01_genome_selection_Prokka/d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files -aax faa -g /home-user/wzsong/DateArTree/gnm_group.txt -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30_demo -e 30 -t 12 -pl /home-user/wzsong/Scripts/catfasta2phyml.pl -g gnm_group.txt -skip_align_trim -jst 6 -qsub

cd /Users/songweizhi/Desktop/demo
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/MarkerRef2Tree.py -m Marker_set_2_Betts_2018_29_arCOG -mx fa -aa d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files -aax faa -g gnm_group.txt -o Marker_set_2_Betts_2018_29_arCOG_Marker2Tree_e30_demo -e 30 -t 10 -pl /Users/songweizhi/Scripts/catfasta2phyml.pl -g gnm_group.txt -skip_align_trim -jst 6 -qsub

'''
