import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
from distutils.spawn import find_executable


Marker2Tree_usage = '''
============================= Marker2Tree example commands =============================

Dependencies: blastp, mafft-einsi, trimal, iqtree

# example commands
BioSAK Marker2Tree -m marker_seq -mx fa -aa gnm_faa_files -aax faa -o op_dir -e 30 -t 6 -pl /home-user/wzsong/Scripts/catfasta2phyml.pl

========================================================================================
'''


def check_dependencies(program_list):
    # check whether executables exist
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


def Marker2Tree(args):

    marker_seq_dir    = args['m']
    marker_seq_ext    = args['mx']
    faa_file_dir      = args['aa']
    faa_file_ext      = args['aax']
    op_dir            = args['o']
    e_value           = args['e']
    num_of_threads    = args['t']
    force_overwrite   = args['f']
    catfasta2phyml_pl = args['pl']

    # check_dependencies
    check_dependencies(['blastp', 'mafft-einsi', 'trimal', 'iqtree'])

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
    blastp_cmd_txt                       = '%s/blastp_cmds_%s.txt'                                  % (op_dir, (len(gnm_id_list_sorted)*len(marker_id_set)))
    pwd_combined_protein                 = '%s/combined.faa'                                        % op_dir
    blast_op_dir                         = '%s/s01_blast_op'                                        % op_dir
    best_hit_id_by_marker_dir            = '%s/s02_identified_marker_id'                            % op_dir
    best_hit_seq_by_marker_dir           = '%s/s03_identified_marker_seq'                           % op_dir
    best_hit_seq_by_marker_dir_renamed   = '%s/s04_identified_marker_seq_renamed'                   % op_dir
    best_hit_aln_by_marker_dir           = '%s/s05_identified_marker_aln'                           % op_dir
    best_hit_aln_by_marker_dir_trimmed   = '%s/s06_identified_marker_aln_trimmed'                   % op_dir
    best_hit_aln_by_marker_dir_trimmed_c = '%s/s07_identified_marker_aln_trimmed_concatenated'      % op_dir
    get_guide_tree_dir                   = '%s/s08_get_guide_tree'                                  % op_dir
    iqtree_dir                           = '%s/s09_iqtree_C60_PMSF'                                 % op_dir
    pwd_concatenate_phy                  = '%s/concatenated.phy'                                    % best_hit_aln_by_marker_dir_trimmed_c
    pwd_concatenate_phy_partition        = '%s/concatenated_partitions.txt'                         % best_hit_aln_by_marker_dir_trimmed_c
    pwd_guide_tree                       = '%s/guide_tree.treefile'                                 % get_guide_tree_dir

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
            #blastp_cmd = 'blastp -subject /home-user/wzsong/DateArTree/02_Williams_2017_45_arCOG/%s.fa -evalue 1e-30 -outfmt 6 -query /home-user/wzsong/DateArTree/01_genome_selection_Prokka/d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files/%s.faa -out %s/%s_vs_%s_blastp.txt' % (each_cog, gnm_id, blast_op, gnm_id, each_cog)
            pwd_blast_op = '%s/%s_vs_%s_blastp.txt'                                                     % (blast_op_dir, gnm_id, each_cog)
            blastp_cmd   = 'blastp -subject %s/%s.fa -evalue 1e-%s -outfmt 6 -query %s/%s.faa -out %s'  % (marker_seq_dir, each_cog, e_value, faa_file_dir, gnm_id, pwd_blast_op)
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
    if force_overwrite is True:
        os.system('mkdir %s' % best_hit_id_by_marker_dir)
        os.system('mkdir %s' % best_hit_seq_by_marker_dir)
        os.system('mkdir %s' % best_hit_seq_by_marker_dir_renamed)
        os.system('mkdir %s' % best_hit_aln_by_marker_dir)
        os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed)
    else:
        if os.path.isdir(best_hit_id_by_marker_dir) is False:
            os.system('mkdir %s' % best_hit_id_by_marker_dir)
        if os.path.isdir(best_hit_seq_by_marker_dir) is False:
            os.system('mkdir %s' % best_hit_seq_by_marker_dir)
        if os.path.isdir(best_hit_seq_by_marker_dir_renamed) is False:
            os.system('mkdir %s' % best_hit_seq_by_marker_dir_renamed)
        if os.path.isdir(best_hit_aln_by_marker_dir) is False:
            os.system('mkdir %s' % best_hit_aln_by_marker_dir)
        if os.path.isdir(best_hit_aln_by_marker_dir_trimmed) is False:
            os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed)

    # write out best hits and extract their sequences
    processing_index = 1
    for each_marker in best_hit_dict_by_marker:

        print('Processing %s/%s: %s' % (processing_index, len(best_hit_dict_by_marker), each_marker))
        processing_index += 1

        current_m_hit_list      = best_hit_dict_by_marker[each_marker]
        marker_hits_txt         = ('%s/%s.txt' % (best_hit_id_by_marker_dir,          each_marker)).replace(':', '')
        marker_hits_seq         = ('%s/%s.fa'  % (best_hit_seq_by_marker_dir,         each_marker)).replace(':', '')
        marker_hits_seq_renamed = ('%s/%s.fa'  % (best_hit_seq_by_marker_dir_renamed, each_marker)).replace(':', '')
        marker_hits_aln         = ('%s/%s.aln' % (best_hit_aln_by_marker_dir,         each_marker)).replace(':', '')
        marker_hits_aln_trimmed = ('%s/%s.aln' % (best_hit_aln_by_marker_dir_trimmed, each_marker)).replace(':', '')

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

        # run mafft-einsi
        mafft_cmd = 'mafft-einsi --quiet %s > %s' % (marker_hits_seq_renamed, marker_hits_aln)
        #print('running: ' + mafft_cmd)
        os.system(mafft_cmd)

        # trim msa
        trimal_cmd = 'trimal -in %s -out %s -automated1' % (marker_hits_aln, marker_hits_aln_trimmed)
        #print('running: ' + trimal_cmd)
        os.system(trimal_cmd)

    # concatenate trimmed msa
    os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed_c)
    catfasta2phyml_cmd = 'perl %s --sequential --concatenate %s/*.aln > %s 2> %s' % (catfasta2phyml_pl, best_hit_aln_by_marker_dir_trimmed, pwd_concatenate_phy, pwd_concatenate_phy_partition)
    print('running: ' + catfasta2phyml_cmd)
    os.system(catfasta2phyml_cmd)

    # get guide tree
    os.system('mkdir %s' % get_guide_tree_dir)
    get_guide_tree_cmd = 'iqtree -s %s --prefix %s/guide_tree --seqtype AA -m LG -T 10 -B 1000 --alrt 1000 --quiet' % (pwd_concatenate_phy, get_guide_tree_dir)
    print('running: ' + get_guide_tree_cmd)
    os.system(get_guide_tree_cmd)

    # run C60 + PMSF
    os.system('mkdir %s' % iqtree_dir)
    iqtree_c60_pmsf_cmd = 'iqtree -s %s --prefix %s/concatenated --seqtype AA -m LG+G+F+C60 -T 10 -B 1000 --alrt 1000 --quiet -ft %s' % (pwd_concatenate_phy, iqtree_dir, pwd_guide_tree)
    print('running: ' + iqtree_c60_pmsf_cmd)
    os.system(iqtree_c60_pmsf_cmd)

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',           required=True,                          help='marker seq dir')
    parser.add_argument('-mx',          required=True,                          help='marker seq ext')
    parser.add_argument('-aa',          required=True,                          help='faa file dir')
    parser.add_argument('-aax',         required=True,                          help='faa file ext')
    parser.add_argument('-o',           required=True,                          help='output dir')
    parser.add_argument('-e',           required=True,                          help='e-value')
    parser.add_argument('-t',           required=True,  type=int,               help='num of threads')
    parser.add_argument('-f',           required=False, action="store_true",    help='force overwrite')
    parser.add_argument('-pl',          required=True,                          help='path to catfasta2phyml.pl')
    args = vars(parser.parse_args())
    Marker2Tree(args)
