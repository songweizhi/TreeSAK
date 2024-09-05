import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO
import multiprocessing as mp
from distutils.spawn import find_executable


MarkerRef2Tree_usage = '''
============================= MarkerRef2Tree example commands =============================

Dependencies: java, blastp, mafft-einsi, trimal, iqtree2

# example commands
TreeSAK MarkerRef2Tree -i faa_files -x faa -m marker_seq -mx fa -o output_dir -bmge -e 10 -t 6 -c 85
TreeSAK MarkerRef2Tree -i faa_files -x faa -m marker_seq -mx fa -o output_dir -bmge -e 10 -t 6 -c 75,100

# file extension need to be faa

===========================================================================================
'''


def check_dependencies(program_list):

    # check whether executables exist
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    # report
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


def catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, concatenated_msa_fasta, partition_file):

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

    # convert fasta to phy
    AlignIO.convert(concatenated_msa_fasta, 'fasta', concatenated_msa_phy, 'phylip-relaxed')


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


def AssessMarkerPA(trimmed_aln_dir, gnm_set, group_to_gnm_dict, present_pct_cutoff_list, op_dir):

    group_to_gnm_num_dict = dict()
    gnm_to_group_dict     = dict()
    for each_g in group_to_gnm_dict:
        gnm_member_list = group_to_gnm_dict[each_g]
        group_to_gnm_num_dict[each_g] = len(gnm_member_list)
        for each_gnm in gnm_member_list:
            gnm_to_group_dict[each_gnm] = each_g
    group_id_list_sorted = sorted(list(group_to_gnm_dict.keys()))

    # exit program if group information is missing
    gnms_without_group_info = set()
    for gnm in gnm_set:
        if gnm not in gnm_to_group_dict:
            gnms_without_group_info.add(gnm)

    if len(gnms_without_group_info) > 0:
        print('Group information for the following genomes are missing, program exited!')
        print(','.join(gnms_without_group_info))
        print('Group information for the above genomes are missing, program exited!')
        exit()

    # read in provided cutoffs
    assess_summary_1_txt    = '%s/assess_by_PA.txt'                % op_dir
    assess_summary_2_txt    = '%s/assess_by_PA_summary.txt'        % op_dir
    itol_binary_txt         = '%s/assess_by_PA_iTOL_binary.txt'    % op_dir

    trimmed_aln_file_re = '%s/*.aln' % (trimmed_aln_dir)
    trimmed_aln_file_list = [os.path.basename(file_name) for file_name in glob.glob(trimmed_aln_file_re)]

    assess_summary_1_txt_handle = open(assess_summary_1_txt, 'w')
    assess_summary_1_txt_handle.write('Marker\t%s\n' % '\t'.join([str(i) for i in group_id_list_sorted]))
    assess_summary_2_txt_handle = open(assess_summary_2_txt, 'w')
    assess_summary_2_txt_handle.write('Marker\t%s\n' % '\t'.join([str(i) for i in present_pct_cutoff_list]))
    cutoff_to_qualified_marker_dict = dict()
    gnm_to_identified_marker_dict = dict()
    marker_id_list = []
    for each_aln in trimmed_aln_file_list:

        marker_id = each_aln.split(('.aln'))[0]
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
                if gnm_group not in current_marker_num_by_group_dict:
                    current_marker_num_by_group_dict[gnm_group] = 1
                else:
                    current_marker_num_by_group_dict[gnm_group] += 1
            else:
                print('Not all genomes used to generate the MSA being found in -aa, program exited!')
                exit()

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
        qualified_marker_set        = cutoff_to_qualified_marker_dict[each_cutoff]
        pwd_qualified_marker_dir    = '%s/marker_PA%s'        % (op_dir, each_cutoff)
        pwd_qualified_marker_id_txt = '%s/marker_PA%s_id.txt' % (op_dir, each_cutoff)

        os.system('mkdir %s' % pwd_qualified_marker_dir)
        for each_marker in qualified_marker_set:
            pwd_marker_aln = '%s/%s.aln' % (trimmed_aln_dir, each_marker)
            cp_cmd = 'cp %s %s/' % (pwd_marker_aln, pwd_qualified_marker_dir)
            os.system(cp_cmd)

        # write out id
        with open(pwd_qualified_marker_id_txt, 'w') as pwd_qualified_marker_id_txt_handle:
            pwd_qualified_marker_id_txt_handle.write('%s\n' % '\n'.join(qualified_marker_set))

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


def get_gap_stats(msa_in_fa, stats_txt):

    gap_pct_dict = dict()
    for each_seq in SeqIO.parse(msa_in_fa, 'fasta'):
        seq_id = each_seq.id
        seq_str = str(each_seq.seq)
        gap_pct = seq_str.count('-')*100/len(seq_str)
        gap_pct = float("{0:.2f}".format(gap_pct))
        gap_pct_dict[seq_id] = gap_pct

    gap_pct_sorted = sorted(gap_pct_dict.items(), key=lambda x:x[1])

    stats_txt_handle = open(stats_txt, 'w')
    stats_txt_handle.write('Sequence\tGap\n')
    for each_seq in gap_pct_sorted:
        stats_txt_handle.write('%s\t%s\n' % (each_seq[0], each_seq[1]))
    stats_txt_handle.close()


def BMGE(msa_in, op_prefix, trim_model, entropy_score_cutoff):

    # define file name
    msa_out_phylip = '%s.BMGE.phylip' % op_prefix
    msa_out_fasta  = '%s.BMGE.fasta'  % op_prefix
    msa_out_nexus  = '%s.BMGE.nexus'  % op_prefix
    msa_out_html   = '%s.BMGE.html'   % op_prefix

    # specify path to BMGE.jar
    current_file_path   = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar        = '%s/BMGE.jar' % current_file_path

    # run BMGE
    bmge_cmd = 'java -jar %s -i %s -m %s -t AA -h %s -op %s -of %s -on %s -oh %s' % (pwd_bmge_jar, msa_in, trim_model, entropy_score_cutoff, msa_out_phylip, msa_out_fasta, msa_out_nexus, msa_out_html)
    print('Running %s' % bmge_cmd)
    os.system(bmge_cmd)


def MarkerRef2Tree(args):

    faa_file_dir                = args['i']
    faa_file_ext                = args['x']
    marker_seq_dir              = args['m']
    marker_seq_ext              = args['mx']
    gnm_group_txt               = args['g']
    op_dir                      = args['o']
    e_value                     = args['e']
    num_of_threads              = args['t']
    force_overwrite             = args['f']
    pa_cutoff_str               = args['c']
    minimal_marker_number       = args['mmn']
    run_psiblast                = args['psiblast']
    trim_with_bmge              = args['bmge']
    bmge_trim_model             = args['bmge_m']
    bmge_entropy_score_cutoff   = args['bmge_esc']

    # check dependencies
    check_dependencies(['java', 'psiblast', 'blastp', 'mafft-einsi', 'trimal', 'iqtree'])

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar      = '%s/BMGE.jar' % current_file_path

    # get marker id set
    marker_seq_re   = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
    marker_seq_list = [os.path.basename(file_name) for file_name in glob.glob(marker_seq_re)]
    marker_id_set = set()
    for each_marker_seq_file in marker_seq_list:
        _, marker_seq_basename, _ = sep_path_basename_ext(each_marker_seq_file)
        marker_id_set.add(marker_seq_basename)

    # get gnm id list
    faa_file_re   = '%s/*.%s' % (faa_file_dir, faa_file_ext)
    faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
    gnm_set = set()
    for each_faa_file in faa_file_list:
        faa_path, faa_basename, faa_ext = sep_path_basename_ext(each_faa_file)
        gnm_set.add(faa_basename)

    gnm_id_list_sorted = sorted([i for i in gnm_set])

    #################### check genome grouping files ####################

    group_to_gnm_dict = dict()
    if gnm_group_txt is None:
        group_to_gnm_dict['group_1'] = set()
        for each_gnm in gnm_id_list_sorted:
            group_to_gnm_dict['group_1'].add(each_gnm)
    else:
        if os.path.isfile(gnm_group_txt) is False:
            print('Specified %s not found, program exited!' % gnm_group_txt)
            exit()
        else:
            for each_gnm in open(gnm_group_txt):
                each_gnm_split = each_gnm.strip().split('\t')
                gnm_id = each_gnm_split[0]
                domain_name = each_gnm_split[1]
                if gnm_id in gnm_set:
                    if domain_name not in group_to_gnm_dict:
                        group_to_gnm_dict[domain_name] = {gnm_id}
                    else:
                        group_to_gnm_dict[domain_name].add(gnm_id)

    ############################################## check file/folder name ##############################################

    blastp_cmd_txt                          = '%s/s01_blast_cmds_%s.txt'            % (op_dir, (len(marker_seq_list)*len(faa_file_list)))
    pwd_combined_protein                    = '%s/combined.faa'                     % op_dir
    blast_op_dir                            = '%s/s01_blast'                        % op_dir
    best_hit_id_by_marker_dir               = '%s/s02_marker_id'                    % op_dir
    best_hit_seq_by_marker_dir              = '%s/s03_marker_seq'                   % op_dir
    best_hit_seq_by_marker_dir_renamed      = '%s/s04_marker_seq_by_genome_name'    % op_dir
    best_hit_aln_by_marker_dir              = '%s/s05_marker_aln'                   % op_dir
    best_hit_aln_by_marker_dir_trimmed      = '%s/s06_marker_aln_trimal'            % op_dir
    if trim_with_bmge is True:
        best_hit_aln_by_marker_dir_trimmed  = '%s/s06_marker_aln_BMGE'              % op_dir
    assess_marker_pa_dir                    = '%s/s07_assess_marker_PA'             % op_dir
    trimmed_msa_PA_concatenated_dir         = '%s/s08_iqtree_wd'                    % op_dir
    iqtree_dir                              = '%s/s08_iqtree_wd'                    % op_dir
    iqtree_cmd_txt                          = '%s/iqtree_cmds.txt'                  % iqtree_dir

    ####################################################################################################################

    # create folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % blast_op_dir)

    # get blastp command
    blast_cmd_list = []
    blast_op_to_cmd_dict = dict()
    blastp_cmd_txt_handle = open(blastp_cmd_txt, 'w')
    for gnm_id in gnm_id_list_sorted:
        for each_cog in marker_id_set:
            pwd_blast_op = '%s/%s_vs_%s_blastp.txt'                                                  % (blast_op_dir, gnm_id, each_cog)

            blast_cmd   = 'blastp -subject %s/%s.%s -evalue %s -outfmt 6 -query %s/%s.%s -out %s' % (marker_seq_dir, each_cog, marker_seq_ext, e_value, faa_file_dir, gnm_id, faa_file_ext, pwd_blast_op)
            if run_psiblast is True:
                blast_cmd = 'psiblast -subject %s/%s.%s -evalue %s -outfmt 6 -query %s/%s.%s -out %s' % (marker_seq_dir, each_cog, marker_seq_ext, e_value, faa_file_dir, gnm_id, faa_file_ext, pwd_blast_op)

            blast_op_to_cmd_dict[pwd_blast_op] = blast_cmd
            blastp_cmd_txt_handle.write(blast_cmd + '\n')
            blast_cmd_list.append(blast_cmd)
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

    # create dir
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

    os.system('cat %s/*.%s > %s' % (faa_file_dir, faa_file_ext, pwd_combined_protein))

    processing_index = 1
    for each_marker in best_hit_dict_by_marker:
        print('Processing (extract sequence, align and trim) marker %s/%s: %s' % (processing_index, len(best_hit_dict_by_marker), each_marker))
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
        mafft_cmd = 'mafft-einsi --thread %s --quiet %s  > %s' % (num_of_threads, marker_hits_seq_renamed, marker_hits_aln)
        os.system(mafft_cmd)

        # trim msa
        trim_cmd = 'trimal -in %s -out %s -automated1' % (marker_hits_aln, marker_hits_aln_trimmed)
        if trim_with_bmge is True:
            trim_cmd = 'java -jar %s -i %s -m %s -t AA -h %s -of %s' % (pwd_bmge_jar, marker_hits_aln, bmge_trim_model, bmge_entropy_score_cutoff, marker_hits_aln_trimmed)
        os.system(trim_cmd)

    ########## Assess marker by PA ##########

    present_pct_cutoff_list = [int(i) for i in pa_cutoff_str.split(',')]

    # create output dir
    if os.path.isdir(assess_marker_pa_dir) is True:
        os.system('rm -r %s' % assess_marker_pa_dir)
    os.system('mkdir %s' % assess_marker_pa_dir)

    AssessMarkerPA(best_hit_aln_by_marker_dir_trimmed, gnm_set, group_to_gnm_dict, present_pct_cutoff_list, assess_marker_pa_dir)

    ########## concatenate marker ##########

    # create output dir
    if os.path.isdir(trimmed_msa_PA_concatenated_dir) is True:
        os.system('rm -r %s' % trimmed_msa_PA_concatenated_dir)
    os.system('mkdir %s' % trimmed_msa_PA_concatenated_dir)

    qualified_cutoff_list = []
    for each_c in present_pct_cutoff_list:
        current_cutoff_trimmed_msa_dir = '%s/marker_PA%s' % (assess_marker_pa_dir, each_c)
        trimmed_msa_re   = '%s/*.aln' % current_cutoff_trimmed_msa_dir
        trimmed_msa_list = [os.path.basename(file_name) for file_name in glob.glob(trimmed_msa_re)]

        if len(trimmed_msa_list) < minimal_marker_number:
            print('The number of qualified marker under PA cutoff %s: %s, skipped!' % (each_c, len(trimmed_msa_list)))
        else:
            qualified_cutoff_list.append(each_c)
            pwd_concatenated_marker_phy       = '%s/marker_pa%s.phy'           % (trimmed_msa_PA_concatenated_dir, each_c)
            pwd_concatenated_marker_fasta     = '%s/marker_pa%s.fasta'         % (trimmed_msa_PA_concatenated_dir, each_c)
            pwd_concatenated_marker_partition = '%s/marker_pa%s_partition.txt' % (trimmed_msa_PA_concatenated_dir, each_c)
            catfasta2phy(current_cutoff_trimmed_msa_dir, 'aln', pwd_concatenated_marker_phy, pwd_concatenated_marker_fasta, pwd_concatenated_marker_partition)

    ########## get guide tree and C60+PMSF tree for each set of marker set ##########

    iqtree_cmd_txt_handle = open(iqtree_cmd_txt, 'w')
    iqtree_cmd_list = []
    for each_c in qualified_cutoff_list:
        os.system('mkdir %s/PA%s_guide_tree'    % ((iqtree_dir, each_c)))
        os.system('mkdir %s/PA%s_PMSF_C60_tree' % ((iqtree_dir, each_c)))

        msa_to_use = 'marker_pa%s.fasta' % each_c

        get_guide_tree_cmd = 'iqtree2 --seqtype AA -B 1000 --alrt 1000 --quiet -T %s -s %s --prefix PA%s_guide_tree/PA%s_guide_tree -m LG'                                                       % (num_of_threads, msa_to_use, each_c, each_c)
        get_c60_tree_cmd   = 'iqtree2 --seqtype AA -B 1000 --alrt 1000 --quiet -T %s -s %s --prefix PA%s_PMSF_C60_tree/PA%s_PMSF_C60 -m LG+C60+F+G -ft PA%s_guide_tree/PA%s_guide_tree.treefile' % (num_of_threads, msa_to_use, each_c, each_c, each_c, each_c)
        cmds_in_one_line   = '%s; %s' % (get_guide_tree_cmd, get_c60_tree_cmd)
        iqtree_cmd_txt_handle.write('%s\n' % cmds_in_one_line)
        iqtree_cmd_list.append(cmds_in_one_line)
    iqtree_cmd_txt_handle.close()

    # run iqtree
    os.chdir(iqtree_dir)
    for each_cmd in iqtree_cmd_list:
        print('Running: %s' % each_cmd)
        os.system(each_cmd)

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',               required=True,                              help='faa dir')
    parser.add_argument('-x',               required=False, default='faa',              help='faa file extension, default: faa')
    parser.add_argument('-m',               required=True,                              help='marker seq dir, file extension need to be faa')
    parser.add_argument('-mx',              required=False, default='faa',              help='marker seq file extension, default: faa')
    parser.add_argument('-g',               required=False, default=None,               help='genome group')
    parser.add_argument('-c',               required=False, default='85',               help='presence-absence cutoffs, default: 85')
    parser.add_argument('-o',               required=True,                              help='output dir')
    parser.add_argument('-e',               required=False, default='1e-30',            help='e-value cutoff, default: 1e-30')
    parser.add_argument('-t',               required=True,  type=int,                   help='num of threads')
    parser.add_argument('-mmn',             required=False, default=1, type=int,        help='minimal marker number, default: 1')
    parser.add_argument('-psiblast',        required=False, action="store_true",        help='run psiblast')
    parser.add_argument('-bmge',            required=False, action="store_true",        help='trim MSA with BMGE, default is trimal')
    parser.add_argument('-bmge_m',          required=False, default='BLOSUM30',         help='BMGE trim model, default: BLOSUM30')
    parser.add_argument('-bmge_esc',        required=False, default='0.55',             help='BMGE entropy score cutoff, default: 0.55')
    parser.add_argument('-f',               required=False, action="store_true",        help='force overwrite')
    args = vars(parser.parse_args())
    MarkerRef2Tree(args)
