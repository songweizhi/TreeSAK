import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO
import multiprocessing as mp
from distutils.spawn import find_executable


MarkerRef2Tree_usage = '''
============================= MarkerRef2Tree example commands =============================

Dependencies: blastp, mafft-einsi, trimal, iqtree

# example commands
TreeSAK MarkerRef2Tree -m marker_seq -mx fa -aa gnm_faa_files -aax faa -o op_dir -e 30 -t 6

===========================================================================================
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


def catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, partition_file):

    concatenated_msa_fasta = '%s.fasta' % concatenated_msa_phy
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

    # convert msa in fasta to phy
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


def AssessMarkerPA(trimmed_aln_dir, gnm_set, gnm_group_txt, present_pct_cutoff_list, op_dir):

    # read in genome metadata
    group_to_gnm_dict     = dict()
    group_to_gnm_num_dict = dict()
    gnm_to_group_dict     = dict()
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

    # read in provided cutoffs
    assess_summary_1_txt    = '%s/assessment_PA.txt'                % op_dir
    assess_summary_2_txt    = '%s/assessment_PA_summary.txt'        % op_dir
    itol_binary_txt         = '%s/assessment_PA_iTOL_binary.txt'    % op_dir

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
        pwd_qualified_marker_dir    = '%s/qualified_marker_PA_%s'        % (op_dir, each_cutoff)
        pwd_qualified_marker_id_txt = '%s/qualified_marker_PA_%s_id.txt' % (op_dir, each_cutoff)

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


def MarkerRef2Tree(args):

    marker_seq_dir          = args['m']
    marker_seq_ext          = args['mx']
    faa_file_dir            = args['aa']
    faa_file_ext            = args['aax']
    gnm_group_txt           = args['g']
    op_dir                  = args['o']
    e_value                 = args['e']
    num_of_threads          = args['t']
    force_overwrite         = args['f']
    js_cpu_num              = args['jst']
    pa_cutoff_str           = args['pac']
    skip_align_trim         = args['skip_align_trim']
    submit_job              = args['qsub']
    minimal_marker_number   = args['mmn']
    print(marker_seq_ext)
    # check dependencies
    check_dependencies(['blastp', 'mafft-einsi', 'trimal', 'iqtree'])

    # check input files
    if os.path.isfile(gnm_group_txt) is False:
        print('%s not found, program exited!' % gnm_group_txt)
        exit()

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

    # define output dir
    blastp_cmd_txt                       = '%s/blastp_cmds_%s.txt'                                  % (op_dir, (len(gnm_id_list_sorted)*len(marker_id_set)))
    pwd_combined_protein                 = '%s/combined.faa'                                        % op_dir
    blast_op_dir                         = '%s/s01_blast_op'                                        % op_dir
    best_hit_id_by_marker_dir            = '%s/s02_identified_marker_id'                            % op_dir
    best_hit_seq_by_marker_dir           = '%s/s03_identified_marker_seq'                           % op_dir
    best_hit_seq_by_marker_dir_renamed   = '%s/s04_identified_marker_seq_renamed'                   % op_dir
    best_hit_aln_by_marker_dir           = '%s/s05_identified_marker_aln'                           % op_dir
    best_hit_aln_by_marker_dir_trimmed   = '%s/s06_identified_marker_aln_trimmed'                   % op_dir
    assess_marker_pa_dir                 = '%s/s07_assess_marker_PA'                                % op_dir
    trimmed_msa_PA_concatenated_dir      = '%s/s08_marker_after_PA_concatenated'                    % op_dir
    iqtree_dir                           = '%s/s09_iqtree_for_deltaLL'                              % op_dir
    deltall_dir                          = '%s/s10_assess_marker_deltaLL'                           % op_dir

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
            pwd_blast_op = '%s/%s_vs_%s_blastp.txt'                                                 % (blast_op_dir, gnm_id, each_cog)
            blastp_cmd   = 'blastp -subject %s/%s.%s -evalue %s -outfmt 6 -query %s/%s.faa -out %s' % (marker_seq_dir, each_cog, marker_seq_ext, e_value, faa_file_dir, gnm_id, pwd_blast_op)
            print(marker_seq_ext)
            print(blastp_cmd)
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
    if os.path.isdir(best_hit_aln_by_marker_dir) is False:
        os.system('mkdir %s' % best_hit_aln_by_marker_dir)
    if os.path.isdir(best_hit_aln_by_marker_dir_trimmed) is False:
        os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed)

    # write out best hits and extract sequences
    if skip_align_trim is True:
        print('Skipping the extraction, alignment and trimming of markers')
    else:
        processing_index = 1
        for each_marker in best_hit_dict_by_marker:
            print('Processing (extract, align and trim) marker %s/%s: %s' % (processing_index, len(best_hit_dict_by_marker), each_marker))
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
            #print('running: ' + mafft_cmd)
            os.system(mafft_cmd)

            # trim msa
            trimal_cmd = 'trimal -in %s -out %s -automated1' % (marker_hits_aln, marker_hits_aln_trimmed)
            #print('running: ' + trimal_cmd)
            os.system(trimal_cmd)

    ########## Assess marker by PA ##########

    present_pct_cutoff_list = [int(i) for i in pa_cutoff_str.split('-')]

    # create output dir
    if os.path.isdir(assess_marker_pa_dir) is True:
        os.system('rm -r %s' % assess_marker_pa_dir)
    os.system('mkdir %s' % assess_marker_pa_dir)

    # Assess marker by PA
    AssessMarkerPA(best_hit_aln_by_marker_dir_trimmed, gnm_set, gnm_group_txt, present_pct_cutoff_list, assess_marker_pa_dir)

    ########## concatenate marker set ##########

    # create output dir
    if os.path.isdir(trimmed_msa_PA_concatenated_dir) is True:
        os.system('rm -r %s' % trimmed_msa_PA_concatenated_dir)
    os.system('mkdir %s' % trimmed_msa_PA_concatenated_dir)

    qualified_cutoff_list = []
    for each_cutoff in present_pct_cutoff_list:
        current_cutoff_trimmed_msa_dir = '%s/qualified_marker_PA_%s' % (assess_marker_pa_dir, each_cutoff)
        trimmed_msa_re   = '%s/*.aln' % current_cutoff_trimmed_msa_dir
        trimmed_msa_list = [os.path.basename(file_name) for file_name in glob.glob(trimmed_msa_re)]

        if len(trimmed_msa_list) < minimal_marker_number:
            print('The number of qualified marker under PA cutoff %s: %s, skipped!' % (each_cutoff, len(trimmed_msa_list)))
        else:
            qualified_cutoff_list.append(each_cutoff)
            pwd_concatenated_marker_phy       = '%s/qualified_marker_PA_%s_concatenated.phy'                % (trimmed_msa_PA_concatenated_dir, each_cutoff)
            pwd_concatenated_marker_partition = '%s/qualified_marker_PA_%s_concatenated_partition.txt'      % (trimmed_msa_PA_concatenated_dir, each_cutoff)
            #catfasta2phyml_cmd                = 'perl %s --sequential --concatenate %s/*.aln > %s 2> %s'    % (catfasta2phyml_pl, current_cutoff_trimmed_msa_dir, pwd_concatenated_marker_phy, pwd_concatenated_marker_partition)
            #print('running: ' + catfasta2phyml_cmd)
            #os.system(catfasta2phyml_cmd)
            catfasta2phy(current_cutoff_trimmed_msa_dir, 'aln', pwd_concatenated_marker_phy, pwd_concatenated_marker_partition)

    ########## get guide tree and C60+PMSF tree for each set of marker set ##########

    # create output dir
    if os.path.isdir(iqtree_dir) is True:
        os.system('rm -r %s' % iqtree_dir)
    os.system('mkdir %s' % iqtree_dir)

    current_dir = os.getcwd()
    for each_qualified_cutoff in qualified_cutoff_list:

        # create output dir
        get_guide_tree_wd    = '%s/PA_%s_guide_tree'    % (iqtree_dir, each_qualified_cutoff)
        get_c60_pmsf_tree_wd = '%s/PA_%s_C60_PMSF_tree' % (iqtree_dir, each_qualified_cutoff)
        os.system('mkdir %s' % get_guide_tree_wd)
        os.system('mkdir %s' % get_c60_pmsf_tree_wd)

        # define file name
        pwd_concatenated_marker_phy = '%s/%s/qualified_marker_PA_%s_concatenated.phy'                                                                   % (os.getcwd(), trimmed_msa_PA_concatenated_dir, each_qualified_cutoff)
        pwd_guide_tree              = '%s/%s/guide_tree.treefile'                                                                                       % (os.getcwd(), get_guide_tree_wd)
        get_guide_tree_cmd          = 'iqtree -s %s --prefix %s/%s/guide_tree --seqtype AA -m LG -T %s -B 1000 --alrt 1000 --quiet'                     % (pwd_concatenated_marker_phy, os.getcwd(), get_guide_tree_wd, js_cpu_num)
        get_c60_tree_cmd            = 'iqtree -s %s --prefix %s/%s/concatenated --seqtype AA -m LG+G+F+C60 -T %s -B 1000 --alrt 1000 --quiet -ft %s'    % (pwd_concatenated_marker_phy, os.getcwd(), get_c60_pmsf_tree_wd, js_cpu_num, pwd_guide_tree)

        # write out job script
        js_iqtree = '%s/js_iqtree_%s.sh' % (iqtree_dir, each_qualified_cutoff)
        with open(js_iqtree, 'w') as js_iqtree_handle:
            js_iqtree_handle.write('#!/bin/bash\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task %s\n' % js_cpu_num)
            js_iqtree_handle.write('cd %s/%s\n'    % (os.getcwd(), get_guide_tree_wd))
            js_iqtree_handle.write('%s\n' % get_guide_tree_cmd)
            js_iqtree_handle.write('cd %s/%s\n'    % (os.getcwd(), get_c60_pmsf_tree_wd))
            js_iqtree_handle.write('%s\n' % get_c60_tree_cmd)

        # submit job script
        print('Commands for running iqtree exported to: %s' % js_iqtree)
        if submit_job is True:
            os.chdir(iqtree_dir)
            os.system('qsub js_iqtree_%s.sh' % each_qualified_cutoff)
            os.chdir(current_dir)

    ########## provide commands for running DeltaLL ##########

    # create output dir
    if os.path.isdir(deltall_dir) is True:
        os.system('rm -r %s' % deltall_dir)
    os.system('mkdir %s' % deltall_dir)

    print('Suggested commands for running DeltaLL')
    for each_qualified_cutoff in qualified_cutoff_list:
        pwd_js          = '%s/%s/js_deltall_PA%s.sh'                        % (os.getcwd(), deltall_dir, each_qualified_cutoff)
        deltaLL_wd      = '%s/%s/PA_%s'                                     % (os.getcwd(), deltall_dir, each_qualified_cutoff)
        msa_dir         = '%s/%s/qualified_marker_PA_%s'                    % (os.getcwd(), assess_marker_pa_dir, each_qualified_cutoff)
        pwd_tree_file   = '%s/%s/PA_%s_C60_PMSF_tree/concatenated.treefile' % (os.getcwd(), iqtree_dir, each_qualified_cutoff)
        deltaLL_stdout  = 'PA_%s_stdout.txt'                                % (each_qualified_cutoff)
        deltall_cmd     = 'ruby %s --force --cpu %s -T %s --outdir %s --indir %s -t %s --outgrp_file %s --taxa %s > %s' % ('$deltaLL_rb', js_cpu_num, js_cpu_num, deltaLL_wd, msa_dir, pwd_tree_file, '$outgrp_file', '$taxa_list_txt', deltaLL_stdout)

        with open(pwd_js, 'w') as pwd_js_handle:
            pwd_js_handle.write('#!/bin/bash\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task %s\n\n' % js_cpu_num)
            pwd_js_handle.write('source activate ruby\n')
            pwd_js_handle.write('RUBYLIB=$RUBYLIB:/home-user/wzsong/Software/ruby_lib_sswang/\n')
            pwd_js_handle.write('export RUBYLIB\n')
            pwd_js_handle.write('cd %s/%s\n' % (os.getcwd(), deltall_dir))
            pwd_js_handle.write('outgrp_file="deltaLL_outgroup.txt"\n')
            pwd_js_handle.write('taxa_list_txt="deltaLL_eu_taxa_list.txt"\n')
            pwd_js_handle.write('deltaLL_rb="/home-user/wzsong/Scripts/deltaLL.rb"\n')
            pwd_js_handle.write(deltall_cmd + '\n')

    ##########################################################

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',               required=True,                          help='marker seq dir')
    parser.add_argument('-mx',              required=True,                          help='marker seq ext')
    parser.add_argument('-aa',              required=True,                          help='faa file dir')
    parser.add_argument('-aax',             required=True,                          help='faa file ext')
    parser.add_argument('-g',               required=True,                          help='genome group')
    parser.add_argument('-pac',             required=False, default='0-50-75-100',  help='cutoffs, default: 0-50-75-100')
    parser.add_argument('-o',               required=True,                          help='output dir')
    parser.add_argument('-e',               required=False, default='1e-30',        help='e-value cutoff, default: 1e-30')
    parser.add_argument('-t',               required=True,  type=int,               help='num of threads')
    parser.add_argument('-skip_align_trim', required=False, action="store_true",    help='skip extracting, aligning and trimming markers')
    parser.add_argument('-mmn',             required=False, default=10, type=int,   help='minimal marker number, default: 10')
    parser.add_argument('-jst',             required=False, default='6',            help='threads to request in job script, for running iqtree')
    parser.add_argument('-qsub',            required=False, action="store_true",    help='submit job scripts')
    parser.add_argument('-f',               required=False, action="store_true",    help='force overwrite')
    args = vars(parser.parse_args())
    MarkerRef2Tree(args)


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
