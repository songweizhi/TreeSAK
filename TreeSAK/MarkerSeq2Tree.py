import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO
from distutils.spawn import find_executable


MarkerSeq2Tree_usage = '''
======================== MarkerSeq2Tree example commands ========================

Dependencies: mafft, trimal, bmge, perl and iqtree2

TreeSAK MarkerSeq2Tree -i best_25 -x fa -o op_dir -t 12 -f -bmge -prune 10,20,30

=================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


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


def pruneMSA(msa_in, conserved_cutoffs):

    msa_path, msa_base, msa_ext = sep_path_basename_ext(msa_in)

    current_file_path   = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    alignment_pruner_pl = '%s/alignment_pruner.pl'  % current_file_path
    cutoff_list         = conserved_cutoffs.split(',')

    op_file_list = []
    for each_cutoff in cutoff_list:
        cutoff_formatted = str(float(each_cutoff)/100).replace('0.', '.')
        current_msa_out = '%s/%s.pruner%s%s' % (msa_path, msa_base, each_cutoff, msa_ext)
        perl_cmd            = 'perl %s --file %s --conserved_threshold %s > %s' % (alignment_pruner_pl,   msa_in, cutoff_formatted, current_msa_out)
        perl_cmd_for_report = 'perl %s --file %s --conserved_threshold %s > %s' % ('alignment_pruner.pl', msa_in, cutoff_formatted, current_msa_out)
        op_file_list.append(current_msa_out)
        print(perl_cmd_for_report)
        os.system(perl_cmd)

    # report
    print('Pruned MSA exported to:')
    print('\n'.join(op_file_list))


def MarkerSeq2Tree(args):

    marker_seq_dir              = args['i']
    marker_seq_ext              = args['x']
    op_dir                      = args['o']
    num_of_threads              = args['t']
    run_bmge                    = args['bmge']
    trim_with_bmge              = args['bmge']
    bmge_trim_model             = args['bmge_m']
    bmge_entropy_score_cutoff   = args['bmge_esc']
    force_overwrite             = args['f']
    alignment_pruner_cutoffs    = args['prune']

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar      = '%s/BMGE.jar' % current_file_path

    # check dependencies
    not_detected_programs = []
    for needed_program in ['mafft-einsi', 'trimal', 'iqtree2']:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)
    if not_detected_programs != []:
        print('%s not detected, program exited!' % ', '.join(not_detected_programs))
        exit()

    # get marker id set
    marker_seq_re   = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
    marker_seq_list = sorted(glob.glob(marker_seq_re))

    # define output dir
    renamed_marker_seq_dir              = '%s/renamed_markers'                      % op_dir
    renamed_marker_aln_dir              = '%s/renamed_markers_aln'                  % op_dir
    if trim_with_bmge is False:
        cmds_2_trim_txt                 = '%s/cmds_2_trimal.txt'                    % op_dir
        renamed_marker_aln_dir_trimmed  = '%s/renamed_markers_aln_trimal'           % op_dir
    else:
        cmds_2_trim_txt                 = '%s/cmds_2_BMGE.txt'                      % op_dir
        renamed_marker_aln_dir_trimmed  = '%s/renamed_markers_aln_BMGE'             % op_dir
    concatenated_phy                    = '%s/concatenated.phy'                     % op_dir
    concatenated_phy_fasta              = '%s/concatenated.phy.fasta'               % op_dir
    concatenated_phy_partition          = '%s/concatenated_partition.txt'           % op_dir
    iqtree_dir                          = '%s/iqtree_wd'                            % op_dir
    cmds_1_mafft_txt                    = '%s/cmds_1_mafft.txt'                     % op_dir
    cmds_3_iqtree_txt                   = '%s/cmds_3_iqtree2.txt'                   % op_dir
    pwd_guide_tree                      = '%s/iqtree_wd/guide_tree.treefile'        % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)
    os.mkdir(renamed_marker_seq_dir)
    os.mkdir(renamed_marker_aln_dir)
    os.mkdir(renamed_marker_aln_dir_trimmed)

    # write out best hits and extract sequences
    for marker_seq_file in marker_seq_list:

        f_path, f_base, f_ext = sep_path_basename_ext(marker_seq_file)
        pwd_renamed_marker_seq          = '%s/%s.%s'  % (renamed_marker_seq_dir, f_base, marker_seq_ext)
        pwd_renamed_marker_aln          = '%s/%s.aln' % (renamed_marker_aln_dir, f_base)
        pwd_renamed_marker_aln_trimmed  = '%s/%s.aln' % (renamed_marker_aln_dir_trimmed, f_base)

        # rename sequences
        marker_hits_seq_renamed_handle = open(pwd_renamed_marker_seq, 'w')
        for each_seq in SeqIO.parse(marker_seq_file, 'fasta'):
            seq_id = each_seq.id
            seq_gnm = '_'.join(seq_id.split('_')[:-1])
            marker_hits_seq_renamed_handle.write('>%s\n' % seq_gnm)
            marker_hits_seq_renamed_handle.write('%s\n' % str(each_seq.seq))
        marker_hits_seq_renamed_handle.close()

        # align
        mafft_cmd  = 'mafft-einsi --thread %s --quiet %s > %s'          % (num_of_threads, pwd_renamed_marker_seq, pwd_renamed_marker_aln)

        # trim
        trim_cmd = 'trimal -in %s -out %s -automated1'                  % (pwd_renamed_marker_aln, pwd_renamed_marker_aln_trimmed)
        if trim_with_bmge is True:
            trim_cmd = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'    % (pwd_bmge_jar, pwd_renamed_marker_aln, bmge_trim_model, bmge_entropy_score_cutoff, pwd_renamed_marker_aln_trimmed)

        # write out mafft cmds
        with open(cmds_1_mafft_txt, 'a') as cmds_1_mafft_txt_handle:
            cmds_1_mafft_txt_handle.write(mafft_cmd + '\n')

        # write out trimal cmds
        with open(cmds_2_trim_txt, 'a') as cmds_2_trim_txt_handle:
            cmds_2_trim_txt_handle.write(trim_cmd + '\n')

        # run cmds
        os.system(mafft_cmd)
        os.system(trim_cmd)

    # concatenate alignments
    catfasta2phy(renamed_marker_aln_dir_trimmed, 'aln', concatenated_phy, concatenated_phy_partition)

    # run iqtree2
    os.mkdir(iqtree_dir)
    get_guide_tree_cmd  = 'iqtree2 --seqtype AA -T %s -B 1000 --alrt 1000 --quiet -s %s --prefix %s/guide_tree -m LG '                  % (num_of_threads, concatenated_phy, iqtree_dir, )
    get_c60_tree_cmd    = 'iqtree2 --seqtype AA -T %s -B 1000 --alrt 1000 --quiet -s %s --prefix %s/concatenated -m LG+C60+G+F -ft %s'  % (num_of_threads, concatenated_phy, iqtree_dir, pwd_guide_tree)

    # write out iqtree2 cmds
    with open(cmds_3_iqtree_txt, 'a') as cmds_3_iqtree_txt_handle:
        cmds_3_iqtree_txt_handle.write(get_guide_tree_cmd + '\n')
        cmds_3_iqtree_txt_handle.write(get_c60_tree_cmd + '\n')

    # run alignment_pruner.pl
    if alignment_pruner_cutoffs is not None:
        pruneMSA(concatenated_phy_fasta, alignment_pruner_cutoffs)

    # run cmds
    print('Running iqtree')
    # os.system(get_guide_tree_cmd)
    # os.system(get_c60_tree_cmd)

    print('You may want to submit the following commands to infer tree')
    print('To be added...')
    print('Done!')


if __name__ == '__main__':

    MarkerSeq2Tree_parser = argparse.ArgumentParser()
    MarkerSeq2Tree_parser.add_argument('-i',           required=True,                          help='marker seq dir')
    MarkerSeq2Tree_parser.add_argument('-x',           required=True,                          help='marker seq ext')
    MarkerSeq2Tree_parser.add_argument('-o',           required=True,                          help='output dir')
    MarkerSeq2Tree_parser.add_argument('-t',           required=False, type=int, default=1,    help='num of threads')
    MarkerSeq2Tree_parser.add_argument('-bmge',        required=False, action="store_true",    help='perform BMGE trimming on concatenated MSA')
    MarkerSeq2Tree_parser.add_argument('-bmge_m',      required=False, default='BLOSUM30',     help='BMGE trim model, default: BLOSUM30')
    MarkerSeq2Tree_parser.add_argument('-bmge_esc',    required=False, default='0.55',         help='BMGE entropy score cutoff, default: 0.55')
    MarkerSeq2Tree_parser.add_argument('-prune',       required=False, default=None,           help='conservation cutoffs for alignment_pruner.pl')
    MarkerSeq2Tree_parser.add_argument('-f',           required=False, action="store_true",    help='force overwrite')
    args = vars(MarkerSeq2Tree_parser.parse_args())
    MarkerSeq2Tree(args)
