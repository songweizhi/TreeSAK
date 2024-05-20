from __future__ import print_function
import os
import glob
import argparse
from Bio import SeqIO


SplitScore1_usage = '''
======================== SplitScore1 example commands ========================

TreeSAK SplitScore1 -i marker_seq -x fa -o SplitScore1_op_dir -jst 9 -f

# Format of gene id
APA_bin56_00001
APA_bin56_00002
APA_bin56_00003

==============================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def SplitScore1(args):

    oma_op_fasta        = args['i']
    fasta_file_ext      = args['x']
    interested_gnm_txt  = args['u']
    iqtree_model        = args['m']
    cov_cutoff          = args['c']
    force_overwrite     = args['f']
    num_of_js_threads   = args['jst']
    op_dir              = args['o']

    ################################################################################

    interested_gnm_set = set()
    if interested_gnm_txt is not None:
        if os.path.isfile(interested_gnm_txt):
            for each_gnm in open(interested_gnm_txt):
                interested_gnm_set.add(each_gnm.strip())
        else:
            print('%s not found, program exited' % interested_gnm_txt)
            exit()

    ################################################################################

    fa_file_re   = '%s/*.%s' % (oma_op_fasta, fasta_file_ext)
    fa_file_list = glob.glob(fa_file_re)
    if len(fa_file_list) == 0:
        print('No file found in %s, progeam exited!' % oma_op_fasta)
        exit()

    og_to_gene_dict = dict()
    for each_fa in fa_file_list:
        _, f_base, _ = sep_path_basename_ext(each_fa)
        seq_id_set = set()
        for each_seq in SeqIO.parse(each_fa, 'fasta'):
            seq_id_set.add(each_seq.id)
        og_to_gene_dict[f_base] = seq_id_set

    ################################################################################

    gnm_to_process = set()
    for each_og in og_to_gene_dict:
        gene_set = og_to_gene_dict[each_og]
        gnm_set = set()
        for each_gene in gene_set:
            gnm_id = '_'.join(each_gene.split('_')[:-1])
            gnm_set.add(gnm_id)
            if interested_gnm_txt is None:
                gnm_to_process.add(gnm_id)
            else:
                if gnm_id in interested_gnm_set:
                    gnm_to_process.add(gnm_id)

        if len(gene_set) != len(gnm_set):
            print('Program exited!')
            exit()

    ################################################################################

    # define file name
    qualified_og_dir    = '%s/qualified_OGs'        % op_dir
    cmd_1_mafft_txt     = '%s/cmd_1_mafft.txt'      % op_dir
    cmd_2_trimal_txt    = '%s/cmd_2_trimal.txt'     % op_dir
    cmd_3_iqtree_txt    = '%s/cmd_3_iqtree.txt'     % op_dir
    ignored_marker_txt  = '%s/ignored_markers.txt'  % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)
    os.mkdir(qualified_og_dir)

    ################################################################################

    cmd_1_mafft_txt_handle = open(cmd_1_mafft_txt, 'w')
    cmd_2_trimal_txt_handle = open(cmd_2_trimal_txt, 'w')
    cmd_3_iqtree_txt_handle = open(cmd_3_iqtree_txt, 'w')
    ignored_og_dict = dict()
    for each_og in sorted(list(og_to_gene_dict.keys())):
        seq_file_in          = '%s/%s.%s'       % (oma_op_fasta, each_og, fasta_file_ext)
        file_out_seq         = '%s/%s.%s'       % (qualified_og_dir, each_og, fasta_file_ext)
        file_out_aln         = '%s.aln'         % each_og
        file_out_aln_trimmed = '%s_trimmed.aln' % each_og

        seq_file_out_handle = open(file_out_seq, 'w')
        current_gnm_set = set()
        for each_seq in SeqIO.parse(seq_file_in, 'fasta'):
            seq_id = each_seq.id
            gnm_id = '_'.join(seq_id.split('_')[:-1])
            if gnm_id in gnm_to_process:
                current_gnm_set.add(gnm_id)
                seq_file_out_handle.write('>%s\n' % each_seq.id)
                seq_file_out_handle.write('%s\n'  % each_seq.seq)
        seq_file_out_handle.close()

        cov_value = len(current_gnm_set)*100/len(gnm_to_process)
        cov_value = float("{0:.2f}".format(cov_value))

        if cov_value < cov_cutoff:
            report_str = 'Ignored %s, contains proteins from %s (%s%s) genomes, < %s%s.' % (each_og, len(current_gnm_set), cov_value, '%', cov_cutoff, '%')
            ignored_og_dict[each_og] = report_str
            os.system('rm %s' % file_out_seq)
        else:
            # align, trim and iqtree
            mafft_cmd     = 'mafft-einsi --thread %s --quiet %s.%s > %s'                                % (num_of_js_threads, each_og, fasta_file_ext, file_out_aln)
            trimal_cmd    = 'trimal -in %s -out %s -automated1'                                         % (file_out_aln, file_out_aln_trimmed)
            iqtree_cmd    = 'iqtree2 -s %s --seqtype AA -m %s -B 1000 --wbtl --bnni --prefix %s -T %s --quiet' % (file_out_aln_trimmed, iqtree_model, each_og, num_of_js_threads)
            # Undinarchaeota illuminate DPANN phylogeny and the impact of gene transfer on archaeal evolution, settings: -m LG+G -bb 1000 -wbtl -bnni
            cmd_1_mafft_txt_handle.write(mafft_cmd + '\n')
            cmd_2_trimal_txt_handle.write(trimal_cmd + '\n')
            cmd_3_iqtree_txt_handle.write(iqtree_cmd + '\n')
    cmd_1_mafft_txt_handle.close()
    cmd_2_trimal_txt_handle.close()
    cmd_3_iqtree_txt_handle.close()

    # report ignored markers
    if len(ignored_og_dict) > 0:
        print('The following %s markers were ignored due to low genome coverage, see details in %s:' % (len(ignored_og_dict), ignored_marker_txt))
        print('\n'.join(sorted(list(ignored_og_dict.keys()))))
        ignored_marker_txt_handle = open(ignored_marker_txt, 'w')
        for each_ignored_marker in sorted(list(ignored_og_dict.keys())):
            ignored_marker_txt_handle.write(ignored_og_dict[each_ignored_marker] + '\n')
        ignored_marker_txt_handle.close()

    # report
    print('You will need to execute the commands exported to the following three files before moving to SplitScore2')
    print(cmd_1_mafft_txt)
    print(cmd_2_trimal_txt)
    print(cmd_3_iqtree_txt)
    print('Done!')


if __name__ == '__main__':

    SplitScore1_parser = argparse.ArgumentParser()
    SplitScore1_parser.add_argument('-i',   required=True,                          help='orthologous gene sequence')
    SplitScore1_parser.add_argument('-x',   required=True,                          help='fasta file extension')
    SplitScore1_parser.add_argument('-o',   required=True,                          help='output directory')
    SplitScore1_parser.add_argument('-u',   required=False, default=None,           help='interested genomes, no file extension')
    SplitScore1_parser.add_argument('-m',   required=False, default='LG+G+I',       help='iqtree_model, default: LG+G+I')
    SplitScore1_parser.add_argument('-c',   required=False, type=int, default=85,   help='coverage cutoff, default: 85')
    SplitScore1_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    SplitScore1_parser.add_argument('-jst', required=False, type=int, default=1,    help='num of threads for iqtree2, default: 1')
    args = vars(SplitScore1_parser.parse_args())
    SplitScore1(args)
