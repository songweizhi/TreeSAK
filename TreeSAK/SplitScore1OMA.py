from __future__ import print_function
import os
import argparse
from Bio import SeqIO


SplitScore1OMA_usage = '''
======================== SplitScore1OMA example commands ========================

# SplitScore1
TreeSAK SplitScore1OMA -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f
TreeSAK SplitScore1OMA -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f -u interested_gnm.txt
# Please ensure that all the commands in iqtree_cmds.txt have been executed before proceeding to step 2.

=================================================================================
'''


def select_seq(seq_file, seq_id_list, output_file):
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        seq_id = seq_record.id
        if seq_id in seq_id_list:
            output_file_handle.write('>%s\n' % seq_id)
            output_file_handle.write('%s\n' % str(seq_record.seq))
    output_file_handle.close()


def get_gene_tree(oma_op_txt, oma_op_fasta, interested_gnm_txt, cov_cutoff, oma_op_fasta_qualified, iqtree_model, num_of_js_threads, force_overwrite, get_gene_tree_cmd_txt):

    # get the total number of genome
    genome_id_set = set()
    for each_group in open(oma_op_txt):
        if not each_group.startswith('#'):
            for each_gene in each_group.strip().split('\t')[1:]:
                gnm_id = '_'.join(each_gene.split(':')[1].split(' ')[0].split('_')[:-1])
                genome_id_set.add(gnm_id)

    interested_gnm_set = set()
    if interested_gnm_txt is not None:
        for each_gnm in open(interested_gnm_txt):
            interested_gnm_set.add(each_gnm.strip())
    else:
        interested_gnm_set = genome_id_set

    # create output folder
    if os.path.isdir(oma_op_fasta_qualified) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % oma_op_fasta_qualified)
        else:
            print('%s already exist, program exited!' % oma_op_fasta_qualified)
            exit()
    os.system('mkdir %s' % oma_op_fasta_qualified)

    # filter OMA output
    qualified_grp_to_gene_dict = dict()
    for each_group in open(oma_op_txt):
        if not each_group.startswith('#'):
            each_group_split  = each_group.strip().split('\t')
            group_id          = each_group_split[0]
            gene_list_by_gnm  = each_group_split[1:]
            current_gene_list = [i.split(':')[1].split(' ')[0] for i in gene_list_by_gnm]
            current_gnm_list_interested = []
            current_gene_list_interested = []
            for gene in current_gene_list:
                gnm = '_'.join(gene.split('_')[:-1])
                if gnm in interested_gnm_set:
                    current_gnm_list_interested.append(gnm)
                    current_gene_list_interested.append(gene)

            current_cov = len(current_gnm_list_interested) * 100 / len(interested_gnm_set)
            if current_cov >= cov_cutoff:
                qualified_grp_to_gene_dict[group_id] = current_gene_list_interested

    print('The number of orthologous groups with coverage >= %s is %s.' % (cov_cutoff, len(qualified_grp_to_gene_dict)))

    # prepare commands for getting gene tree
    get_gene_tree_cmd_txt_handle = open(get_gene_tree_cmd_txt, 'w')
    for qualified_grp in sorted(list(qualified_grp_to_gene_dict.keys())):
        group_id_only_num = qualified_grp.replace('OMA', '')
        while group_id_only_num[0] == '0':
            group_id_only_num = group_id_only_num[1:]

        # define file name
        og_id               = 'OG%s'                % group_id_only_num
        pwd_seq_file_in     = '%s/%s.fa'            % (oma_op_fasta, og_id)
        pwd_og_seq          = '%s/%s.fa'            % (oma_op_fasta_qualified, og_id)
        pwd_og_aln          = '%s/%s.aln'           % (oma_op_fasta_qualified, og_id)
        pwd_og_aln_trimmed  = '%s/%s_trimmed.aln'   % (oma_op_fasta_qualified, og_id)

        # get sequence
        if len(interested_gnm_set) == len(genome_id_set):
            cp_cmd = 'cp %s %s' % (pwd_seq_file_in, pwd_og_seq)
            os.system(cp_cmd)
        else:
            select_seq(pwd_seq_file_in, qualified_grp_to_gene_dict[qualified_grp], pwd_og_seq)

        # align, trim and iqtree
        mafft_cmd     = 'mafft-einsi --thread %s --quiet %s > %s'                                      % (num_of_js_threads, pwd_og_seq, pwd_og_aln)
        trimal_cmd    = 'trimal -in %s -out %s -automated1'                                            % (pwd_og_aln, pwd_og_aln_trimmed)
        iqtree_cmd    = 'iqtree2 -s %s --seqtype AA -m %s -T %s -B 1000 --quiet --wbtl --prefix %s/%s' % (pwd_og_aln_trimmed, iqtree_model, num_of_js_threads, oma_op_fasta_qualified, og_id)
        cmds_one_line = '%s; %s; %s'                                                                   % (mafft_cmd, trimal_cmd, iqtree_cmd)
        get_gene_tree_cmd_txt_handle.write(cmds_one_line.replace((oma_op_fasta_qualified + '/'), '') + '\n')
    get_gene_tree_cmd_txt_handle.close()


def SplitScore1OMA(args):

    oma_op_txt              = args['i']
    oma_op_fasta            = args['s']
    interested_gnm_txt      = args['u']
    iqtree_model            = args['m']
    cov_cutoff              = args['c']
    force_overwrite         = args['f']
    num_of_js_threads       = args['jst']
    step_1_op_dir           = args['o']

    # define file name
    qualified_og_dir    = '%s/qualified_OGs'    % step_1_op_dir
    iqtree_cmds_txt     = '%s/iqtree_cmds.txt'  % step_1_op_dir

    # create output folder
    if os.path.isdir(step_1_op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % step_1_op_dir)
        else:
            print('%s exist, program exited!' % step_1_op_dir)
            exit()
    os.mkdir(step_1_op_dir)
    os.mkdir(qualified_og_dir)

    # get get_gene_tree
    get_gene_tree(oma_op_txt, oma_op_fasta, interested_gnm_txt, cov_cutoff, qualified_og_dir, iqtree_model, num_of_js_threads, force_overwrite, iqtree_cmds_txt)


if __name__ == '__main__':

    SplitScore1OMA_parser = argparse.ArgumentParser()
    SplitScore1OMA_parser.add_argument('-i',   required=True,                        help='OrthologousGroups.txt, produced by OMA')
    SplitScore1OMA_parser.add_argument('-s',   required=True,                        help='OrthologousGroupsFasta, produced by OMA')
    SplitScore1OMA_parser.add_argument('-u',   required=False, default= None,        help='ID of interested genomes, no file extension')
    SplitScore1OMA_parser.add_argument('-o',   required=True,                        help='output directory')
    SplitScore1OMA_parser.add_argument('-m',   required=False, default='LG+G+I',     help='iqtree_model, default: LG+G+I')
    SplitScore1OMA_parser.add_argument('-c',   required=False, type=int, default=80, help='coverage cutoff, default: 80')
    SplitScore1OMA_parser.add_argument('-f',   required=False, action="store_true",  help='force overwrite')
    SplitScore1OMA_parser.add_argument('-jst', required=False, type=int, default=1,  help='num of threads for inferring gene tree, default: 1')
    args = vars(SplitScore1OMA_parser.parse_args())
    SplitScore1OMA(args)
